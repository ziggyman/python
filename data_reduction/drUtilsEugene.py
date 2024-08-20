import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from PyAstronomy import pyasl
from scipy.optimize import curve_fit
from scipy import exp
import os




def getImageData(fname,hduNum=1):
    hdulist = pyfits.open(fname)
    scidata = hdulist[hduNum].data
    hdulist.close()
    return np.array(scidata)

def getHeader(fName, hduNum=0):
    hdulist = pyfits.open(fName)
    header = hdulist[hduNum].header
    hdulist.close()
    return header

def getHeaderValue(fname, keyword, hduNum=0):
    hdulist = pyfits.open(fname)
    header = hdulist[hduNum].header
    hdulist.close()
    try:
        return header[keyword]
    except:
        return None

def getWavelengthArr(fname, hduNum=0):
    hdulist = pyfits.open(fname)
    header = hdulist[hduNum].header
    hdulist.close()
    if 'CDELT1' in header.keys():
        cdelt = header['CDELT1']
        wLen = ((np.arange(header['NAXIS1']) + 1.0) - header['CRPIX1']) * cdelt + header['CRVAL1']
    elif 'CD1_1' in header.keys():
        cdelt = header['CD1_1']
        wLen = ((np.arange(header['NAXIS1']) + 1.0) - header['CRPIX1']) * cdelt + header['CRVAL1']
    else:
        print('WARNING: neither CDELT1 nor CD1_1 found in header of file <'+fname+'>')
        wLen = np.arange(header['NAXIS1']) + 1.0
    return wLen

def gauss(x,a,x0,sigma,yBackground=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground

def gaussNorm(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))+1.

#def gauss_linOnly(x,a,x0,sigma,linear=0.):
#    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground+(linear*x)

def gauss_lin(x,a,x0,sigma,yBackground=0.,linear=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground+(linear*x)

def normalizeX(x):
    xZero = x - x[0]
    return 2.0 * xZero / xZero[-1] - 1.0

# @brief : fit background and subtract from y
# @param x : 1D array of x-values
# @param y : 1D array of y-values
# @param deg : int: degree for Legendre Polynomial
# @param indicesToIgnore : 1D array of indices to ignore for fitting the background
def sigmaRej(values,
             sigLow,
             sigHigh,
             replace=False,
             adjustSigLevels=False,
             useMean=False,
             keepFirstAndLastX=False):
    print('sigmaRej: sigLow = ',sigLow,', sigHigh = ',sigHigh,', adjustSigLevels = ',adjustSigLevels,', replace = ',replace,', useMean = ',useMean,', keepFirstAndLastX = ',keepFirstAndLastX)
    ySkyMedian = None
    if useMean:
        ySkyMedian = np.mean(values)
    else:
        ySkyMedian = np.median(values)
    sigma = np.std(values)
    nRej = 0
    outArr = []
    outIndices = []
    if adjustSigLevels:
        if sigma > 3. * ySkyMedian:
            sigLow = 0.5
            sigHigh = 0.05
    for i in range(len(values)):
        if (values[i] < ySkyMedian - (sigLow * sigma)) or (values[i] > ySkyMedian + (sigHigh * sigma)):
            if ((i == 0) and keepFirstAndLastX):
                outArr.append(np.mean(values[i:i+3]))
                outIndices.append(i)
            elif ((i == len(values)-1) and keepFirstAndLastX):
                outArr.append(np.mean(values[i-2:i+1]))
                outIndices.append(i)
            else:
                if replace:
                    outArr.append(ySkyMedian)
                    outIndices.append(i)
                nRej += 1
                print('sigmaRej: median = ',ySkyMedian,', sigma = ',sigma,': removed pixel ',i,' from sky with value ',values[i])
                if (values[i] < ySkyMedian - (sigLow * sigma)):
                    print('sigmaRej: values[',i,'](=',values[i],') < ySkyMedian(=',ySkyMedian,') - (sigLow(=',sigLow,') * sigma(=',sigma,'))(=',sigLow*sigma,')=',ySkyMedian - (sigLow * sigma))
                else:
                    print('sigmaRej: values[',i,'](=',values[i],') > ySkyMedian(=',ySkyMedian,') + (sigHigh(=',sigHigh,') * sigma(=',sigma,'))(=',sigHigh*sigma,')=',ySkyMedian + (sigHigh * sigma))
        else:
            outArr.append(values[i])
            outIndices.append(i)
    if nRej > 0:
        print('sigmaRej: rejected ',nRej,' out of ',len(values),' pixels')
    return [np.array(outArr),np.array(outIndices)]

def sigmaReject(y,
                nIter,
                lowReject,
                highReject,
                replace=False,
                adjustSigLevels=False,
                useMean=False,
                keepFirstAndLastX=True):
    indices = np.arange(0,len(y),1)
    yGood, goodIndices = sigmaRej(y, lowReject, highReject, replace=replace, adjustSigLevels=adjustSigLevels, useMean=useMean, keepFirstAndLastX=keepFirstAndLastX)
    indices = indices[goodIndices]
    #print('sigmaReject: iter = 0: indices = ',indices.shape,': ',indices)
    for iter in np.arange(1,nIter,1):
        yGood, goodIndices = sigmaRej(yGood, lowReject, highReject, replace=replace, adjustSigLevels=adjustSigLevels, useMean=useMean, keepFirstAndLastX=keepFirstAndLastX)
        indices = indices[goodIndices]
        #print('sigmaReject: iter = ',iter,': yGood = ',yGood.shape,': ',yGood)
        #print('sigmaReject: iter = ',iter,': indices = ',indices.shape,': ',indices)
    return [yGood,indices]

# NOTE: requires x to be normalized
def sfit(x,
         y,
         fittingFunction,
         solveFunction,
         order,
         nIterReject,
         nIterFit,
         lowReject,
         highReject,
         adjustSigLevels=False,
         useMean=False,
         display=False):
    coeffs = fittingFunction(x,y,order)
    fittedValues = np.array(solveFunction(x,coeffs))
    fittedIndices = np.arange(0,fittedValues.shape[0],1)
    print('fittedValues = ',fittedValues)
    print('fittedIndices = ',fittedIndices)
    print('sfit: iter = -1: x.shape = ',x.shape,', y.shape = ',y.shape,', fittedValues.shape = ',fittedValues.shape,', fittedIndices = ',fittedIndices.shape,': ',fittedIndices)
    fittedValuesTemp = fittedValues
    for i in range(nIterFit):
        #print('sfit: iter = ',i,': x.shape = ',x.shape,', y.shape = ',y.shape,', fittedValues.shape = ',fittedValues.shape,', fittedIndices = ',fittedIndices.shape,': ',fittedIndices)
        #print('sfit: y[fittedIndices] = ',y[fittedIndices])
        if display:
            plt.plot(x,fittedValues,label='previous fit')
            plt.plot(x[fittedIndices],y[fittedIndices], label='input')
            plt.plot(x[fittedIndices],y[fittedIndices] - fittedValuesTemp, label='difference')
        fittedValuesNotRejected, fittedIndicesTemp = sigmaReject(y[fittedIndices] - fittedValuesTemp,
                                                                 nIter=nIterReject,
                                                                 lowReject=lowReject,
                                                                 highReject=highReject,
                                                                 replace=False,
                                                                 adjustSigLevels=adjustSigLevels,
                                                                 useMean=useMean,
                                                                 keepFirstAndLastX=False,
                                                                )
        fittedValuesNotRejected += fittedValuesTemp[fittedIndicesTemp]
        fittedIndices = fittedIndices[fittedIndicesTemp]
        #print('sfit: y[fittedIndices] = ',y[fittedIndices].shape,': ',y[fittedIndices])
        #print('sfit: fittedValuesNotRejected = ',fittedValuesNotRejected.shape,': ',fittedValuesNotRejected)
        """x is already required to be normalized at function call"""
        coeffs = fittingFunction(x[fittedIndices],y[fittedIndices],order)
        """keep first and last x values to get the normalization right!!!"""
        fittedValues = solveFunction(x,coeffs)
        #print('sfit: sfit: iter = ',i,': x.shape = ',x.shape,', y.shape = ',y.shape,', fittedValues.shape = ',fittedValues.shape,', fittedIndices = ',fittedIndices.shape,': ',fittedIndices)
        fittedValuesTemp = fittedValues[fittedIndices]
        if display:
            #print('sfit: sfit: iter = ',i,': y[fittedIndices] = ',y[fittedIndices],', fittedValuesNotRejected = ',fittedValuesNotRejected)
            plt.plot(x[fittedIndices],y[fittedIndices] - fittedValuesTemp,'b*', label='fitted points')
            plt.plot(x[fittedIndices],y[fittedIndices],'r+', label='fitted points')
            plt.plot(x,fittedValues,label='new fit')
            plt.legend()
            plt.show()
    return [coeffs, fittedValues]

# NOTE: only takes the median of the lower half of the sorted values in order to allow the sky subtraction for extended objects
def writeFits1D(flux, outFileName, wavelength=None, header=None, CRVAL1=None, CRPIX1=None, CDELT1=None):
    head = None
    if not header is None:
        if isinstance(header,str):
#            refFileName = header
            hdulist = pyfits.open(header)
            head = hdulist[0].header
            hdulist.close()
        else:
            head = header
        if 'NAXIS2' in head.keys():
            del head['NAXIS2']
        print('writeFits1D: dir(head) = ',dir(head))
        head['NAXIS1'] = flux.shape[0]
        for key in head:
            print('writeFits1D: ',key+': <',head[key],'>: ',type(head[key]))
            if isinstance(head[key],str):
                if '\n' in head[key]:
                    head[key].replace('\n','')
            elif (key == 'COMMENT'):# and (head[key] == ''):
                del head[key]
            elif '\n' in str(head[key]):
                del head[key]
#            elif isinstance(head[key],'astropy.io.fits.header._HeaderCommentaryCards'):
#                del head[key]


    waveParams = None
    if not CRVAL1 is None:
        waveParams = {'CRVAL1': CRVAL1,
                      'CRPIX1': CRPIX1,
                      'CDELT1': CDELT1}
    #print('writeFits1D: flux = ',flux)
    #print('writeFits1D: wavelength = ',wavelength)
    #print('writeFits1D: waveParams = ',waveParams)
    #print('writeFits1D: head = ',head)
    print('writeFits1D: writing file <'+outFileName+'>')
    pyasl.write1dFitsSpec(outFileName, flux, wvl=wavelength, waveParams=waveParams, fluxErr=None, header=head, clobber=True, refFileName=None, refFileExt=0)


def markEmissionLines(xRange, yRange):
    lines = [['[OII]',3727.],
            ['[NeIII]',3969.],
            ['HeI',4026.],
            ['[SII]',4072.],
            ['Hδ',4102.],
            ['CII',4267.],
            ['Hγ',4340.],
            ['[OIII]',4363.],
            ['HeI',4388.],
            ['HeI',4472.],
            ['HeII',4542.],
            ['[MgI]',4571.],
            ['HeII',4686.],
            ['[ArIV]',4740.],
            ['Hβ',4861.],
            ['HeI',4922.],
            ['[OIII]',4959.],
            ['[OIII]',5007.],
            ['[NI]',5199.],
            ['HeII',5412.],
            ['[ClIII]',5518.],
            ['[ClIII]',5538.],
            ['[OI]',5577.],
            ['[NII]',5754.],
            ['HeI',5876.],
            ['[OI]',6300.],
            ['[SIII]',6312.],
            ['[OI]',6364.],
            ['[ArV]',6435.],
            ['[NII]',6548.],
            ['Hα',6563.],
            ['[NII]',6583.],
            ['HeI',6678.],
            ['[SII]',6716.],
            ['[SII]',6731.],
            ['HeII',6891.],
            ['[ArV]',7006.],
            ['HeI',7065.],
            ['[ArIII]',7136.],
            ['HeII',7176.],
            ['[ArIV]',7237.],
            ['[ArIV]',7263.],
            ['HeI',7281.],
            ['[OII]',7325.],
            ['[SIII]',9069.],
            ['[SIII]',9532.],
            ]

#    if xRange is None:
#        xRange = [0,10000]
    for line in lines:
        if (line[1] >= xRange[0]) and (line[1] <= xRange[1]):
            plt.plot([line[1],line[1]],[yRange[0],yRange[1]])
            plt.text(line[1],yRange[1]+((yRange[1]-yRange[0])/50.),line[0],rotation='vertical')
    plt.xlim(xRange[0],xRange[1])
    plt.ylim(yRange[0],yRange[1])


#@params:
# spectrumFileNameIn: string (name of fits file to read) / data array
# spectrumFileNameOut: string
# fittingFunction: np.polynomial.legendre.legfit - function to fit continuum
# evalFunction: np.polynomial.legendre.legval - function to evaluate continuum
# order: fitting function polynomial order
# nIterReject: number of rejectsion iterations
# nIterFit: number of fitting iterations (how many times to attempt to reject bad values and fit the continuum)
# lowReject: reject values < fit - (lowReject * sigma)
# highReject: reject values > fit + (highReject * sigma)
# wLen: if set take wavelength vector from this, otherwise from fits header. Set this if spectrumFileNameIn is a data vector
# xLim: only fit from xLim[0] to xLim[1]
# regions: ignore values outside the regions (in wavelength)
# type: "difference" - subtract fitted continuum from spectrum
#       "ratio"      - divide spectrum by fitted continuum
#       "fit"        - write fitted continuum to spectrumFileNameOut
# adjustSigLevels: try to estimate lowReject and highReject parameters based on spectrum
# useMean: use mean instead of median for sigma rejection
# display: plot intermediate and final results
# returnCoeffs: return fitting function coefficients
def continuum(spectrumFileNameIn,
            spectrumFileNameOut,
            fittingFunction = np.polynomial.legendre.legfit,
            evalFunction = np.polynomial.legendre.legval,
            order = 5,
            nIterReject = 2,
            nIterFit = 2,
            lowReject = 2.,
            highReject = 3.5,
            wLen = None,
            xLim=None,#limits are both included
            regions=None,
            type='difference',
            adjustSigLevels=False,
            useMean=False,
            display = False,
            returnCoeffs = False):
    if isinstance(spectrumFileNameIn,str):
        print('continuum: reading file '+spectrumFileNameIn)
        specOrig = getImageData(spectrumFileNameIn,0)
    else:
        specOrig = spectrumFileNameIn
    if wLen is None:
        wLen = getWavelengthArr(spectrumFileNameIn,0)
    if display:
        plt.plot(wLen, specOrig,label='original spectrum')
#        plt.title(spectrumFileNameIn[spectrumFileNameIn.rfind('/')+1:spectrumFileNameIn.rfind('.')])
        plt.legend()
        plt.show()

    if xLim is None:
        xLim = [wLen[0],wLen[len(wLen)-1]]

    xrange = np.arange(0,specOrig.shape[0],1)
    specFit = specOrig
    wLenFit = wLen
    if xLim is not None:
        idx = np.where((wLen >= xLim[0]) & (wLen <= xLim[1]))[0]
        xrange = idx
        specFit = specOrig[idx]
        wLenFit = wLen[idx]
    if regions is not None:
        wave = []
        spectrum = []
        xrange = []
        #print('cotinuum: wLenFit = ',wLenFit)
        #print('cotinuum: wLenFit.shape = ',wLenFit.shape)
        minX = np.min(np.array([region[0] for region in regions]))
        maxX = np.max(np.array([region[1] for region in regions]))
        print('minX = ',minX)
        print('maxX = ',maxX)
        for i in range(wLenFit.shape[0]):
            for region in regions:
                if (wLenFit[i] >= xLim[0]) & (wLenFit[i] <= xLim[1]):
                    if (wLenFit[i] >= region[0]) & (wLenFit[i] <= region[1]):
                        wave.append(wLenFit[i])
                        spectrum.append(specFit[i])
            if (wLenFit[i] >= minX) and (wLenFit[i] <= maxX):
                xrange.append(i)
        print('wave = ',wave)
        print('spectrum = ',spectrum)
        print('xrange = ',xrange)
    else:
        wave = wLen
        spectrum = specOrig

#    if wLen is not None:
#        range = wLen
    #print('continuum: wave = ',len(wave),': ',wave)
    #print('continuum: spectrum = ',len(spectrum),': ',spectrum)

#    xToNorm = [xLim[0]]
#    for wav in wave:
#        xToNorm.append(wav)
#    xToNorm.append(xLim[1])
    #print('continuum: xToNorm = ',len(xToNorm),': ',xToNorm)
    xNorm = normalizeX(np.array(wave))
#    xNorm = xNorm[1:]
#    xNorm = xNorm[:len(xNorm)-1]
    print('continuum: xNorm = ',len(xNorm),': ',xNorm)
    print('continuum: spectrum = ',len(spectrum),': ',spectrum)
    #print('continuum: spectrum = ',len(spectrum),': ',spectrum)
    #STOP
    #sfit(x, y, fittingFunction, solveFunction, order, nIterReject, nIterFit, lowReject, highReject, adjustSigLevels=False, useMean=False, display=False)
    coeffs, resultFit = sfit(np.asarray(xNorm),
                             np.asarray(spectrum),
                             fittingFunction,#np.polynomial.legendre.legfit,
                             evalFunction,#np.polynomial.legendre.legval,
                             order=order,
                             nIterReject=nIterReject,
                             nIterFit=nIterFit,
                             lowReject=lowReject,
                             highReject=highReject,
                             adjustSigLevels=adjustSigLevels,
                             useMean=useMean,
                             display=display)
    print('continuum: resultFit = ',resultFit)
    fittedSpectrum = evalFunction(normalizeX(np.array(xrange)),coeffs)
    print('fitted continuum = ',fittedSpectrum)
    specDone = specOrig
    if type == 'difference':
        print('xrange = ',xrange)
        specDone[xrange] = specOrig[xrange] - fittedSpectrum
    elif type == 'ratio':
        specDone[xrange] = specOrig[xrange] / fittedSpectrum
    elif type == 'fit':
        specDone[xrange] = fittedSpectrum
    else:
        print('continuum: ERROR: type <'+type+'> not recognised')
        STOP

    if display:
        print('len(wLen) = ',len(wLen),', len(specOrig) = ',len(specOrig))
        plt.plot(wLen, specOrig, label='original')
        print('len(wLen) = ',len(wLen),', len(specDone) = ',len(specDone))
        x = np.arange(xLim[0],xLim[1]+1,1)
        fit = evalFunction(normalizeX(x),coeffs)
        plt.plot(x,fit,label='fit')
        plt.plot(wLen, specDone, label = 'continuum corrected')
        plt.legend()
        xRange = [wLen[0],wLen[len(wLen)-1]]
        yRange = [np.min([np.min(specOrig), np.min(specDone)]),np.max([np.max(specOrig), np.max(specDone)])]
        if (spectrumFileNameOut is not None) and (isinstance(spectrumFileNameOut,str)):
            plt.text(xRange[1],yRange[0],spectrumFileNameOut[spectrumFileNameOut.rfind('/')+1:spectrumFileNameOut.rfind('.')],rotation='vertical')
        markEmissionLines(xRange, yRange)
        plt.show()

    if (spectrumFileNameOut is not None) and isinstance(spectrumFileNameOut,str):
        writeFits1D(specDone,
                    spectrumFileNameOut,
                    wavelength=None,
                    header=spectrumFileNameIn,
                    CRVAL1=getHeaderValue(spectrumFileNameIn,'CRVAL1'),
                    CRPIX1=getHeaderValue(spectrumFileNameIn,'CRPIX1'),
                    CDELT1=getHeaderValue(spectrumFileNameIn,'CDELT1'),
                )
    elif spectrumFileNameOut is not None:
        spectrumFileNameOut = specDone
    for i in range(len(specOrig)-1):
        print('specOrig[',i,'] = ',specOrig[i],', specDone[',i,'] = ',specDone[i])
    if returnCoeffs:
        return coeffs
    return specDone

def fitLines(inputSpec1D,outputSpec):
    from matplotlib.widgets import RadioButtons, TextBox, Button
    from matplotlib.patches import Polygon

    global wLen
    global spec
    global xRange
    global continuum_order
    global continuum_high_reject
    global continuum_low_reject
    global spec_bak
    global wlen_bak
    global cleanType
    global newRegion
    global normRegion
    global regions
    spec_bak = None
    wlen_bak = None

    def submit_continuum_order(text):
        global continuum_order
        continuum_order = int(text)

    def submit_continuum_low_reject(text):
        global continuum_low_reject
        continuum_low_reject = float(text)

    def submit_continuum_high_reject(text):
        global continuum_high_reject
        continuum_high_reject = float(text)

    def normalize(event):
        global spec
        global regions
        xLim = axMain.get_xlim()
        writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])
        fittingFunction = np.polynomial.legendre.legfit
        evalFunction = np.polynomial.legendre.legval
        order = continuum_order
        nIterReject = 2
        nIterFit = 3
        lowReject = continuum_low_reject
        highReject = continuum_high_reject
        useMean = True
#        if len(regions) > 0:
#            wave = []
#            spectrum = []
#            for i in range(len(wLen)):
#                for region in regions:
#                    if (wLen[i] >= xLim[0]) & (wLen[i] <= xLim[1]):
#                        if (wLen[i] >= region[0]) & (wLen[i] <= region[1]):
#                            wave.append(wLen[i])
#                            spectrum.append(spec[i])
#        else:
#            wave = wLen
#            spectrum = spec
        normSpec = continuum(outputSpec,
                            outputSpec,
                            fittingFunction,
                            evalFunction,
                            order,
                            nIterReject,
                            nIterFit,
                            lowReject,
                            highReject,
                            xLim = xLim,
                            regions=regions,
                            #wLen=np.asarray(wave),
                            type='ratio',
                            adjustSigLevels=False,
                            useMean=useMean,
                            display=False)
#        if len(regions) > 0:
#            for i in range(len(wave)):
#                idx = np.where(wLen == wave[i])[0]
#                if idx >= 0:
#                    print('wave[',i,'] = ',wave[i],': idx=',idx,', wLen[',idx,'] = ',wLen[idx],', spec[',idx,'] = ',spec[idx])
#                    spec[idx] = normSpec[i]
#                    print('wave[',i,'] = ',wave[i],': idx=',idx,', wLen[',idx,'] = ',wLen[idx],', spec[',idx,'] set to ',spec[idx])
#        else:
        spec = normSpec
#        spec = getImageData(outputSpec,0)
        axMain.plot(wLen,spec)
        fig.canvas.draw_idle()

    def undo(event):
        global spec_bak
        global wlen_bak
        global spec
        global wLen
        if spec_bak is not None:
            spec = spec_bak.copy()
        if wlen_bak is not None:
            wLen = wlen_bak.copy()
        axMain.plot(wLen,spec)
        fig.canvas.draw_idle()


    fig = plt.figure(figsize=(15,9))
    axMainRect = [0.04,0.26,0.95,0.7]
    axMain = plt.axes(axMainRect)
    axTrimClean = plt.axes([0.01,0.01,0.08,0.1])
    undo_axbox = plt.axes([0.9,0.11,0.1,0.05])
    normalize_axbox = plt.axes([0.01,0.13,0.1,0.05])
    axNorm = plt.axes([0.15,0.11,0.1,0.1])
    continuum_order_axbox = plt.axes([0.44,0.11,0.1,0.05])
    continuum_high_reject_axbox = plt.axes([0.6,0.11,0.1,0.05])
    continuum_low_reject_axbox = plt.axes([0.76,0.11,0.1,0.05])
    min_val=0

    normalize_box = Button(normalize_axbox, 'normalise')
    normalize_box.on_clicked(normalize)

    undo_box = Button(undo_axbox, "undo")
    undo_box.on_clicked(undo)

    continuum_order = 3
    continuum_order_box = TextBox(continuum_order_axbox, 'order', initial=continuum_order)
    continuum_order_box.on_submit(submit_continuum_order)

    continuum_low_reject = 3.
    continuum_low_reject_box = TextBox(continuum_low_reject_axbox, 'low reject', initial=continuum_low_reject)
    continuum_low_reject_box.on_submit(submit_continuum_low_reject)

    continuum_high_reject = 3.
    continuum_high_reject_box = TextBox(continuum_high_reject_axbox, 'high reject', initial=continuum_high_reject)
    continuum_high_reject_box.on_submit(submit_continuum_high_reject)

    spec = getImageData(inputSpec1D,0)
    wLen = getWavelengthArr(inputSpec1D,0)
    spec_bak = np.array(spec)
    wlen_bak = np.array(wLen)
    xRange = []
    regions = []
    cleanType = 'fit'
    normRegion = 'add region'
    newRegion = []

    norm = RadioButtons(axNorm, ('add region','remove region'), active=0)
    def setNormRegion(label):
        global normRegion
        normRegion = label
    norm.on_clicked(setNormRegion)

    radio = RadioButtons(axTrimClean, ('fit', 'fit&remove', 'clean', 'normalise', 'integrate'), active=0)
    def setCleanType(label):
        global cleanType
        cleanType = label

    radio.on_clicked(setCleanType)

    def plotRegions():
        xLim = axMain.get_xlim()
        yLim = axMain.get_ylim()
        axMain.plot(wLen,spec,'b-')
        y = [yLim[0]+(yLim[1]-yLim[0])/30.,yLim[0]+(yLim[1]-yLim[0])/30.]
        axMain.plot(xLim,y,'w-')
        for region in regions:
            x = [region[0],region[1]]
            print('plotting x=',x,', y=',y)
            axMain.plot(x,y,'r-')
        axMain.set_xlim(xLim)
        axMain.set_ylim(yLim)
        #axMain.show()
        fig.canvas.draw_idle()

    def onClick(event):
        global wLen
        global spec
        global xRange
        global newRegion
        global normRegion
        global regions

        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))
        if event.inaxes is axMain:
            print('fig.canvas.toolbar.mode = ',fig.canvas.toolbar.mode)
            if fig.canvas.toolbar.mode == '':
                if cleanType == 'normalise':
                    if normRegion == 'add region':
                        if len(newRegion) == 1:
                            print('old newRegion = ',newRegion)
                            newRegion.append(event.xdata)
                            newRegion = [newRegion[1],newRegion[0]] if (newRegion[1] < newRegion[1]) else [newRegion[0],newRegion[1]]
                            regions.append(newRegion)
                            newRegion = []
                            plotRegions()
                        else:
                            newRegion = [event.xdata]
                        print('newRegion = ',newRegion)
                    else:
                        for iReg in np.arange(len(regions)-1,-1,-1):
                            print('event.xdata = ',event.xdata,', regions[',iReg,'] = ',regions[iReg])
                            if (event.xdata >= regions[iReg][0]) & (event.xdata <= regions[iReg][1]):
                                del regions[iReg]
                        plotRegions()
                    print('regions = ',len(regions),': ',regions)

                elif cleanType == 'trim':
                    if event.xdata < wLen[int(len(wLen)/2.)]:
                        idx = np.where(wLen > event.xdata)
                    else:
                        idx = np.where(wLen < event.xdata)
                    wLen = wLen[idx]
                    spec = spec[idx]
                    axMain.plot(wLen,spec)
                    fig.canvas.draw_idle()
                elif cleanType in ['fit','fit&remove']:
                    if len(xRange) == 1:
                        xRange.append(event.xdata)
                        idx = np.where((wLen >= xRange[0]) & (wLen <= xRange[1]))[0]
                        print('fit: idx = ',idx)
                        p0 = [np.max(spec[idx])-np.min(spec[idx]),
                              wLen[idx[0]]+((wLen[idx[-1]]-wLen[idx[0]])/2.),
                              1.,1.,0.]
                        for p in p0:
                            print('type(',p,') = ',type(p))
                        popt,pcov = curve_fit(gauss_lin,wLen[idx],spec[idx],p0=p0)
                        print('fit: popt=',popt)
                        axMain.plot(wLen[idx],gauss_lin(wLen[idx],*popt))
                        #spec[idx] = spec[idx[0]-1]
                        axMain.plot(wLen,spec)
                        fig.canvas.draw_idle()
                        if cleanType == 'fit&remove':
                            spec[idx] -= gauss_lin(wLen[idx],popt[0],popt[1],popt[2])#,popt[3])
                            spec[idx] += popt[3] * popt[4]*wLen[idx]
                            axMain.plot(wLen,spec)
                            fig.canvas.draw_idle()
                    else:
                        xRange = [event.xdata]
                elif cleanType == 'clean':
                    xLim = axMain.get_xlim()
                    yLim = axMain.get_ylim()
                    print('dir(event) = ',dir(event))
                    print('event = ',event)
#                    if len(xRange) == 1:
#                        xRange.append(event.xdata)
                    idx = np.where(np.abs(wLen - event.xdata) < (wLen[1]-wLen[0])/2.)[0]
                    print('clean: idx = ',idx)
                    spec[idx] = event.ydata
                    axMain.plot(wLen,spec)
                    axMain.set_xlim(xLim)
                    axMain.set_ylim(yLim)
                    fig.canvas.draw_idle()
#                    else:
#                        xRange = [event.xdata]
                elif cleanType == 'integrate':
                    if len(xRange) == 1:
                        xRange.append(event.xdata)
                        idx = np.where((wLen >= xRange[0]) & (wLen <= xRange[1]))[0]
                        ew = 0.
                        for i in range(len(idx)-1):
                            ew += (wLen[idx[i+1]] - wLen[idx[i]]) * (1.-(spec[idx[i+1]]+spec[idx[i]])/2.)
                            coordsPoly = [(wLen[idx[i]],1.),(wLen[idx[i]],spec[idx[i]]),(wLen[idx[i+1]],spec[idx[i+1]]),(wLen[idx[i+1]],1.)]
                            p=Polygon(coordsPoly,facecolor='g')
                            axMain.add_patch(p)
                            fig.canvas.draw_idle()
                        print('equivalent width = ',ew)
                        xRange = []
                    else:
                        xRange = [event.xdata]
                else:
                    print('ERROR: cleanType <'+cleanType+'> not supported')

    cid = fig.canvas.mpl_connect('button_press_event', onClick)

    axMain.plot(wLen,spec)
    axMain.plot(wLen,np.ones(len(wLen)))
    plt.title(inputSpec1D[inputSpec1D.rfind('/')+1:])
    plt.show()
#    writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])


def cleanSpec(inputSpec1D, outputSpec, inputSpec2D=None):
    from matplotlib.widgets import AxesWidget, RadioButtons, Slider, TextBox, Button
    import matplotlib.colors as colors

    global cleanType
    global wLen
    global spec
    global xRange
    global continuum_order
    global continuum_high_reject
    global continuum_low_reject
    global spec_bak
    global wlen_bak
    spec_bak = None
    wlen_bak = None

    def submit_continuum_order(text):
        global continuum_order
        continuum_order = int(text)

    def submit_continuum_low_reject(text):
        global continuum_low_reject
        continuum_low_reject = float(text)

    def submit_continuum_high_reject(text):
        global continuum_high_reject
        continuum_high_reject = float(text)

    def normalize(event):
        global spec
        global wLen
        writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])
        fittingFunction = np.polynomial.legendre.legfit
        evalFunction = np.polynomial.legendre.legval
        order = continuum_order
        nIterReject = 2
        nIterFit = 3
        lowReject = continuum_low_reject
        highReject = continuum_high_reject
        useMean = True
        continuum(outputSpec,
                  outputSpec,
                  fittingFunction,
                  evalFunction,
                  order,
                  nIterReject,
                  nIterFit,
                  lowReject,
                  highReject,
                  type='ratio',
                  adjustSigLevels=False,
                  useMean=useMean,
                  display=False)
        spec = getImageData(outputSpec,0)
        axMain.plot(wLen,spec)
        fig.canvas.draw_idle()

    def remove_continuum(event):
        global spec
        writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])
        fittingFunction = np.polynomial.legendre.legfit
        evalFunction = np.polynomial.legendre.legval
        order = continuum_order
        nIterReject = 2
        nIterFit = 3
        lowReject = continuum_low_reject
        highReject = continuum_high_reject
        useMean = True
        continuum(outputSpec,
                  outputSpec,
                  fittingFunction,
                  evalFunction,
                  order,
                  nIterReject,
                  nIterFit,
                  lowReject,
                  highReject,
                  type='difference',
                  adjustSigLevels=False,
                  useMean=useMean,
                  display=False)
        spec = getImageData(outputSpec,0)
        axMain.plot(wLen,spec)
        fig.canvas.draw_idle()

    def undo(event):
        global spec_bak
        global wlen_bak
        global spec
        global wLen
        if spec_bak is not None:
            spec = spec_bak.copy()
        if wlen_bak is not None:
            wLen = wlen_bak.copy()
        axMain.plot(wLen,spec)
        fig.canvas.draw_idle()


    fig = plt.figure(figsize=(15,9))
    axMainRect = [0.04,0.26,0.95,0.5]
    ax2DRect = plt.axes([0.04,0.77,0.95,0.2])
    axMain = plt.axes(axMainRect)
    axTrimClean = plt.axes([0.01,0.01,0.08,0.1])
    axVMin = plt.axes([0.15,0.06,0.84,0.04])
    axVMax = plt.axes([0.15,0.01,0.84,0.04])
    do_continuum_axbox = plt.axes([0.01,0.11,0.1,0.05])
    undo_axbox = plt.axes([0.9,0.11,0.1,0.05])
    normalize_axbox = plt.axes([0.2,0.11,0.1,0.05])
    continuum_order_axbox = plt.axes([0.44,0.11,0.1,0.05])
    continuum_high_reject_axbox = plt.axes([0.6,0.11,0.1,0.05])
    continuum_low_reject_axbox = plt.axes([0.76,0.11,0.1,0.05])
    max_val=0
    min_val=0

    do_continuum_box = Button(do_continuum_axbox, 'remove continuum')
    do_continuum_box.on_clicked(remove_continuum)

    normalize_box = Button(normalize_axbox, 'normalise')
    normalize_box.on_clicked(normalize)

    undo_box = Button(undo_axbox, "undo")
    undo_box.on_clicked(undo)

    continuum_order = 9
    continuum_order_box = TextBox(continuum_order_axbox, 'order', initial=continuum_order)
    continuum_order_box.on_submit(submit_continuum_order)

    continuum_low_reject = 3.
    continuum_low_reject_box = TextBox(continuum_low_reject_axbox, 'low reject', initial=continuum_low_reject)
    continuum_low_reject_box.on_submit(submit_continuum_low_reject)

    continuum_high_reject = 3.
    continuum_high_reject_box = TextBox(continuum_high_reject_axbox, 'high reject', initial=continuum_high_reject)
    continuum_high_reject_box.on_submit(submit_continuum_high_reject)

    spec = getImageData(inputSpec1D,0)
    wLen = getWavelengthArr(inputSpec1D,0)
    spec_bak = np.array(spec)
    wlen_bak = np.array(wLen)
    if os.path.exists(inputSpec2D):
        image = getImageData(inputSpec2D,0)
        vmax = np.max([2. * np.mean(image), 1.])
    else:
        image = None
        vmax=1.
    vmin = 0.#np.min([1.5 * np.mean(image), 1.])
    if image is not None:
        im = ax2DRect.imshow(image, origin='lower', norm=colors.SymLogNorm(linthresh=0.5, linscale=1,
                             vmin=vmin, vmax=vmax, base=10))#,vmin=vmin, vmax=vmax)#,cmap='gist_rainbow')
    ax2DRect.set_title(inputSpec1D[inputSpec1D.rfind('/')+1:inputSpec1D.rfind('.')])
    xRange = []
    cleanType = 'trim'

    radio = RadioButtons(axTrimClean, ('trim', 'clean'), active=0)
    def setCleanType(label):
        global cleanType
        cleanType = label

    radio.on_clicked(setCleanType)

    def update_max(val):
        max_val=val
    #    axMain.clear()
        if image is not None:
            im.set_clim(vmax=max_val)
            fig.canvas.draw_idle()

    def update_min(val):
        if image is not None:
            min_val=val
            im.set_clim(vmin=min_val)
            fig.canvas.draw_idle()


    d_val=(vmax-vmin)/1000

    #axcolor = 'lightgoldenrodyellow'
    #axfreq = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)
    vMaxS = Slider(axVMax, 'vmax', vmin, vmax,valfmt='% .2f', valinit=0, valstep=d_val)
    vMaxS.on_changed(update_max)
    vMaxS.reset()
    #axfreq1 = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor=axcolor)
    vMinS = Slider(axVMin, 'vmin', vmin, vmax,valfmt='% .2f', valinit=0, valstep=d_val)
    vMinS.on_changed(update_min)
    vMinS.reset()
    max_val=vmax
    min_val=vmin
    vMaxS.set_val(vmax)
    vMinS.set_val(min_val)

    def onClick(event):
        global wLen
        global spec
        global xRange
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))
        if event.inaxes is axMain:
            print('fig.canvas.toolbar.mode = ',fig.canvas.toolbar.mode)
            if fig.canvas.toolbar.mode == '':

                if cleanType == 'trim':
                    if event.xdata < wLen[int(len(wLen)/2.)]:
                        idx = np.where(wLen > event.xdata)
                    else:
                        idx = np.where(wLen < event.xdata)
                    wLen = wLen[idx]
                    spec = spec[idx]
                    axMain.plot(wLen,spec)
                    fig.canvas.draw_idle()
                elif cleanType == 'clean':
                    if len(xRange) == 1:
                        xRange.append(event.xdata)
                        idx = np.where((wLen > xRange[0]) & (wLen < xRange[1]))[0]
                        print('clean: idx = ',idx)
                        spec[idx] = spec[idx[0]-1]
                        axMain.plot(wLen,spec)
                        fig.canvas.draw_idle()
                    else:
                        xRange = [event.xdata]
                else:
                    print('ERROR: cleanType <'+cleanType+'> not supported')

    cid = fig.canvas.mpl_connect('button_press_event', onClick)

    axMain.plot(wLen,spec)
    plt.show()
    if outputSpec is not None:
        writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])


# --- iSpec=[0...nSpec-1]
def getWavelengthMultiSpec(header, iSpec, axis=1):
    dispStr = header['WAT2_001']
    iWat = 2
    while(True):
        try:
            watStr = 'WAT2_%03i' % (iWat)
#            print('getWavelengthMultiSpec: watStr = <'+watStr+'>')
            dispStr += header[watStr]
            iWat += 1
        except:
            break
#    print('getWavelengthMultiSpec: dispStr = <'+dispStr+'>')
#    a='82 1 0 4275.4858499095 2.1134804589752 1333 0. 49.83 56.22'
#    b='2 0 0 4275.4858499095 2.1134804589752 1333 0. 603.94 610.34'
    startStr = 'spec2'
    endStr = '"'
    specDispStr = dispStr[dispStr.find(startStr):]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDispStr = specDispStr[specDispStr.find(endStr)+1:]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDispStr = specDispStr[:specDispStr.find(endStr)]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDisp = specDispStr.split(' ')
#    print('getWavelengthMultiSpec: specDisp = ',specDisp)
    nDigitsBehindDot = []
    for dispVal in specDisp:
#        print(dispVal+'.find(".") = ',dispVal.find('.'))
        if dispVal.find('.') >= 0:
#            print('getWavelengthMultiSpec: dispVal = <'+dispVal+'>: found .')
            nDigitsBehindDot.append(len(dispVal) - dispVal.find('.') - 1)
        else:
#            print('getWavelengthMultiSpec: dispVal = <'+dispVal+'>: did not find .')
            if len(nDigitsBehindDot) == 0:
                nDigitsBehindDot.append(-2)
            else:
                nDigitsBehindDot.append(0-len(dispVal))
#    print('getWavelengthMultiSpec: nDigitsBehindDot = ',nDigitsBehindDot)

    startStr = 'spec'+str(iSpec+1)
    endStr = '"'
    specDispStr = dispStr[dispStr.find(startStr):]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDispStr = specDispStr[specDispStr.find(endStr)+1:]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDispStr = specDispStr[:specDispStr.find(endStr)]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDisp = specDispStr.split(' ')
#    print('getWavelengthMultiSpec: specDisp = ',specDisp)
    if len(specDisp) != len(nDigitsBehindDot):
        specDispNew = []
        iNDig = 0
        for iDispVal in range(len(specDisp)):
            dispVal = specDisp[iDispVal]
#            print('getWavelengthMultiSpec: iNDig = ',iNDig)
#            print('dispVal.find(".") = ',dispVal.find('.'))
#            print('nDigitsBehindDot[',iNDig,'] < 0 = ',nDigitsBehindDot[iNDig] < 0)
            if (dispVal.find('.') < 0) or ((dispVal.find('.') > 0) and (nDigitsBehindDot[iNDig] < 0)):
                if len(dispVal) > (0-nDigitsBehindDot[iNDig]):
                    specDispNew.append(dispVal[0:0-nDigitsBehindDot[iNDig]])
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
                    specDispNew.append(dispVal[0-nDigitsBehindDot[iNDig]:])
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
                    iNDig += 1
#                    print('getWavelengthMultiSpec: iNDig = ',iNDig)
                else:
                    specDispNew.append(dispVal)
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
            else:
#                print('getWavelengthMultiSpec: nDigitsBehindDot[',iNDig,'] = ',nDigitsBehindDot[iNDig])
                if nDigitsBehindDot[iNDig] < (len(dispVal) - dispVal.find('.') - 1):
#                    print('getWavelengthMultiSpec: wrong length: nDigitsBehindDot[iNDig] = ',nDigitsBehindDot[iNDig],', (len(dispVal)(=',len(dispVal),') - dispVal.find(".")(=',dispVal.find('.'),' - 1) = ',len(dispVal) - dispVal.find('.') - 1)
                    specDispNew.append(dispVal[0:dispVal.find('.')+1+nDigitsBehindDot[iNDig]])
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
                    specDispNew.append(dispVal[dispVal.find('.')+1+nDigitsBehindDot[iNDig]:])
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
                    iNDig += 1
#                    print('getWavelengthMultiSpec: iNDig = ',iNDig)
                else:
                    specDispNew.append(dispVal)
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
            iNDig += 1
        specDisp = specDispNew
#    print('getWavelengthMultiSpec: specDisp = ',specDisp)
    nPix = int(header['NAXIS'+str(axis)])
    crVal = float(specDisp[3])
    cDelt = float(specDisp[4])
    lam = np.arange(crVal, crVal + (nPix*cDelt), cDelt, dtype=np.float32)
    if lam.shape[0] > nPix:
        lam = lam[0:nPix]
#    print 'getWavelength: lam = ',len(lam),': ',lam
    return lam


def fitEmissionLines(inputSpec1D,outputSpec):
    from matplotlib.widgets import RadioButtons, TextBox, Button
    from matplotlib.patches import Polygon

    global wLen
    global spec
    global xRange
    global continuum_order
    global continuum_high_reject
    global continuum_low_reject
    global spec_bak
    global wlen_bak
    global cleanType
    global newRegion
    global normRegion
    global regions
    spec_bak = None
    wlen_bak = None

    def submit_continuum_order(text):
        global continuum_order
        continuum_order = int(text)

    def submit_continuum_low_reject(text):
        global continuum_low_reject
        continuum_low_reject = float(text)

    def submit_continuum_high_reject(text):
        global continuum_high_reject
        continuum_high_reject = float(text)

    def normalize(event):
        global spec
        global regions
        xLim = axMain.get_xlim()
        writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])
        fittingFunction = np.polynomial.legendre.legfit
        evalFunction = np.polynomial.legendre.legval
        order = continuum_order
        nIterReject = 2
        nIterFit = 3
        lowReject = continuum_low_reject
        highReject = continuum_high_reject
        useMean = True
#        if len(regions) > 0:
#            wave = []
#            spectrum = []
#            for i in range(len(wLen)):
#                for region in regions:
#                    if (wLen[i] >= xLim[0]) & (wLen[i] <= xLim[1]):
#                        if (wLen[i] >= region[0]) & (wLen[i] <= region[1]):
#                            wave.append(wLen[i])
#                            spectrum.append(spec[i])
#        else:
#            wave = wLen
#            spectrum = spec
        print('regions = ',regions)
        normSpec = continuum(outputSpec,
                            outputSpec,
                            fittingFunction,
                            evalFunction,
                            order,
                            nIterReject,
                            nIterFit,
                            lowReject,
                            highReject,
                            xLim = xLim,
                            regions=regions,
                            #wLen=np.asarray(wave),
                            type='ratio',
                            adjustSigLevels=False,
                            useMean=useMean,
                            display=False)
#        if len(regions) > 0:
#            for i in range(len(wave)):
#                idx = np.where(wLen == wave[i])[0]
#                if idx >= 0:
#                    print('wave[',i,'] = ',wave[i],': idx=',idx,', wLen[',idx,'] = ',wLen[idx],', spec[',idx,'] = ',spec[idx])
#                    spec[idx] = normSpec[i]
#                    print('wave[',i,'] = ',wave[i],': idx=',idx,', wLen[',idx,'] = ',wLen[idx],', spec[',idx,'] set to ',spec[idx])
#        else:
        spec = normSpec
#        spec = getImageData(outputSpec,0)
        axMain.plot(wLen,spec)
        fig.canvas.draw_idle()

    def undo(event):
        global spec_bak
        global wlen_bak
        global spec
        global wLen
        if spec_bak is not None:
            spec = spec_bak.copy()
        if wlen_bak is not None:
            wLen = wlen_bak.copy()
        axMain.plot(wLen,spec)
        fig.canvas.draw_idle()


    fig = plt.figure(figsize=(15,9))
    axMainRect = [0.04,0.26,0.95,0.7]
    axMain = plt.axes(axMainRect)
    axTrimClean = plt.axes([0.01,0.01,0.1,0.1])
    undo_axbox = plt.axes([0.9,0.11,0.1,0.05])
    normalize_axbox = plt.axes([0.01,0.13,0.1,0.05])
    axNorm = plt.axes([0.15,0.11,0.1,0.1])
    continuum_order_axbox = plt.axes([0.44,0.11,0.1,0.05])
    continuum_high_reject_axbox = plt.axes([0.6,0.11,0.1,0.05])
    continuum_low_reject_axbox = plt.axes([0.76,0.11,0.1,0.05])
    min_val=0

    normalize_box = Button(normalize_axbox, 'normalise')
    normalize_box.on_clicked(normalize)

    undo_box = Button(undo_axbox, "undo")
    undo_box.on_clicked(undo)

    continuum_order = 3
    continuum_order_box = TextBox(continuum_order_axbox, 'order', initial=continuum_order)
    continuum_order_box.on_submit(submit_continuum_order)

    continuum_low_reject = 3.
    continuum_low_reject_box = TextBox(continuum_low_reject_axbox, 'low reject', initial=continuum_low_reject)
    continuum_low_reject_box.on_submit(submit_continuum_low_reject)

    continuum_high_reject = 3.
    continuum_high_reject_box = TextBox(continuum_high_reject_axbox, 'high reject', initial=continuum_high_reject)
    continuum_high_reject_box.on_submit(submit_continuum_high_reject)

    spec = getImageData(inputSpec1D,0)
    wLen = getWavelengthArr(inputSpec1D,0)
    spec_bak = np.array(spec)
    wlen_bak = np.array(wLen)
    xRange = []
    regions = []
    cleanType = 'fit'
    normRegion = 'add region'
    newRegion = []

    norm = RadioButtons(axNorm, ('add region','remove region'), active=0)
    def setNormRegion(label):
        global normRegion
        normRegion = label
    norm.on_clicked(setNormRegion)

    radio = RadioButtons(axTrimClean, ('fit', 'fit&remove', 'clean', 'normalise', 'integrate', 'integrate norm'), active=0)
    def setCleanType(label):
        global cleanType
        cleanType = label

    radio.on_clicked(setCleanType)

    def plotRegions():
        xLim = axMain.get_xlim()
        yLim = axMain.get_ylim()
        axMain.plot(wLen,spec,'b-')
        y = [yLim[0]+(yLim[1]-yLim[0])/30.,yLim[0]+(yLim[1]-yLim[0])/30.]
        axMain.plot(xLim,y,'w-')
        for region in regions:
            x = [region[0],region[1]]
            print('plotting x=',x,', y=',y)
            axMain.plot(x,y,'r-')
        axMain.set_xlim(xLim)
        axMain.set_ylim(yLim)
        #axMain.show()
        fig.canvas.draw_idle()

    def onClick(event):
        global wLen
        global spec
        global xRange
        global newRegion
        global normRegion
        global regions

        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))
        if event.inaxes is axMain:
            print('fig.canvas.toolbar.mode = ',fig.canvas.toolbar.mode)
            if fig.canvas.toolbar.mode == '':
                if cleanType == 'normalise':
                    if normRegion == 'add region':
                        if len(newRegion) == 1:
                            print('old newRegion = ',newRegion)
                            newRegion.append(event.xdata)
                            newRegion = [newRegion[1],newRegion[0]] if (newRegion[1] < newRegion[0]) else [newRegion[0],newRegion[1]]
                            regions.append(newRegion)
                            newRegion = []
                            plotRegions()
                        else:
                            newRegion = [event.xdata]
                        print('newRegion = ',newRegion)
                    else:
                        for iReg in np.arange(len(regions)-1,-1,-1):
                            print('event.xdata = ',event.xdata,', regions[',iReg,'] = ',regions[iReg])
                            if (event.xdata >= regions[iReg][0]) & (event.xdata <= regions[iReg][1]):
                                del regions[iReg]
                        plotRegions()
                    print('regions = ',len(regions),': ',regions)

                elif cleanType == 'trim':
                    if event.xdata < wLen[int(len(wLen)/2.)]:
                        idx = np.where(wLen > event.xdata)
                    else:
                        idx = np.where(wLen < event.xdata)
                    wLen = wLen[idx]
                    spec = spec[idx]
                    axMain.plot(wLen,spec)
                    fig.canvas.draw_idle()
                elif cleanType in ['fit','fit&remove']:
                    if len(xRange) == 1:
                        xRange.append(event.xdata)
                        idx = np.where((wLen >= xRange[0]) & (wLen <= xRange[1]))[0]
                        print('fit: idx = ',idx)
                        p0 = [np.max(spec[idx])-np.min(spec[idx]),
                              wLen[idx[0]]+((wLen[idx[-1]]-wLen[idx[0]])/2.),
                              1.,1.,0.]
                        for p in p0:
                            print('type(',p,') = ',type(p))
                        popt,pcov = curve_fit(gauss_lin,wLen[idx],spec[idx],p0=p0)
                        print('fit: popt=',popt)
                        axMain.plot(wLen[idx],gauss_lin(wLen[idx],*popt))
                        #spec[idx] = spec[idx[0]-1]
                        axMain.plot(wLen,spec)
                        fig.canvas.draw_idle()
                        if cleanType == 'fit&remove':
                            spec[idx] -= gauss_lin(wLen[idx],popt[0],popt[1],popt[2])#,popt[3])
                            spec[idx] += popt[3] * popt[4]*wLen[idx]
                            axMain.plot(wLen,spec)
                            fig.canvas.draw_idle()
                    else:
                        xRange = [event.xdata]
                elif cleanType == 'clean':
                    xLim = axMain.get_xlim()
                    yLim = axMain.get_ylim()
                    print('dir(event) = ',dir(event))
                    print('event = ',event)
#                    if len(xRange) == 1:
#                        xRange.append(event.xdata)
                    idx = np.where(np.abs(wLen - event.xdata) < (wLen[1]-wLen[0])/2.)[0]
                    print('clean: idx = ',idx)
                    spec[idx] = event.ydata
                    axMain.plot(wLen,spec)
                    axMain.set_xlim(xLim)
                    axMain.set_ylim(yLim)
                    fig.canvas.draw_idle()
#                    else:
#                        xRange = [event.xdata]
                elif cleanType == 'integrate':
                    if len(xRange) == 1:
                        xRange.append(event.xdata)
                        idx = np.where((wLen >= xRange[0]) & (wLen <= xRange[1]))[0]
                        p0 = [np.max(spec[idx])-np.min(spec[idx]),
                              wLen[idx[0]]+((wLen[idx[-1]]-wLen[idx[0]])/2.),
                              (wLen[idx[-1]]-wLen[idx[0]])/2.]
                        popt,pcov = curve_fit(gauss,wLen[idx],spec[idx],p0=p0)
                        gaussFit = gauss(wLen[idx],popt[0],popt[1],popt[2])
                        ew = 0.
                        for i in range(len(idx)-1):
                            ew += (wLen[idx[i+1]] - wLen[idx[i]]) * (gaussFit[i+1]+gaussFit[i])/2.
                            coordsPoly = [(wLen[idx[i]],0.),(wLen[idx[i]],gaussFit[i]),(wLen[idx[i+1]],gaussFit[i+1]),(wLen[idx[i+1]],0.)]
                            p=Polygon(coordsPoly,facecolor='g')
                            axMain.add_patch(p)
                            fig.canvas.draw_idle()
                        print('equivalent width = ',ew)
                        xRange = []
                    else:
                        xRange = [event.xdata]
                elif cleanType == 'integrate norm':
                    if len(xRange) == 1:
                        xRange.append(event.xdata)
                        idx = np.where((wLen >= xRange[0]) & (wLen <= xRange[1]))[0]
                        p0 = [np.max(spec[idx])-np.min(spec[idx]),
                              wLen[idx[0]]+((wLen[idx[-1]]-wLen[idx[0]])/2.),
                              (wLen[idx[-1]]-wLen[idx[0]])/2.]
                        popt,pcov = curve_fit(gaussNorm,wLen[idx],spec[idx],p0=p0)
                        gaussFit = gauss(wLen[idx],popt[0],popt[1],popt[2]) + 1.
                        ew = 0.
                        for i in range(len(idx)-1):
                            ew += (wLen[idx[i+1]] - wLen[idx[i]]) * (gaussFit[i+1]+gaussFit[i])/2.
                            coordsPoly = [(wLen[idx[i]],1.),(wLen[idx[i]],gaussFit[i]+1),(wLen[idx[i+1]],gaussFit[i+1]+1),(wLen[idx[i+1]],1.)]
                            p=Polygon(coordsPoly,facecolor='g')
                            axMain.add_patch(p)
                            fig.canvas.draw_idle()
                        print('equivalent width = ',ew)
                        xRange = []
                    else:
                        xRange = [event.xdata]
                else:
                    print('ERROR: cleanType <'+cleanType+'> not supported')

    cid = fig.canvas.mpl_connect('button_press_event', onClick)

    axMain.plot(wLen,spec)
#    axMain.plot(wLen,np.ones(len(wLen)))
    plt.title(inputSpec1D[inputSpec1D.rfind('/')+1:])
    plt.show()
#    writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])
