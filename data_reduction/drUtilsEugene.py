import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from PyAstronomy import pyasl
from scipy.optimize import curve_fit
from scipy import exp





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

def continuum(spectrumFileNameIn,
            spectrumFileNameOut,
            fittingFunction,
            evalFunction,
            order,
            nIterReject,
            nIterFit,
            lowReject,
            highReject,
            wLen = None,
            xLim=None,
            regions=None,
            type='difference',
            adjustSigLevels=False,
            useMean=False,
            display = False):
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
        xLim = [wLen[0],wLen(len(wLen)-1)]

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
        print('cotinuum: wLenFit = ',wLenFit)
        print('cotinuum: wLenFit.shape = ',wLenFit.shape)
        for i in range(wLenFit.shape[0]):
            for region in regions:
                if (wLenFit[i] >= xLim[0]) & (wLenFit[i] <= xLim[1]):
                    if (wLenFit[i] >= region[0]) & (wLenFit[i] <= region[1]):
                        wave.append(wLenFit[i])
                        spectrum.append(specFit[i])
    else:
        wave = wLen
        spectrum = spec

#    if wLen is not None:
#        range = wLen

    xNorm = normalizeX(wave)
    print('continuum: xNorm = ',xNorm)
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
    fittedSpectrum = evalFunction(normalizeX(xrange),coeffs)
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
        plt.plot(wLen, specDone, label = 'continuum corrected')
        plt.legend()
        xRange = [wLen[0],wLen[len(wLen)-1]]
        yRange = [np.min([np.min(specOrig), np.min(specDone)]),np.max([np.max(specOrig), np.max(specDone)])]
        if spectrumFileNameOut is not None:
            plt.text(xRange[1],yRange[0],spectrumFileNameOut[spectrumFileNameOut.rfind('/')+1:spectrumFileNameOut.rfind('.')],rotation='vertical')
        markEmissionLines(xRange, yRange)
        plt.show()

    if spectrumFileNameOut is not None:
        writeFits1D(specDone,
                    spectrumFileNameOut,
                    wavelength=None,
                    header=spectrumFileNameIn,
                    CRVAL1=getHeaderValue(spectrumFileNameIn,'CRVAL1'),
                    CRPIX1=getHeaderValue(spectrumFileNameIn,'CRPIX1'),
                    CDELT1=getHeaderValue(spectrumFileNameIn,'CDELT1'),
                )
    for i in range(len(specOrig)-1):
        print('specOrig[',i,'] = ',specOrig[i],', specDone[',i,'] = ',specDone[i])
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
