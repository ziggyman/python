import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
from astropy.coordinates import SkyCoord

def getImageData(fname,
                 hduNum=0):
    hdulist = pyfits.open(fname)
    scidata = hdulist[hduNum].data
    hdulist.close()
    return scidata

def getHeader(fName,
              hduNum=0):
    hdulist = pyfits.open(fName)
    header = hdulist[hduNum].header
    hdulist.close()
    return header

def getHeaderValue(fname,
                   keyword,
                   hduNum=0):
    hdulist = pyfits.open(fname)
    header = hdulist[hduNum].header
    hdulist.close()
    try:
        return header[keyword]
    except:
        return None

def setHeaderValue(fitsFileName,
                   keyword,
                   value,
                   hduNum=0):
    hdulist = pyfits.open(fitsFileName)
    header = hdulist[hduNum].header
    header[keyword] = value
    hdulist.writeto(fitsFileName,overwrite=True)
    hdulist.close()

def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

def getWavelengthArr(fname,
                     hduNum=0):
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

# NOTE: currently it is required for the x values to be accending
def normalizeX(x):
    xZero = x - x[0]
    return 2.0 * xZero / xZero[-1] - 1.0

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
    print('sigmaReject: iter = 0: indices = ',indices.shape,': ',indices)
    for iter in np.arange(1,nIter,1):
        yGood, goodIndices = sigmaRej(yGood, lowReject, highReject, replace=replace, adjustSigLevels=adjustSigLevels, useMean=useMean, keepFirstAndLastX=keepFirstAndLastX)
        indices = indices[goodIndices]
        print('sigmaReject: iter = ',iter,': yGood = ',yGood.shape,': ',yGood)
        print('sigmaReject: iter = ',iter,': indices = ',indices.shape,': ',indices)
    return [yGood,indices]

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
    fittedValues = solveFunction(x,coeffs)
    fittedIndices = np.arange(0,fittedValues.shape[0],1)
    print('sfit: iter = -1: x.shape = ',x.shape,', y.shape = ',y.shape,', fittedValues.shape = ',fittedValues.shape,', fittedIndices = ',fittedIndices.shape,': ',fittedIndices)
    fittedValuesTemp = fittedValues
    for i in range(nIterFit):
        print('sfit: iter = ',i,': x.shape = ',x.shape,', y.shape = ',y.shape,', fittedValues.shape = ',fittedValues.shape,', fittedIndices = ',fittedIndices.shape,': ',fittedIndices)
        print('sfit: y[fittedIndices] = ',y[fittedIndices])
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
        print('sfit: y[fittedIndices] = ',y[fittedIndices].shape,': ',y[fittedIndices])
        print('sfit: fittedValuesNotRejected = ',fittedValuesNotRejected.shape,': ',fittedValuesNotRejected)
        """x is already required to be normalized at function call"""
        coeffs = fittingFunction(x[fittedIndices],y[fittedIndices],order)
        """keep first and last x values to get the normalization right!!!"""
        fittedValues = solveFunction(x,coeffs)
        print('sfit: sfit: iter = ',i,': x.shape = ',x.shape,', y.shape = ',y.shape,', fittedValues.shape = ',fittedValues.shape,', fittedIndices = ',fittedIndices.shape,': ',fittedIndices)
        fittedValuesTemp = fittedValues[fittedIndices]
        if display:
            print('sfit: sfit: iter = ',i,': y[fittedIndices] = ',y[fittedIndices],', fittedValuesNotRejected = ',fittedValuesNotRejected)
            plt.plot(x[fittedIndices],y[fittedIndices] - fittedValuesTemp,'b*', label='fitted points')
            plt.plot(x[fittedIndices],y[fittedIndices],'r+', label='fitted points')
            plt.plot(x,fittedValues,label='new fit')
            plt.legend()
            plt.show()
    return [coeffs, fittedValues]

#            fittingFunction = np.polynomial.legendre.legfit
#            evalFunction = np.polynomial.legendre.legval
#            order = 5
#            nIterReject = 2
#            nIterFit = 3
#            lowReject = 2.
#            highReject = 1.8
#            useMean = True

def continuum(spectrumFileNameIn,
              spectrumFileNameOut,
              fittingFunction = np.polynomial.legendre.legfit,
              evalFunction = np.polynomial.legendre.legval,
              order = 5,
              nIterReject = 2,
              nIterFit = 3,
              lowReject = 2.,
              highReject = 1.8,
              type = 'difference',
              adjustSigLevels = False,
              useMean = True,
              display = True):
    if isinstance(spectrumFileNameIn,str):
        print('continuum: reading file '+spectrumFileNameIn)
        specOrig = getImageData(spectrumFileNameIn,0)
    else:
        specOrig = spectrumFileNameIn

    if display:
        wLen = getWavelengthArr(spectrumFileNameIn,0)
        plt.plot(wLen, specOrig,label='original spectrum')
        plt.title(spectrumFileNameIn[spectrumFileNameIn.rfind('/')+1:spectrumFileNameIn.rfind('.')])
        plt.legend()
        plt.show()

    xNorm = normalizeX(np.arange(0,specOrig.shape[0],1.))
    #sfit(x, y, fittingFunction, solveFunction, order, nIterReject, nIterFit, lowReject, highReject, adjustSigLevels=False, useMean=False, display=False)
    coeffs, resultFit = sfit(xNorm,
                             specOrig,
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
    if type == 'difference':
        spec = specOrig - resultFit
    elif type == 'ratio':
        spec = specOrig / resultFit
    elif type == 'fit':
        spec = resultFit
    else:
        print('continuum: ERROR: type <'+type+'> not recognised')
        STOP

    if display:
        plt.plot(wLen, specOrig, label='original')
        plt.plot(wLen, spec, label = 'continuum corrected')
        plt.legend()
        xRange = [wLen[0],wLen[len(wLen)-1]]
        yRange = [np.min([np.min(specOrig), np.min(spec)]),np.max([np.max(specOrig), np.max(spec)])]
        plt.text(xRange[1],yRange[0],spectrumFileNameOut[spectrumFileNameOut.rfind('/')+1:spectrumFileNameOut.rfind('.')],rotation='vertical')
#        markEmissionLines(xRange, yRange)
        plt.show()

    writeFits1D(spec,
                spectrumFileNameOut,
                wavelength=None,
                header=spectrumFileNameIn,
                CRVAL1=getHeaderValue(spectrumFileNameIn,'CRVAL1'),
                CRPIX1=getHeaderValue(spectrumFileNameIn,'CRPIX1'),
                CDELT1=getHeaderValue(spectrumFileNameIn,'CDELT1'),
               )
    print(spectrumFileNameOut,' written')


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

def plotSpec(fitsFileName):
    wLen = getWavelengthArr(fitsFileName)
    spec = getImageData(fitsFileName,0)
    plt.plot(wLen,spec)
    plt.show()

# read header from inputFileName, add metaKeys and metaData to that header,
# adjust the size according to ccdData, and write ccdData and header to outputFileName
def writeFits(spectrum,
              inputFileName,
              outputFileName,
              metaKeys=None,
              metaData=None,
              overwrite=True):
    hdulist = pyfits.open(inputFileName)
    print('writeFits: len(hdulist) = ',len(hdulist),', inputFileName = ',inputFileName,', outputFileName = ',outputFileName,', old mean = ',hdulist[len(hdulist)-1].data,', mean(ccdData) = ',np.mean(spectrum))
    print('len(hdulist) = ',len(hdulist))
    hdulist[len(hdulist)-1].data = spectrum
    hdulist[len(hdulist)-1].header['NAXIS1'] = np.asarray(spectrum).shape[0]
    if metaKeys is not None:
        for iKey in np.arange(0,len(metaKeys),1):
            hdulist[len(hdulist)-1].header[metaKeys[iKey]] = metaData[iKey]
#    print('hdulist[len(hdulist)-1].header.keys() = ',hdulist[len(hdulist)-1].header.keys())
    if 'OBSERVER' in hdulist[len(hdulist)-1].header.keys():
        hdulist[len(hdulist)-1].header['OBSERVER'] = 'Quentin Parker + Travis Stenborg'
        print("hdulist[",len(hdulist)-1,"].header['OBSERVER'] = ",hdulist[len(hdulist)-1].header['OBSERVER'])
    print('writeFits: mean(hdulist[',len(hdulist)-1,'].data) = ',np.mean(hdulist[len(hdulist)-1].data))
    hdulist.writeto(outputFileName, overwrite=overwrite)
    print('writeFits: new mean after writing image = ',np.mean(getImageData(outputFileName,len(hdulist)-1)))

def trimSpec(inputFileName,
             wLenStart,
             wLenEnd,
             outputFileName):
    spectrum = getImageData(inputFileName)
    wLen = getWavelengthArr(inputFileName)

    goodPoints = np.where((wLen >= wLenStart) & (wLen <= wLenEnd))
    print('goodPoints = ',len(goodPoints),': ',goodPoints)
    wLenNew = wLen[goodPoints]
    print('wLen = ',len(wLen),': ',wLen)
    print('wLenNew = ',len(wLenNew),': ',wLenNew)
    spectrumNew = spectrum[goodPoints]
    writeFits(spectrumNew,
              inputFileName=inputFileName,
              outputFileName=outputFileName,
              metaKeys=['CRVAL1'],
              metaData=[wLenNew[0]],
              overwrite=True)

def removeInterval(inputFileName,
                   wLenStart,
                   wLenEnd,
                   outputFileName):
    spectrum = getImageData(inputFileName)
    wLen = getWavelengthArr(inputFileName)

    badPoints = np.where((wLen >= wLenStart) & (wLen <= wLenEnd))
    print('badPoints = ',len(badPoints),': ',badPoints)
    spectrum[badPoints] = spectrum[badPoints[0][0]-1]
    writeFits(spectrum,
              inputFileName=inputFileName,
              outputFileName=outputFileName,
              overwrite=True)

# all angles must be in degrees
#@return: angular distance in degrees
def angularDistancePyAsl(ra1, dec1, ra2, dec2):
    from PyAstronomy import pyasl
#    print('ra1 = ',ra1,', dec1 = ',dec1,', ra2 = ',ra2,', dec2 = ',dec2)
    return pyasl.getAngDist(ra1, dec1, ra2, dec2)


#@param ra1: RA in degrees
#@param dec1: DEC in degrees
#@param ra2: RA in degrees
#@param dec2: DEC in degrees
#@return: angular distance in arc seconds
def angularDistance(ra1, dec1, ra2, dec2):
    if ':' in str(ra1):
        print('values must be in degrees!')
        STOP
#    print('ra1 = ',ra1,', dec1 = ',dec1,', ra2 = ',ra2,', dec2 = ',dec2)
    mm1 = SkyCoord(ra=ra1, dec=dec1, frame='icrs', unit="deg")
    mm2 = SkyCoord(ra=ra2, dec=dec2, frame='icrs', unit="deg")
    return mm1.separation(mm2).arcsecond


# RA string = xx:yy:zz.zzz
def hmsToDeg(string):
    h, m, s = [float(i) for i in string.split(':')]
    return (15. * s / 3600.) + (15. * m / 60.) + (h * 15.)

def degToHMS(degrees):
    h = int(degrees / 15.)
    m = int((degrees - (h * 15.)) * 4.)
    s = (degrees - (m/4.) - (h*15.)) * 240.
    sStr = '%.3f' % (s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (h,m,sStr)

# DEC string = xx:yy:zz.zzz
def dmsToDeg(string):
    d, m, s = [float(i) for i in string.split(':')]
#    print('dmsToDeg: string = <'+string+'>: d = ',d,', m = ',m,', s = ',s)
    if string[0] == '-':
        d = 0. - d
        return 0. - (s / 3600. + m / 60. + d)
    return s / 3600. + m / 60. + d

def degToDMS(degrees):
    d = int(degrees)
    m = int((degrees - d) * 60)
    s = (degrees - d - (m/60.)) * 3600.
    sStr = '%.3f' % abs(s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (d,abs(m),sStr)
