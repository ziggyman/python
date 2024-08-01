from astropy import units as u
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.legendre import legval
import os
from astropy.coordinates import SkyCoord
import re
import subprocess
from PyAstronomy import pyasl
from scipy.optimize import curve_fit
from specutils import Spectrum1D
from specutils.manipulation.resample import FluxConservingResampler,LinearInterpolatedResampler
from scipy import exp,ndimage

c0 = 299792.458 # km/s

def getRaDecXY(string):
    strs = re.sub( '\s+', ' ', string ).strip()
#    print('getRaDecXY: strs = ',strs)
    strs = strs.rstrip().split(' ')
#    print('getRaDecXY: strs = ',strs)
    return [strs[0], strs[1], float(strs[len(strs)-2]), float(strs[len(strs)-1])]

#@brief convert RA and DEC to Galactic Longitude and Latitude
#@param ra: RA in degrees
#@param dec: DEC in degrees
#@return [lon, lat]
def raDecToLonLat(ra, dec):
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
#    print('c.galactic = ',type(c.galactic),': ',c.galactic)
#    print('dir(c) = ',type(c),': ',dir(c))
#    print('dir(c.galactic) = ',type(c.galactic),': ',dir(c.galactic))
#    print('c.galactic.l = ',type(c.galactic.l),': ',dir(c.galactic.l))
#    print('c.galactic.l.value = ',type(c.galactic.l.value),': ',dir(c.galactic.l.value))
    return [float(c.galactic.l.value),float(c.galactic.b.value)]

def lonLatToRaDec(l, b):
    c = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
    return [float(c.icrs.ra.value),float(c.icrs.dec.value)]

#def lonLatToRaDec(lon, lat):
#    c = SkyCoord(frame="galactic", l="1h12m43.2s", b="+1d12m43s")
#    return [float(c.ra.value),float(c.dec.value)]

#def getPixel(hammerX, hammerY):
#    ham = hammer.Hammer()
#    for pix in ham.getPixels():
#        if ham.isInside(pix, hammerX, hammerY):
#            return pix
#    return None

# string = xx:yy:zz.zzz
def hmsToDeg(string):
    try:
        h, m, s = [float(i) for i in string.split(':')]
        return (15. * s / 3600.) + (15. * m / 60.) + (h * 15.)
    except:
        print('hmsToDeg: string = <'+string+'>')
        STOP

def degToHMS(degrees):
    h = int(degrees / 15.)
    m = int((degrees - (h * 15.)) * 4.)
    s = (degrees - (m/4.) - (h*15.)) * 240.
    sStr = '%.3f' % (s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (h,m,sStr)

# string = xx:yy:zz.zzz
def dmsToDeg(string):
    try:
        d, m, s = [float(i) for i in string.split(':')]
    #    print('dmsToDeg: string = <'+string+'>: d = ',d,', m = ',m,', s = ',s)
        if string[0] == '-':
            d = 0. - d
            return 0. - (s / 3600. + m / 60. + d)
        return s / 3600. + m / 60. + d
    except:
        return float(string)

def degToDMS(degrees):
    d = int(degrees)
    m = int((degrees - d) * 60)
    s = (degrees - d - (m/60.)) * 3600.
    sStr = '%.3f' % abs(s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (d,abs(m),sStr)

def hmsToCoordStr(string):
    h, m, s = [i for i in string.split(':')]
    return h+'h'+m+'m'+s+'s'

def dmsToCoordStr(string):
    d, m, s = [i for i in string.split(':')]
    return d+'d'+m+'m'+s+'s'

def getPNG(lon,lat):
    return str(lon)[:str(lon).find('.')+2].zfill(5)+('+'if lat >= 0. else '-')+str(abs(lat))[:str(abs(lat)).find('.')+2].zfill(4)

# all angles must be in degrees
#@return: angular distance in degrees
def angularDistancePyAsl(ra1, dec1, ra2, dec2):
    from PyAstronomy import pyasl
#    print('ra1 = ',ra1,', dec1 = ',dec1,', ra2 = ',ra2,', dec2 = ',dec2)
    return pyasl.getAngDist(ra1, dec1, ra2, dec2)

def angularDistanceFromXYPyAsl(fitsName, x1, y1, x2, y2):
    result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x1), str(y1)])
    result2 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x2), str(y2)])
#    print('xy2sky(',str(x1),', ',str(y1),') = ',result1)
#    print('xy2sky(',str(x2),', ',str(y2),') = ',result2)
    raHMS1, decHMS1, x1, y1 = getRaDecXY(result1.decode('utf-8'))
    raHMS2, decHMS2, x2, y2 = getRaDecXY(result2.decode('utf-8'))
#    print('raHMS1 = ',raHMS1,', decHMS1 = ',decHMS1)
#    print('raHMS2 = ',raHMS2,', decHMS2 = ',decHMS2)
    ra1 = hmsToDeg(raHMS1)
    dec1 = dmsToDeg(decHMS1)
    ra2 = hmsToDeg(raHMS2)
    dec2 = dmsToDeg(decHMS2)
#    print('ra1 = ',ra1,', dec1 = ',dec1)
#    print('ra2 = ',ra2,', dec2 = ',dec2)
    return angularDistance(ra1, dec1, ra2, dec2)

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

def angularDistanceFromXY(fitsName, x1, y1, x2, y2):
    result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x1), str(y1)])
    result2 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x2), str(y2)])
#    print('xy2sky(',str(x1),', ',str(y1),') = ',result1)
#    print('xy2sky(',str(x2),', ',str(y2),') = ',result2)
    raHMS1, decHMS1, x1, y1 = getRaDecXY(result1.decode('utf-8'))
    raHMS2, decHMS2, x2, y2 = getRaDecXY(result2.decode('utf-8'))
#    print('raHMS1 = ',raHMS1,', decHMS1 = ',decHMS1)
#    print('raHMS2 = ',raHMS2,', decHMS2 = ',decHMS2)
    ra1 = hmsToCoordStr(raHMS1)
    dec1 = dmsToCoordStr(decHMS1)
    ra2 = hmsToCoordStr(raHMS2)
    dec2 = dmsToCoordStr(decHMS2)
#    print('ra1 = ',ra1,', dec1 = ',dec1)
#    print('ra2 = ',ra2,', dec2 = ',dec2)
    return angularDistance(ra1, dec1, ra2, dec2)

def getXYFromRaDec(fitsName, raHMS, decDMS):
    print('getXYFromRaDec: fitsName = ',fitsName,', raHMS = ',raHMS,', decDMS = ',decDMS)
    result = subprocess.check_output(['sky2xy', '-j', fitsName, raHMS, decDMS])
    print('getXYFromRaDec: result = <',result,'>')
    strs = re.sub( '\s+', ' ', result.decode('utf-8') ).strip()
    print('getXYFromRaDec: strs = ',strs)
    strs = strs.rstrip().split(' ')
    print('getXYFromRaDec: strs = ',strs)
    return [float(strs[4]), float(strs[5])]

def getRaDecFromXY(fitsName, x, y):
    result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x), str(y)])
    raHMS, decHMS, x1, y1 = getRaDecXY(result1)
    return [raHMS, decHMS]

def getArcsecDistance(fitsName, x1, y1, x2, y2):
    result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x1), str(y1)])
    result2 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x2), str(y2)])
    print('xy2sky(',str(x1),', ',str(y1),') = ',result1)
    print('xy2sky(',str(x2),', ',str(y2),') = ',result2)
    raHMS1, decHMS1, x1, y1 = getRaDecXY(result1)
    raHMS2, decHMS2, x2, y2 = getRaDecXY(result2)
    #print('raHMS1 = ',raHMS1,', decHMS1 = ',decHMS1)
    #print('raHMS2 = ',raHMS2,', decHMS2 = ',decHMS2)
    mm1 = SkyCoord(ra=raHMS1, dec=decHMS1, unit=(u.hourangle, u.deg))
    mm2 = SkyCoord(ra=raHMS2, dec=decHMS2, unit=(u.hourangle, u.deg))
    return mm1.separation(mm2).arcsecond

def degToArcsec(deg):
    return deg * 3600.


def readFileToArr(fname):
    with open(fname, "r") as f:
        lines = f.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut


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

def setHeaderValue(fitsFileName,keyword,value,hduNum=0):
    hdulist = pyfits.open(fitsFileName)
    header = hdulist[hduNum].header
    header[keyword] = value
    hdulist.writeto(fitsFileName,overwrite=True)
    hdulist.close()

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


# read header from inputFileName, add metaKeys and metaData to that header,
# adjust the size according to ccdData, and write ccdData and header to outputFileName
def writeFits(ccdData, inputFileName, outputFileName, metaKeys=None, metaData=None, overwrite=True):
    hdulist = pyfits.open(inputFileName)
    print('writeFits: len(hdulist) = ',len(hdulist),', inputFileName = ',inputFileName,', outputFileName = ',outputFileName,', old mean = ',hdulist[len(hdulist)-1].data,', mean(ccdData) = ',np.mean(ccdData.data))
    print('len(hdulist) = ',len(hdulist))
    hdulist[len(hdulist)-1].data = ccdData.data
    hdulist[len(hdulist)-1].header['NAXIS1'] = np.asarray(ccdData.data).shape[1]
    hdulist[len(hdulist)-1].header['NAXIS2'] = np.asarray(ccdData.data).shape[0]
    if metaKeys is not None:
        for iKey in np.arange(0,len(metaKeys),1):
            hdulist[len(hdulist)-1].header[metaKeys[iKey]] = metaData[iKey]
#    print('hdulist[len(hdulist)-1].header.keys() = ',hdulist[len(hdulist)-1].header.keys())
    if 'OBSERVER' in hdulist[len(hdulist)-1].header.keys():
        hdulist[len(hdulist)-1].header['OBSERVER'] = 'Quentin Parker + Travis Stenborg'
        print("hdulist[",len(hdulist)-1,"].header['OBSERVER'] = ",hdulist[len(hdulist)-1].header['OBSERVER'])
    print('writeFits: mean(hdulist[',len(hdulist)-1,'].data) = ',np.mean(hdulist[len(hdulist)-1].data))
#    print('hdulist[',len(hdulist)-1,'].header = ',hdulist[len(hdulist)-1].header)
#    print('hdulist[len(hdulist)-1].header[TRIMSEC] = ',hdulist[len(hdulist)-1].header['TRIMSEC'])
    hdulist.writeto(outputFileName, overwrite=overwrite)
    print('writeFits: new mean after writing image = ',np.mean(getImageData(outputFileName,len(hdulist)-1)))

def gauss(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))


def gauss_const(x,a,x0,sigma,yBackground=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground


def gauss_lin(x,a,x0,sigma,yBackground=0.,linear=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground+(linear*x)


# method = ['preserveFlux','linear','spline']
def rebin_Spectrum1D(wavelength, spectrum, newWavelength, method='preserveFlux'):#, outFileName = None, header = None):
    from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler, SplineInterpolatedResampler
    wLenOld = wavelength * u.AA
    flux = spectrum * u.Unit('erg cm-2 s-1 AA-1')
    input_spec = Spectrum1D(spectral_axis = wLenOld, flux=flux)
    wLenNew = newWavelength * u.AA
    if method == 'preserveFlux':
        fluxcon = FluxConservingResampler()
        new_spec = fluxcon(input_spec, wLenNew)
    elif method == 'linear':
        linear = LinearInterpolatedResampler()
        new_spec = linear(input_spec, wLenNew)
    elif method == 'spline':
        spline = SplineInterpolatedResampler()
        new_spec = spline(input_spec, wLenNew)
#    print('dir(new_spec) = ',dir(new_spec))
#    print('dir(new_spec.flux) = ',dir(new_spec.flux))
    ret = new_spec.flux.to_value()
#    print('ret = ',ret)
#    STOP
    return ret


def rebin_spec(wave, specin, wavnew):
    from pysynphot import observation
    from pysynphot import spectrum
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')

    return obs.binflux


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


def applyVRadCorrection(wavelength, vRad):
    return np.array([wLen - (vRad * wLen / c0) for wLen in wavelength])

def chisqfunc(fac, object, sky, sigma):
    model = fac * sky
    chisq = np.sum(((object - model) / sigma)**2)
    return chisq

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


# NOTE: currently it is required for the x values to be accending
def normalizeX(x):
    xZero = x - x[0]
    return 2.0 * xZero / xZero[-1] - 1.0

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
        #print('cotinuum: wLenFit = ',wLenFit)
        #print('cotinuum: wLenFit.shape = ',wLenFit.shape)
        for i in range(wLenFit.shape[0]):
            for region in regions:
                if (wLenFit[i] >= xLim[0]) & (wLenFit[i] <= xLim[1]):
                    if (wLenFit[i] >= region[0]) & (wLenFit[i] <= region[1]):
                        wave.append(wLenFit[i])
                        spectrum.append(specFit[i])
    else:
        wave = wLen
        spectrum = specOrig

#    if wLen is not None:
#        range = wLen
    #print('continuum: wave = ',len(wave),': ',wave)
    #print('continuum: spectrum = ',len(spectrum),': ',spectrum)

    xToNorm = [xLim[0]]
    for wav in wave:
        xToNorm.append(wav)
    xToNorm.append(xLim[1])
    #print('continuum: xToNorm = ',len(xToNorm),': ',xToNorm)
    xNorm = normalizeX(np.array(xToNorm))
    xNorm = xNorm[1:]
    xNorm = xNorm[:len(xNorm)-1]
    #print('continuum: xNorm = ',len(xNorm),': ',xNorm)
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


def scombine(fileListName,
            spectrumFileNameOut,
            normalise = False,
            method='median', # mean, median
            lowReject=None,
            highReject=None,
            adjustSigLevels=False,
            useMean=False,
            display=False):

    fileNames = readFileToArr(fileListName)
    fileNames = [fileName if '/' in fileName else os.path.join(fileListName[:fileListName.rfind('/')],fileName) for fileName in fileNames]

    if normalise:
        snrA = measureSNR(fileNames[0])
        snrB = measureSNR(fileNames[1])

        fittingFunction = np.polynomial.legendre.legfit
        evalFunction = np.polynomial.legendre.legval
        order = 19
        nIterReject = 2
        nIterFit = 3
        lowReject = 3
        highReject = 3
        useMean = True

        continuum(fileNames[0],
                  fileNames[0][:fileNames[0].rfind('.')]+'_norm.fits',
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
        fileNames[0] = fileNames[0][:fileNames[0].rfind('.')]+'_norm.fits'


        continuum(fileNames[1],
                  fileNames[1][:fileNames[1].rfind('.')]+'_norm.fits',
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
        fileNames[1] = fileNames[1][:fileNames[1].rfind('.')]+'_norm.fits'

    specA = getImageData(fileNames[0],0)
    wLenA = getWavelengthArr(fileNames[0],0)
    spectrumA = Spectrum1D( flux=np.array(specA) * u.erg / (u.cm * u.cm) / u.s / u.AA,
                                spectral_axis = np.array(wLenA) * u.AA)

    specB = getImageData(fileNames[1],0)
    wLenB = getWavelengthArr(fileNames[1],0)
    spectrumB = Spectrum1D( flux=np.array(specB) * u.erg / (u.cm * u.cm) / u.s / u.AA,
                                spectral_axis = np.array(wLenB) * u.AA)

    new_spectral_axis = np.arange(np.min([wLenA[0],wLenB[0]]),np.max([wLenA[len(wLenA)-1],wLenB[len(wLenB)-1]]),np.min([wLenA[1]-wLenA[0],wLenB[1]-wLenB[0]]))*spectrumA.spectral_axis.unit#np.concatenate([spectrumA.spectral_axis.value, spectrumB.spectral_axis.to_value(spectrumA.spectral_axis.unit)]) * spectrumA.spectral_axis.unit
    print('new_spectral_axis = ',new_spectral_axis)
    resampler = LinearInterpolatedResampler(extrapolation_treatment='zero_fill')

    new_spec1 = resampler(spectrumA, new_spectral_axis)
    new_spec2 = resampler(spectrumB, new_spectral_axis)

    final_spec = new_spec1 + new_spec2
    if display:
        plt.plot(final_spec.spectral_axis, final_spec.flux)
        plt.title(fileNames[0])
        plt.show()

    print('dir(final_spec.spectral_axis) = ',dir(final_spec.spectral_axis))
    print('dir(final_spec.spectral_axis.data) = ',dir(final_spec.spectral_axis.data))

    print('type(final_spec.spectral_axis.data[0]) = ',type(final_spec.spectral_axis.data[0]))

    print('type(final_spec) = ',type(final_spec))
    print('dir(final_spec) = ',dir(final_spec))

    #STOP
#    if os.path.exists(spectrumFileNameOut):
#        os.remove(spectrumFileNameOut)
#    final_spec.write(spectrumFileNameOut,format='tabular-fits')


    #dispersion = np.array(final_spec.spectral_axis.data)
    #plt.plot(np.array(final_spec.spectral_axis.data)[1:]-np.array(final_spec.spectral_axis.data)[0:-1],'b+')
    #plt.show()

    if True:
        if normalise:
            spectrumFileNameOut = spectrumFileNameOut[:spectrumFileNameOut.rfind('.')]+'_SNR%d+%d.fits' % (int(snrA),int(snrB))
        spectrumFileNameOutTemp = spectrumFileNameOut
        i = 'a'
        while os.path.exists(spectrumFileNameOutTemp):
            spectrumFileNameOutTemp = spectrumFileNameOut[:spectrumFileNameOut.rfind('.')]+i+'.fits'
            if i == 'a':
                i='b'
            elif i == 'b':
                i='c'
            elif i == 'c':
                i='d'
            elif i == 'd':
                i='e'
            elif i == 'e':
                i='f'
            else:
                print('scombine: i=f already exists')
        writeFits1D(final_spec.flux.data,
                    spectrumFileNameOutTemp,
                    wavelength=None,
                    header=getHeader(fileNames[0],0),
                    CRVAL1=final_spec.spectral_axis.data[0],
                    CRPIX1=1,
                    CDELT1=final_spec.spectral_axis.data[1]-final_spec.spectral_axis.data[0],
                )
        #STOP
    if False:
        wLenAll = getWavelengthArr(fileNames[0],0)
        print('scombine: wLenAll = ',wLenAll)
        spectra = []
        exptimes = []
        for fileName in fileNames:
            exptimes.append(float(getHeaderValue(fileName,'EXPTIME')))
            if fileName == fileNames[0]:
                spectra.append(getImageData(fileName,0))
            else:
                spectra.append(rebin(getWavelengthArr(fileName,0), getImageData(fileName,0), wLenAll, preserveFlux = True))
            if display:
                plt.plot(wLenAll, spectra[len(spectra)-1], label=fileName[fileName.rfind('/')+1:fileName.rfind('.')])
        exptimes = np.array(exptimes)
        print('scombine: exptimes = ',exptimes)
        goodPix = []
        combinedSpectrum = np.zeros(wLenAll.shape[0])
        exptimesPix = []
        if lowReject is not None:
            for pix in np.arange(0,wLenAll.shape[0],1):
                vPix, indices = sigmaReject([spectrum[pix] for spectrum in spectra],
                                            nIter=1,
                                            lowReject=lowReject,
                                            highReject=highReject,
                                            replace=False,
                                            adjustSigLevels=adjustSigLevels,
                                            useMean=useMean,
                                            keepFirstAndLastX=False)

                print('scombine: vPix = ',vPix)
                goodPix.append(vPix)
                print('scombine: exptimes[indices = ',indices,'] = ',exptimes[indices])
                exptimesPix.append(np.mean(exptimes[indices]))
        else:
            for pix in np.arange(0,wLenAll.shape[0],1):
                goodPix.append([spectrum[pix] for spectrum in spectra])
                exptimesPix.append(np.mean(exptimes))
        for pix in np.arange(0,wLenAll.shape[0],1):
            notNaN = np.where([not np.isnan(a) for a in goodPix[pix]])
            if (len(goodPix[pix]) % 2 == 0) or (method == 'mean'):
                print('scombine: pix = ',pix,': notNaN = ',notNaN)
                combinedSpectrum[pix] = np.mean(np.array(goodPix[pix]))#[notNaN][0]))
            else:
                combinedSpectrum[pix] = np.median(np.array(goodPix[pix]))#[notNaN][0]))
            #print('scombine: pix = ',pix,': wLenAll[',pix,'] = ',wLenAll[pix],': goodPix = ',goodPix[pix],', np.mean(np.array(goodPix[pix])) = ',np.mean(np.array(goodPix[pix])),', np.median(np.array(goodPix[pix])) = ',np.median(np.array(goodPix[pix])),', combinedSpectrum[',pix,'] = ',combinedSpectrum[pix])

        print('scombine: exptimes = ',exptimes)
        header = getHeader(fileNames[0])
        print('scombine: exptimesPix = ',exptimesPix)
        header['EXPTIME'] = np.mean(exptimesPix)
        print('scombine: set header[EXPTIME] to ',header['EXPTIME'])

        if display:
            plt.plot(wLenAll, combinedSpectrum, 'g-', label = 'combined')
            plt.legend()
            plt.show()

        writeFits1D(combinedSpectrum,
                    spectrumFileNameOut,
                    wavelength=None,
                    header=header,
                    CRVAL1=getHeaderValue(fileNames[0],'CRVAL1'),
                    CRPIX1=getHeaderValue(fileNames[0],'CRPIX1'),
                    CDELT1=getHeaderValue(fileNames[0],'CDELT1'),
                )


def plotSpec(fitsFileName):
    if fitsFileName[0] == '@':
        fitsFileNames = readFileToArr(fitsFileName[1:])
        for fName in fitsFileNames:
            if '/' not in fName:
                if '/' in fName:
                    fName = os.path.join(fitsFileName[1:fitsFileName.rfind('/')],fName)
            print('fName = ',fName)
            wLen = getWavelengthArr(fName,0)
            spec = getImageData(fName,0)
            plt.plot(wLen,spec)
            plt.title(fName[fName.rfind('/')+1:])
            plt.show()
    else:
        wLen = getWavelengthArr(fitsFileName)
        spec = getImageData(fitsFileName,0)
        plt.plot(wLen,spec)
        plt.show()

def cleanSpec(inputSpec1D, inputSpec2D, outputSpec):
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
    writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])


def measureSNR(spectrumFileNameIn):
    spec = getImageData(spectrumFileNameIn,0)
    snr = []
    for i in np.arange(10,len(spec)-10,1):
        snr.append(np.mean(spec[i-10:i+10] / np.std(spec[i-10:i+10])))
    plt.plot(snr)
    medianSNR = np.median(snr)
    print(spectrumFileNameIn,': median SNR = ',medianSNR)
    plt.show()
    return medianSNR

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
                else:
                    print('ERROR: cleanType <'+cleanType+'> not supported')

    cid = fig.canvas.mpl_connect('button_press_event', onClick)

    axMain.plot(wLen,spec)
#    axMain.plot(wLen,np.ones(len(wLen)))
    plt.title(inputSpec1D[inputSpec1D.rfind('/')+1:])
    plt.show()
#    writeFits1D(spec,outputSpec,wavelength=None,header=getHeader(inputSpec1D,0), CRVAL1=wLen[0], CRPIX1=1, CDELT1=wLen[1]-wLen[0])

def fitEmissionLine(wLen, spec, xRange, Display=False):
    if Display:
        from matplotlib.patches import Polygon
    idx = np.where((wLen >= xRange[0]) & (wLen <= xRange[1]))[0]
    if Display:
        fig,ax = plt.subplots()
        plt.plot(wLen[idx],spec[idx])
    p0 = [np.max(spec[idx])-np.min(spec[idx]),
            wLen[idx[0]]+((wLen[idx[-1]]-wLen[idx[0]])/2.),
            (wLen[idx[-1]]-wLen[idx[0]])/2.]
    popt,pcov = curve_fit(gauss,wLen[idx],spec[idx],p0=p0)
    gaussFit = gauss(wLen[idx],popt[0],popt[1],popt[2])
    areaUnderGauss = 0.
    for i in range(len(idx)-1):
        areaUnderGauss += (wLen[idx[i+1]] - wLen[idx[i]]) * (gaussFit[i+1]+gaussFit[i])/2.
        if Display:
            coordsPoly = [(wLen[idx[i]],0.),(wLen[idx[i]],gaussFit[i]),(wLen[idx[i+1]],gaussFit[i+1]),(wLen[idx[i+1]],0.)]
            p=Polygon(coordsPoly,facecolor='g')
            ax.add_patch(p)
    print('equivalent width = ',areaUnderGauss)
    if Display:
        plt.show()
    return areaUnderGauss
