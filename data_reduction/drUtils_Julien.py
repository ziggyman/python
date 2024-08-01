import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import exp,ndimage
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, Angle, ICRS, LSR
from astropy.time import Time
from PyAstronomy import pyasl
from specutils import Spectrum1D
from specutils.manipulation.resample import FluxConservingResampler,LinearInterpolatedResampler


c0 = 299792.458 # km/s

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


def plotSpec(fitsFileName):
    wLen = getWavelengthArr(fitsFileName)
    spec = getImageData(fitsFileName,0)
    plt.plot(wLen,spec)
    plt.show()


def gauss(x,a,x0,sigma,yBackground=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground


def getClosestInTime(objectTime, arcTimes):
    if len(arcTimes) < 1:
        return None
    print('getClosestInTime: objectTime = ',type(objectTime),': ',objectTime)
    print('getClosestInTime: arcTimec = ',type(arcTimes[0]),': ',arcTimes)
    timeDiffs = np.absolute(np.array(arcTimes) - objectTime)
    print('getClosestInTime: timeDiffs = ',timeDiffs)
    minTimeDiff = np.min(timeDiffs)
    return [np.where(timeDiffs == minTimeDiff)[0][0],minTimeDiff]

def getClosestArcs(fitsFileName, fitsList):
    arcTimes = []
    for iArc in range(len(fitsList)):
        hdulist = pyfits.open(fitsList[iArc])
        head = hdulist[0].header
        keyWord = 'HJD-OBS'
        try:
            arcTimes.append(float(head[keyWord]))
        except:
            try:
                keyWord = 'MJD-OBS'
                arcTimes.append(float(head[keyWord]))
            except:
                keyWord = 'DATE-OBS'
                arcTimeTemp = Time(head[keyWord], format='isot', scale='utc')
                arcTimes.append(arcTimeTemp.mjd)
    arcTimes = np.array(arcTimes)

    print('arcTimes = ',arcTimes)
    hdulist = pyfits.open(fitsFileName)
    headerSc = hdulist[0].header
    print('dir(headerSc) = ',dir(headerSc))
#    for key in headerSc.keys():
#        print('getClosestArc: headerSc[',key,'] = ',headerSc[key])
    try:
        if keyWord == 'DATE-OBS':
            specTimeTemp = Time(headerSc[keyWord], format='isot', scale='utc')
            specTime = specTimeTemp.mjd
        else:
            specTime = float(headerSc[keyWord])
    except:
        print('getClosestArcs: ERROR: keyWord <'+keyWord+'> not found in '+fitsFileName)
        if keyWord == 'MJD-OBS':
            try:
                specTimeTemp = Time(headerSc['DATE-OBS'], format='isot', scale='utc')
                specTime = specTimeTemp.mjd
            except:
                print('getClosestArcs: ERROR: keyWord <DATE-OBS> not found in '+fitsFileName)
                STOP
        if keyWord == 'HJD-OBS':
            try:
                specTimeTemp = Time(headerSc['DATE-OBS'], format='isot', scale='utc')
                specTime = specTimeTemp.hjd
            except:
                print('getClosestArcs: ERROR: keyWord <DATE-OBS> not found in '+fitsFileName)
                STOP
    hdulist.close()
    whereLT = np.where(arcTimes <= specTime)[0]
    whereGT = np.where(arcTimes >= specTime)[0]
    print('getClosestArcs: whereLT = ',whereLT)
    print('getClosestArcs: whereGT = ',whereGT)
    if len(whereLT) == 0:
        closestBefore = None
    else:
        closestBefore = getClosestInTime(specTime, arcTimes[whereLT])
    closestTemp = getClosestInTime(specTime, arcTimes[whereGT])
    if len(whereGT) == 0:
        closestAfter = None
    else:
        closestAfter = [whereGT[closestTemp[0]],closestTemp[1]]
    print('getClosestArcs: closestBefore = ',closestBefore)
    print('getClosestArcs: closestAfter = ',closestAfter)
    return [closestBefore,closestAfter]

#@brief: apply wavelength to extracted science spectra and resample them to linear dispersion
def dispCor(scienceListIn,
            arcListIn,
            wavelengthsOrigIn,
            scienceListOut,
            observatoryLocation,
            keywordRA,
            keywordDEC,
            keywordObsTime,
            doHelioCor=True):
    print('dispCor: scienceListIn = ',len(scienceListIn),': ',scienceListIn)
    print('dispCor: arcListIn = ',len(arcListIn),': ',arcListIn)
    print('dispCor: wavelengthsOrigIn = ',len(wavelengthsOrigIn))
    for iSpec in range(len(scienceListIn)):
        print('dispCor: running dispCor on '+scienceListIn[iSpec])
        hdulist = pyfits.open(scienceListIn[iSpec])
        headerSc = hdulist[0].header
        print('dispCor: headerSc = ',headerSc)
        closestArcs = getClosestArcs(scienceListIn[iSpec],arcListIn)
        print('dispCor: closestArcs = ',closestArcs)
        if closestArcs[0] is None:
            closestArcs[0] = closestArcs[1]
        if closestArcs[1] is None:
            closestArcs[1] = closestArcs[0]

        if (closestArcs[0][1] + closestArcs[1][1]) == 0:
            closestArcs[0][1] = 0.5
            closestArcs[1][1] = 0.5

        spec = getImageData(scienceListIn[iSpec],0)
        print('dispcor: scienceListIn[',iSpec,'] = ',scienceListIn[iSpec])
        print('dispCor: spec = ',spec.shape,': ',spec)
        arcBefore = arcListIn[closestArcs[0][0]]
        print('dispCor: name of closest Arc before = ',arcBefore)
        wLenSpecBefore = wavelengthsOrigIn[closestArcs[0][0]]
        print('dispcor: wLenSpecBefore = ',wLenSpecBefore.shape,': ',wLenSpecBefore)
        factorBefore = closestArcs[0][1] / (closestArcs[0][1]+closestArcs[1][1])
        print('dispcor: factorBefore = ',factorBefore)

        arcAfter = arcListIn[closestArcs[1][0]]
        print('dispCor: name of closest Arc after = ',arcAfter)
        print('len(wavelengthsOrigIn) = ',len(wavelengthsOrigIn))
        wLenSpecAfter = wavelengthsOrigIn[closestArcs[1][0]]
        print('dispcor: wLenSpecAfter = ',wLenSpecAfter)
        factorAfter = closestArcs[1][1] / (closestArcs[0][1]+closestArcs[1][1])
        print('dispcor: factorAfter = ',factorAfter)

        wLenSpec = (wLenSpecBefore * factorBefore) + (wLenSpecAfter * factorAfter)
        print('dispCor before heliocor: wLenSpec = ',wLenSpec)

        #read science header and append keywords
        headerSc['REFSPEC1'] = arcBefore[arcBefore.rfind('/')+1:]+' %.5f' % (factorBefore)
        headerSc['REFSPEC2'] = arcAfter[arcAfter.rfind('/')+1:]+' %.5f' % (factorAfter)

        #apply heliocentric radial velocity correction
        if doHelioCor:
            print('running heliocor for ',scienceListIn[iSpec])
            vrad = heliocor(observatoryLocation, headerSc, keywordRA, keywordDEC, keywordObsTime)
            headerSc['VHELIO'] = vrad
            wLenSpec = applyVRadCorrection(wLenSpec, vrad)
            print('dispCor: after heliocentric correction for vrad = ',vrad,': wLenSpec = ',wLenSpec.shape,': wLenSpec')
        hdulist.close()
        print('after heliocor: wLenSpec = ',wLenSpec)

        # read wavelength information from reference ARCs
        hdulist = pyfits.open(arcBefore)
        headerArcBefore = hdulist[0].header
        hdulist.close()
        hdulist = pyfits.open(arcAfter)
        headerArcAfter = hdulist[0].header
        hdulist.close()

#        resampledBefore = ((np.arange(headerArcBefore['NAXIS1']) + 1.0) - headerArcBefore['CRPIX1']) * headerArcBefore['CDELT1'] + headerArcBefore['CRVAL1']
#        print('resampledBefore = ',resampledBefore)
#        resampledAfter = ((np.arange(headerArcAfter['NAXIS1']) + 1.0) - headerArcAfter['CRPIX1']) * headerArcAfter['CDELT1'] + headerArcAfter['CRVAL1']
#        print('resampledAfter = ',resampledAfter)
#        resampled = (resampledBefore * factorBefore) + (resampledAfter * factorAfter)
#        print('resampled = ',resampled)
#        resampledSpec = rebin(wLenSpec, spec, resampled, preserveFlux = True)
#        print("headerArcBefore['CRVAL1'] = ",headerArcBefore['CRVAL1'])
#        print("headerArcBefore['CRVAL1'] * factorBefore = ",headerArcBefore['CRVAL1'] * factorBefore)
#        print("headerArcAfter['CRVAL1'] = ",headerArcAfter['CRVAL1'])
#        print("headerArcAfter['CRVAL1'] * factorAfter = ",headerArcAfter['CRVAL1'] * factorAfter)
#        print("(headerArcBefore['CRVAL1'] * factorBefore) + (headerArcAfter['CRVAL1'] * factorAfter) = ",(headerArcBefore['CRVAL1'] * factorBefore) + (headerArcAfter['CRVAL1'] * factorAfter))
#
#        print("headerArcBefore['CRPIX1'] = ",headerArcBefore['CRPIX1'])
#        print("headerArcBefore['CRPIX1'] * factorBefore = ",headerArcBefore['CRPIX1'] * factorBefore)
#        print("headerArcAfter['CRPIX1'] = ",headerArcAfter['CRPIX1'])
#        print("headerArcAfter['CRPIX1'] * factorAfter = ",headerArcAfter['CRPIX1'] * factorAfter)
#        print("(headerArcBefore['CRPIX1'] * factorBefore) + (headerArcAfter['CRPIX1'] * factorAfter) = ",(headerArcBefore['CRPIX1'] * factorBefore) + (headerArcAfter['CRPIX1'] * factorAfter))
#
#        print("headerArcBefore['CDELT1'] = ",headerArcBefore['CDELT1'])
#        print("headerArcBefore['CDELT1'] * factorBefore = ",headerArcBefore['CDELT1'] * factorBefore)
#        print("headerArcAfter['CDELT1'] = ",headerArcAfter['CDELT1'])
#        print("headerArcAfter['CDELT1'] * factorAfter = ",headerArcAfter['CDELT1'] * factorAfter)
#        print("(headerArcBefore['CDELT1'] * factorBefore) + (headerArcAfter['CDELT1'] * factorAfter) = ",(headerArcBefore['CDELT1'] * factorBefore) + (headerArcAfter['CDELT1'] * factorAfter))
        #STOP
        resampled,resampledSpec = resampleSpec(wLenSpec,spec)
        print('scienceListOut = ',len(scienceListOut),scienceListOut)
        print('iSpec = ',iSpec)
        writeFits1D(resampledSpec,
                    scienceListOut[iSpec],
                    wavelength=None,
                    header=headerSc,
                    CRVAL1=resampled[0],#(headerArcBefore['CRVAL1'] * factorBefore) + (headerArcAfter['CRVAL1'] * factorAfter),
                    CRPIX1=1,#(headerArcBefore['CRPIX1'] * factorBefore) + (headerArcAfter['CRPIX1'] * factorAfter),
                    CDELT1=resampled[1]-resampled[0],#(headerArcBefore['CDELT1'] * factorBefore) + (headerArcAfter['CDELT1'] * factorAfter),
                   )
        wLenSpecTest = getWavelengthArr(scienceListOut[iSpec])
        print('dispCor: wLenSpecTest = ',wLenSpecTest)
        #if 'dbs00343' in scienceListIn[iSpec]:
        #    STOP

def resampleSpec(wLen, spec):
    dLam = np.min([np.absolute(wLen[1]-wLen[0]),np.absolute(wLen[wLen.shape[0]-1]-wLen[wLen.shape[0]-2])])
    resampled = np.arange(np.min([wLen[0], wLen[wLen.shape[0]-1]]), np.max([wLen[0], wLen[wLen.shape[0]-1]]), dLam)

    print('resampleSpec: wLen = ',wLen.shape,': ',wLen)
    print('resampleSpec: spec = ',spec.shape,': ',spec)
    print('resampleSpec: resampled = ',resampled.shape,': ',resampled)
    resampledSpec = rebin(wLen, spec, resampled, preserveFlux = False)
    return [resampled,resampledSpec]

def rebin(wavelength, spectrum, newWavelength, preserveFlux = True):#, outFileName = None, header = None):
    if wavelength[1] < wavelength[0]:
        wLen = np.fliplr([wavelength])[0]
        spec = np.fliplr([spectrum])[0]
    else:
        wLen = wavelength
        spec = spectrum
    if newWavelength[1] < newWavelength[0]:
        wLenNew = np.fliplr([newWavelength])[0]
    else:
        wLenNew = newWavelength
    print('spec.shape = ',spec.shape,', spec.shape = ',spec.shape)
    input_spectrum = Spectrum1D( flux=np.array(spec) * u.erg / (u.cm * u.cm) / u.s / u.AA,
                                spectral_axis = np.array(wLen) * u.AA)
    print('rebin: spec = ',spec.shape,': ',spec)
    print('rebin: input_spectrum = ',input_spectrum)
    resample_grid = np.array(wLenNew) *u.AA
#        print('rebin: resample_grid = ',resample_grid)
    if preserveFlux:
        fluxc_resample = FluxConservingResampler(extrapolation_treatment='zero_fill')
    else:
        fluxc_resample = LinearInterpolatedResampler(extrapolation_treatment='zero_fill')
    specInterp = fluxc_resample(input_spectrum, resample_grid) # doctest: +IGNORE_OUTPUT
    print('rebin: specInterp = ',specInterp)
    print('rebin: specInterp.data = ',specInterp.data)
    naNPos = np.argwhere(np.isnan(specInterp.data))
    print('rebin: naNPos = ',naNPos)
    if naNPos.shape[0] > 0:
        for pos in naNPos:
            specInterp.data[pos[0]] = 0.
    naNPos = np.argwhere(np.isnan(specInterp.data))
    print('rebin: naNPos = ',naNPos)
    if naNPos.shape[0] > 0:
        STOP

    return specInterp.data

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

def heliocor(observatoryLocation, header, keywordRA, keywordDEC, keywordObsTime):
    #print('heliocor: EarthLocation.get_site_names() = ',EarthLocation.get_site_names())
    #GTC = EarthLocation.of_site('Roque de los Muchachos')
    #print('heliocor: GTC = ',GTC)
    #print('heliocor: GTC.lon = ',GTC.lon)
    #print('heliocor: GTC.lat = ',GTC.lat)
    #print('heliocor: dir(GTC) = ',dir(GTC))

    obsRA = hmsToDeg(header[keywordRA])
    obsDEC = dmsToDeg(header[keywordDEC])
    sc = SkyCoord(ra=obsRA*u.deg, dec=obsDEC*u.deg)
    obsTime = Time(header[keywordObsTime])
    barycorr = sc.radial_velocity_correction(obstime=obsTime, location=observatoryLocation)
    print('heliocor: barycorr = ',barycorr.to(u.km/u.s))
    heliocorr = sc.radial_velocity_correction('heliocentric', obstime=obsTime, location=observatoryLocation)
    print('heliocor: heliocorr.value = ',heliocorr.value)
    print('heliocor: heliocorr.to_value() = ',heliocorr.to_value())
    print('heliocor: heliocorr.to(u.km/u.s) = ',heliocorr.to(u.km/u.s))
    print('heliocor: heliocorr.to(u.km/u.s).value = ',heliocorr.to(u.km/u.s).value)

    return heliocorr.to(u.km/u.s).value


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


def applyVRadCorrection(wavelength, vRad):
    return np.array([wLen - (vRad * wLen / c0) for wLen in wavelength])


def readFileToArr(fname):
    with open(fname, "r") as f:
        lines = f.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut
