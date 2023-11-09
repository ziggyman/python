import numpy as np
from astropy import units as u
from specutils import Spectrum1D
#from specutils.manipulation.resample import FluxConservingResampler,LinearInterpolatedResampler
from matplotlib import pyplot as plt
from scipy import exp
from scipy.optimize import curve_fit
import astropy.io.fits as pyfits
from PyAstronomy import pyasl
import csv

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

def gauss(x,a,x0,sigma,yBackground=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground

def gauss_lin(x,a,x0,sigma,yBackground=0.,linear=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground+(linear*x)

def xCorFindMinimum(xCorX, xCorY, display = False):
    print('xCorFindMinimum: xCorX = ',xCorX)
    print('xCorFindMinimum: xCorY = ',xCorY)
    y = xCorY - np.amax(xCorY)
    a = np.amin(y)
    print('xCorFindMinimum: a = ',a)
    print('xCorFindMinimum: y = ',y,', np.amin(y) = ',np.amin(y))
    whereMin = np.where(y == a)[0]
    print('xCorFindMinimum: whereMin = ',whereMin)
    x0 = xCorX[whereMin][0]
    print('xCorFindMinimum: x0 = ',type(x0),': ',x0)
    print('xCorX = ',type(xCorX),': ',xCorX)
    print('y = ',type(y),': ',y)
    try:
        p0 = [a,x0,1.,0.,0.]
        for p in p0:
            print('type(',p,') = ',type(p))
        popt,pcov = curve_fit(gauss_lin,xCorX,y,p0=p0)
        where = np.where(y == np.min(y))[0]
        print('popt[1](=',popt[1],') - xCorX[',where[0],'](=',xCorX[where[0]],') = ',popt[1] - xCorX[where[0]])
        if np.abs(popt[1] - xCorX[where[0]]) > (1.5*np.abs(xCorX[1]-xCorX[0])):
            print('difference between bestVRad = ',popt[1],' and vRadRange[',where[0],']=',xCorX[where[0]],' gt dVRad')
            STOP
        if display:
            plt.plot(xCorX,y,label='y')
            plt.plot(xCorX,gauss_lin(xCorX,popt[0],popt[1],popt[2],popt[3],popt[4]),label='gauss')
            plt.legend()
            plt.show()
    except:
        p0 = [a,x0,1.,0.]
        for p in p0:
            print('type(',p,') = ',type(p))
        popt,pcov = curve_fit(gauss,xCorX,y,p0=p0)
        where = np.where(y == np.min(y))[0]
        if np.abs(popt[1] - xCorX[where[0]]) > (1.5*np.abs(xCorX[1]-xCorX[0])):
            print('difference between bestVRad = ',popt[1],' and vRadRange[',where[0],']=',xCorX[where[0]],' gt dVRad')
            STOP
        if display:
            plt.plot(xCorX,y,label='y')
            plt.plot(xCorX,gauss(xCorX,popt[0],popt[1],popt[2],popt[3]),label='gauss')
            plt.legend()
            plt.show()
    return popt

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

def applyVRadCorrection(wavelength, vRad):
    return np.array([wLen - (vRad * wLen / c0) for wLen in wavelength])

# trim spectra to same wavelength range and rebin spectrum with higher dLambda to smaller dLambda
def rebinAndTrimToSameWavelengthRangeAndDispersion(pnSpecWLen, pnSpecData, pnSpecCompWLen, pnSpecCompData, display=False):
    """check wavelength ranges, trim spectra, rebin pnSpec"""
    pnSpecWLenStart = pnSpecWLen[0]
    pnSpecWLenEnd = pnSpecWLen[len(pnSpecWLen)-1]

    pnSpecCompWLenStart = pnSpecCompWLen[0]
    pnSpecCompWLenEnd = pnSpecCompWLen[len(pnSpecCompWLen)-1]

    pnSpecWhere = np.where((pnSpecWLen > np.max([pnSpecWLenStart,pnSpecCompWLenStart])) & (pnSpecWLen < np.min([pnSpecWLenEnd,pnSpecCompWLenEnd])))[0]
    print('rebinAndTrimToSameWavelengthRangeAndDispersion: pnSpecWhere = ',pnSpecWhere)
    pnSpecCompWhere = np.where((pnSpecCompWLen > np.max([pnSpecWLenStart,pnSpecCompWLenStart])) & (pnSpecCompWLen < np.min([pnSpecWLenEnd,pnSpecCompWLenEnd])))[0]
    print('rebinAndTrimToSameWavelengthRangeAndDispersion: pnSpecCompWhere = ',pnSpecCompWhere)

    pnSpecWLen = pnSpecWLen[pnSpecWhere]
    pnSpecData = pnSpecData[pnSpecWhere]
    print('rebinAndTrimToSameWavelengthRangeAndDispersion: pnSpecWLen = ',pnSpecWLen.shape,': [',pnSpecWLen[0],',...,',pnSpecWLen[pnSpecWLen.shape[0]-1],']')

    pnSpecCompWLen = pnSpecCompWLen[pnSpecCompWhere]
    pnSpecCompData = pnSpecCompData[pnSpecCompWhere]
    print('rebinAndTrimToSameWavelengthRangeAndDispersion: pnSpecCompWLen = ',pnSpecCompWLen.shape,': [',pnSpecCompWLen[0],',...,',pnSpecCompWLen[pnSpecCompWLen.shape[0]-1],']')

    pnSpecDeltaLam = np.min(np.array([pnSpecWLen[1] - pnSpecWLen[0],pnSpecWLen[pnSpecWLen.shape[0]-1] - pnSpecWLen[pnSpecWLen.shape[0]-2]]))
    pnSpecCompDeltaLam = np.min(np.array([pnSpecCompWLen[1] - pnSpecCompWLen[0],pnSpecCompWLen[pnSpecCompWLen.shape[0]-1] - pnSpecCompWLen[pnSpecCompWLen.shape[0]-2]]))
    print('rebinAndTrimToSameWavelengthRangeAndDispersion: pnSpecDeltaLam = ',pnSpecDeltaLam)
    print('rebinAndTrimToSameWavelengthRangeAndDispersion: pnSpecCompDeltaLam = ',pnSpecCompDeltaLam)

    wLenSame = np.arange(np.min(pnSpecWLen),np.max(pnSpecWLen),np.min(np.array([pnSpecDeltaLam,pnSpecCompDeltaLam])))
#    pnSpecCompData_rF = rebin(pnSpecCompWLen,pnSpecCompData,wLenSame,preserveFlux=True)
#    pnSpecData_rF = rebin(pnSpecWLen,pnSpecData,wLenSame,preserveFlux=True)

#    pnSpecCompData_linear = rebin(pnSpecCompWLen,pnSpecCompData,wLenSame,preserveFlux=False)
#    pnSpecData_linear = rebin(pnSpecWLen,pnSpecData,wLenSame,preserveFlux=False)

#    pnSpecCompData_rebin1D_preserveFlux = rebin_Spectrum1D(pnSpecCompWLen,pnSpecCompData,wLenSame, method='preserveFlux')
#    pnSpecData_rebin1D_preserveFlux = rebin_Spectrum1D(pnSpecWLen,pnSpecData,wLenSame, method='preserveFlux')

#    pnSpecCompData_rebin1D_linear = rebin_Spectrum1D(pnSpecCompWLen,pnSpecCompData,wLenSame, method='linear')
#    pnSpecData_rebin1D_linear = rebin_Spectrum1D(pnSpecWLen,pnSpecData,wLenSame, method='linear')

    pnSpecCompData_rebin1D_spline = rebin_Spectrum1D(pnSpecCompWLen,pnSpecCompData,wLenSame, method='spline')
    pnSpecData_rebin1D_spline = rebin_Spectrum1D(pnSpecWLen,pnSpecData,wLenSame, method='spline')

#    pnSpecCompData_rebin1D_spec = rebin_spec(pnSpecCompWLen,pnSpecCompData,wLenSame)
#    pnSpecData_rebin1D_spec = rebin_spec(pnSpecWLen,pnSpecData,wLenSame)

    #pnSpecWLen = pnSpecCompWLen
    if display:
#        plt.plot(pnSpecWLen,pnSpecData,label='pnSpec original')
#        plt.plot(pnSpecCompWLen,pnSpecCompData,label='pnSpecComp original')
#        plt.plot(wLenSame,pnSpecData_rF,label='pnSpec rebinned F')
#        plt.plot(wLenSame,pnSpecCompData_rF,label='pnSpecComp rebinned F')
#        plt.legend()
#        plt.show()

#        plt.plot(pnSpecWLen,pnSpecData,label='pnSpec original')
#        plt.plot(pnSpecCompWLen,pnSpecCompData,label='pnSpecComp original')
#        plt.plot(wLenSame,pnSpecData_linear,label='pnSpec rebinned linear')
#        plt.plot(wLenSame,pnSpecCompData_linear,label='pnSpecComp rebinned linear')
#        plt.legend()
#        plt.show()

#        plt.plot(pnSpecWLen,pnSpecData,label='pnSpec original')
#        plt.plot(pnSpecCompWLen,pnSpecCompData,label='pnSpecComp original')
#        plt.plot(wLenSame,pnSpecData_rebin1D_preserveFlux,label='pnSpec rebinned1dF')
#        plt.plot(wLenSame,pnSpecCompData_rebin1D_preserveFlux,label='pnSpecComp rebinned1dF')
#        plt.legend()
#        plt.show()

#        plt.plot(pnSpecWLen,pnSpecData,label='pnSpec original')
#        plt.plot(pnSpecCompWLen,pnSpecCompData,label='pnSpecComp original')
#        plt.plot(wLenSame,pnSpecData_rebin1D_linear,label='pnSpec rebinned1dL')
#        plt.plot(wLenSame,pnSpecCompData_rebin1D_linear,label='pnSpecComp rebinned1dL')
#        plt.legend()
#        plt.show()

        plt.plot(pnSpecWLen,pnSpecData,label='pnSpec original')
        plt.plot(pnSpecCompWLen,pnSpecCompData,label='pnSpecComp original')
        plt.plot(wLenSame,pnSpecData_rebin1D_spline,label='pnSpec rebinned1dS')
        plt.plot(wLenSame,pnSpecCompData_rebin1D_spline,label='pnSpecComp rebinned1dS')
        plt.legend()
        plt.show()

#        plt.plot(pnSpecWLen,pnSpecData,label='pnSpec original')
#        plt.plot(pnSpecCompWLen,pnSpecCompData,label='pnSpecComp original')
#        plt.plot(wLenSame,pnSpecData_rebin1D_spec,label='pnSpec rebinned1d spec')
#        plt.plot(wLenSame,pnSpecCompData_rebin1D_spec,label='pnSpecComp rebinned1d spec')
#        plt.legend()
#        plt.show()
    return [wLenSame,pnSpecData_rebin1D_spline,pnSpecCompData_rebin1D_spline]

def getRadialVelocityFromXCor(pnSpec, pnSpecComp, vRadSpecComp = 0.):
    """read and normalize spectra to maximum of 1.0"""
    pnSpecData = getImageData(pnSpec,0)
    pnSpecData = pnSpecData / np.max(pnSpecData)
    pnSpecCompData = getImageData(pnSpecComp,0)
    pnSpecCompData = pnSpecCompData / np.max(pnSpecCompData)

    pnSpecWLen = getWavelengthArr(pnSpec,0)
    print('getRadialVelocityFromXCor: pnSpecWLen = ',pnSpecWLen.shape,': [',pnSpecWLen[0],',...,',pnSpecWLen[pnSpecWLen.shape[0]-1],']')
    pnSpecCompWLen = getWavelengthArr(pnSpecComp,0)
    print('getRadialVelocityFromXCor: pnSpecCompWLen = ',pnSpecCompWLen.shape,': [',pnSpecCompWLen[0],',...,',pnSpecCompWLen[pnSpecCompWLen.shape[0]-1],']')
    if vRadSpecComp != 0.:
        plt.plot(pnSpecCompWLen,pnSpecCompData,label='original')
        pnSpecCompWLen = applyVRadCorrection(pnSpecCompWLen,0.-vRadSpecComp)
        plt.plot(pnSpecCompWLen,pnSpecCompData,label='vrad=0')
        plt.plot([6562.81,6562.81],[0.,1.])
        plt.legend()
        plt.show()
    plt.plot(pnSpecWLen,pnSpecData,label='original')
    plt.plot([6562.81,6562.81],[0.,1.])
    plt.legend()
    plt.show()

    vRadRange = np.arange(-500.,500.,5.)
    chiSquares = []
    for vRad in vRadRange:
        wLenVRad  = applyVRadCorrection(pnSpecCompWLen, 0.-vRad)
        wLen,pnSpecDataRebinned,pnSpecCompDataRebinned = rebinAndTrimToSameWavelengthRangeAndDispersion(pnSpecWLen, pnSpecData, wLenVRad, pnSpecCompData, True if vRad == vRadRange[int(len(vRadRange)/2.)] else False)

        wLenNaNPos = np.argwhere(np.isnan(wLen))
        print('getRadialVelocityFromXCor: vRad = ',vRad,': wLenNaNPos = ',wLenNaNPos)
        if wLenNaNPos.shape[0] > 0:
            STOP

        pnSpecDataRebinnedNaNPos = np.argwhere(np.isnan(pnSpecDataRebinned))
        print('getRadialVelocityFromXCor: vRad = ',vRad,': pnSpecDataRebinnedNaNPos = ',pnSpecDataRebinnedNaNPos)
        if pnSpecDataRebinnedNaNPos.shape[0] > 0:
            STOP

        pnSpecCompDataRebinnedNaNPos = np.argwhere(np.isnan(pnSpecCompDataRebinned))
        print('getRadialVelocityFromXCor: vRad = ',vRad,': pnSpecCompDataRebinnedNaNPos = ',pnSpecCompDataRebinnedNaNPos)
        if pnSpecCompDataRebinnedNaNPos.shape[0] > 0:
            for pos in pnSpecCompDataRebinnedNaNPos:
                pnSpecCompDataRebinned[pos] = 0.
                print('pnSpecCompDataRebinned[pos=',pos,'] = ',pnSpecCompDataRebinned[pos])

        chiSquare = np.sum(np.square(pnSpecDataRebinned - pnSpecCompDataRebinned)) / pnSpecCompDataRebinned.shape[0]
        print('chiSquare = ',chiSquare)
        if np.isnan(chiSquare):
            STOP
        chiSquares.append(chiSquare)
    print('vRadRange = ',vRadRange)
    plt.plot(vRadRange,chiSquares)
    plt.show()

    print('getRadialVelocityFromXCor: chiSquares = ',chiSquares)
    naNPos = np.argwhere(np.isnan(chiSquares))
    print('getRadialVelocityFromXCor: naNPos = ',naNPos)

    where = np.where(chiSquares==np.min(chiSquares))[0]
    print('where = ',where)
    range = int(len(vRadRange)/10.)
    popt = xCorFindMinimum(vRadRange[where[0]-range:where[0]+range], chiSquares[where[0]-range:where[0]+range], display = True)
    print('getRadialVelocityFromXCor: popt = ',popt)
    bestVRad = popt[1]
    if bestVRad - vRadRange[where[0]] > 1.5*np.abs(vRadRange[1]-vRadRange[0]):
        print('difference between bestVRad = ',bestVRad,' and vRadRange[',where[0],']=',vRadRange[where[0]],' gt dVRad')
        STOP
    print('getRadialVelocityFromXCor: bestVRad = ',bestVRad)
    wLenZero = applyVRadCorrection(pnSpecWLen,bestVRad)

    plt.plot(pnSpecCompWLen, pnSpecCompData, label = 'pnSpecComp')
    plt.plot(wLenZero, pnSpecData, label = 'pnSpec bestVRad')
    plt.legend()
    plt.show()
    return [bestVRad,popt,chiSquares]

def read_csv(filename):
    with open(filename) as f:
        file_data=csv.reader(f)
        headers=next(file_data)
        return [dict(zip(headers,i)) for i in file_data]

def getRadialVelocityFromSyntheticSpectrum(pnSpec):#, templateSpectraList):
    #templateSpectra = read_csv(templateSpectraList)
    #print('templateSpectra = ',templateSpectra)
    #print('dir(templateSpectra) = ',dir(templateSpectra))

    header = ['fitsFile','vrad','err_vrad']
    with open(pnSpec[:pnSpec.rfind('.')]+'_results.csv','w') as f:
        writer = csv.writer(f,delimiter=',')
        writer.writerow(header)
        templateSpectrum = makeTemplateSpec(pnSpec, display=False)
        vrad,popt,chiSquares = getRadialVelocityFromXCor(pnSpec, templateSpectrum)
        writer.writerow([pnSpec,
                        '%.1f' % (vrad),
                        '%.1f' % (popt[2])])

def plotSpec(fitsFileName):
    wLen = getWavelengthArr(fitsFileName)
    spec = getImageData(fitsFileName,0)
    plt.plot(wLen,spec)
    plt.show()

def getPNLines():
    lineList = [[4026.30,21],
                [4068.60,25],
                [4101.70,20],
                [4267.15,10],
                [4340.50,40],
                [4363.20,22],
                [4387.90,5],
                [4471.60,3],
                [4541.60,2],
                [4685.70,61],
                [4740.30,7],
                [4861.30,100],
                [4921.90,1],
                [4958.90,622],
                [5006.80,1790],
                [5198.50,15],
                [5411.50,5],
                [5517.20,5],
                [5537.70,6],
                [5754.80,18],
                [5875.80,12],
                [6300.20,33],
                [6310.20,16],
                [6363.90,10],
                [6548.10,418],
                [6562.80,464],
                [6583.60,1236],
                [6678.10,4],
                [6717.00,8],
                [6731.30,14],
                [7065.20,6],
                [7135.80,37],
                [7236.00,5],
                [7263.30,2],
                [7281.30,10],
                #[7325.00,23]
                ]
    return lineList

def synPNSpec():
    lineList = getPNLines()
    wLen = np.arange(4000.,7400.,0.5)
    spec = np.zeros(wLen.shape[0])
    for line in lineList:
        spec += gauss(wLen,line[1],line[0],4.)
    plt.plot(wLen,spec)
    plt.show()

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

def makeTemplateSpec(fitsFileName, display=False):
    lineList = getPNLines()

    wLen = getWavelengthArr(fitsFileName,0)
    spec = getImageData(fitsFileName,0)

    tempSpec = np.zeros(len(wLen))
    for line in lineList:
        if (line[0] > wLen[0]) & (line[0] < wLen[len(wLen)-1]):
            amp = np.max(spec[np.where((wLen > (line[0] - 4.3)) & (wLen < (line[0] + 4.3)))[0]])
            tempSpec += gauss(wLen,amp,line[0],4.)

    if display:
        plt.plot(wLen,spec,label='original')
        plt.plot(wLen,tempSpec,label='template')
        plt.legend()
        plt.show()

    fitsOutName = fitsFileName[:fitsFileName.rfind('.')]+'_template.fits'
    writeFits1D(tempSpec,
                fitsOutName,
                header=fitsFileName,
                CRVAL1=getHeaderValue(fitsFileName,'CRVAL1'),
                CDELT1=getHeaderValue(fitsFileName,'CDELT1'),
                CRPIX1=getHeaderValue(fitsFileName,'CRPIX1'),
                )
    return fitsOutName
