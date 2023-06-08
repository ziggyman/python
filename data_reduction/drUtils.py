#from astropy import units as u
#from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, Angle, ICRS, LSR
from astropy.nddata import StdDevUncertainty
import astropy.io.fits as pyfits
import astropy.modeling.tests.irafutil as iu
from astropy.nddata import CCDData
from astropy.time import Time
import ccdproc# import Combiner, subtract_overscan
from collections import namedtuple
import matplotlib.colorbar as cbar
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.legendre import legval
#from numpy.polynomial import Legendre as L
#import hammer
#import numpy as np
import os
from PyAstronomy import pyasl
#from pyraf import iraf
from scipy import exp,ndimage
from scipy.integrate import simps
#from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from scipy.interpolate import make_lsq_spline#, BSpline
#from scipy import interpolate
from scipy.interpolate import griddata
#from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit, minimize
from scipy.signal import find_peaks
#from scipy.signal import medfilt
from shutil import copyfile
#from sklearn import linear_model
from specutils import Spectrum1D
from specutils.manipulation.resample import FluxConservingResampler,LinearInterpolatedResampler
from specreduce import fluxcal,calibration_data
#from fluxcal import standard_sensfunc, apply_sensfunc, onedstd, obs_extinction, airmass_cor
from apextract import trace, extract
from fluxcal import standard_sensfunc, apply_sensfunc, onedstd, obs_extinction, airmass_cor
from myUtils import hmsToDeg,dmsToDeg,subtractSky#,sigmaRej
from myUtils import getDate,getDateTime

# TODO: Make it possible to pass in CCDData instead of fits file names

Info = namedtuple('Info', 'start height')
c0 = 299792.458 # km/s
plot = False
fluxStdDirsByPriority = ['spec50cal',
                         'oke1990',
                         'irscal',
                         'iidscal',
                         'spechayescal',
                         'gemini']

def readFileToArr(fname):
    with open(fname, "r") as f:
        lines = f.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

# remove file <fileName> (including directory) if it exists
# also works for </dir/start"*"ending>
def silentRemove(fileName):
    if '*' not in fileName:
        if os.path.exists(fileName): os.remove(fileName)
    else:
        dirName = fileName[:fileName.rfind('/')+1]
        fileList = os.listdir(dirName)
        for item in fileList:
            if item.endswith(fileName[fileName.find('*')+1:]):
                os.remove(os.path.join(dirName, item))

#def writeFits(ccdData, output_file, overwrite=True):
#    hdulist = ccdData.to_hdu()
#    hdulist.writeto(output_file, overwrite=overwrite)


def getImageData(fname,hduNum=1):
    hdulist = pyfits.open(fname)
    scidata = hdulist[hduNum].data
    hdulist.close()
    return scidata

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

def readMyPNLineList(fName = '/Users/azuri/entwicklung/python/data_reduction/pnLineList.dat'):
    with open(fName,'r') as f:
        lines = f.readlines()
    lineWave = []
    lineLabel = []
    for line in lines[1:]:
        label,wave = line.strip().split()
        lineWave.append(float(wave))
        lineLabel.append(label)
    return [lineWave,lineLabel]

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
    hdulist.writeto(outputFileName, overwrite=overwrite)
    print('writeFits: new mean after writing image = ',np.mean(getImageData(outputFileName,len(hdulist)-1)))

# insert <(_)suffix> before '.ending'
# if <suffix> is '' then no '_' is inserted either and the original <fileName>
# is returned
def addSuffixToFileName(fileName, suffix):
    if suffix == '':
        return fileName

    # suffix is not empty:
    fileNameOut = (fileName[:fileName.rfind('.')]
                 + '_'
                 + suffix
                 + fileName[fileName.rfind('.'):]
                 )
    return fileNameOut

# take a list of fits files called <inList> and create separate lists for
# * Biases
# * DomeFlats
# * SkyFlats
# * non Biases
# * non Flats
# * objects
# * each individual object
# all of these lists with the suffixes given in <suffixes>
#
# fits header keys read: * OBJECT
#                        * EXPTYPE
#
# parameters:
# <inList>: string
# <suffixes>: list of suffixes for which to create the individual file lists
#             #e.g. = ['','_z','_zf']
# <exptypes>: list of strings of exposure types for which to create lists,
#             default is None. If not None then <objects> MUST also be given
# <objects>: list of lists of strings for each exposure type for which to create
#            individual lists. If '*' then all files with the exposure type will
#            be added to the list. If 'individual' then one list will be created
#            for all individual objects with exposure type in <exptypes>.
#            Default is None. Must be given if <exptypes>
#            is not None
# <changeNames>: if True copy each <file> in <inList> to <EXPTYPE_OBJECT_><file>
def separateFileList(inList, suffixes, exptypes=None, objects=None, changeNames=False, fluxStandardNames=None):
    def createLists(exptype, object, lines, suffix):
        print('createLists: exptype = <'+exptype+'>, object = <'+object+'>, suffix = <'+suffix+'>')
        listMiddleName = object
        if (object == '*'):
            listMiddleName = ''
        isList = os.path.join(path, exptype+listMiddleName+'.list')
        isntList = os.path.join(path, 'non'+exptype+listMiddleName+'.list')
        isList = addSuffixToFileName(isList, suffix)
        isntList = addSuffixToFileName(isntList, suffix)
        isNames = []
        isntNames = []
        fOutName = ''
        for line in lines:
            print('createLists: reading '+line)
            hdulist = pyfits.open(line)

            fOutName = line
#            fOutName = addSuffixToFileName(fOutName, suffix)
            try:
                expType = hdulist[0].header['EXPTYPE']
            except:
                try:
                    expType = hdulist[0].header['IMAGETYP']
                except:
                    expType = hdulist[0].header['OBJECT']
            if ' ' in expType:
                expType = expType[:expType.find(' ')]
            objectName = hdulist[0].header['OBJECT']
            if ' ' in objectName:
                objectName = objectName[:objectName.find(' ')]
            if changeNames:
                fOutName = os.path.join(fOutName[:fOutName.rfind('/')],
                                        expType+'_'+objectName+'_'+fOutName[fOutName.rfind('/')+1:])
                if suffix == '':
                    copyfile(line, fOutName)
                line = fOutName
            if expType.lower() == exptype.lower():
                if object == '*':
                    print('createLists: exptype = <'+exptype+'>, object = <'+object+'>: expType = <'+expType+'>, objectName = <'+objectName+'>: adding <'+line+'> to isNames')
                    isNames.append(line)
                else:
                    if objectName.lower() == object.lower():
                        print('createLists: exptype = <'+exptype+'>, object = <'+object+'>: expType = <'+expType+'>, objectName = <'+objectName+'>: adding <'+line+'> to isNames')
                        isNames.append(line)
                    else:
                        print('createLists: exptype = <'+exptype+'>, object = <'+object+'>: expType = <'+expType+'>, objectName = <'+objectName+'>: adding <'+line+'> to isntNames')
                        isntNames.append(line)
            elif exptype.lower() == 'FLUXSTDS'.lower():
                if fluxStandardNames is None:
                    print('createLists: ERROR: no fluxStandardNames given')
                    STOP
                print('createLists: objectName = ',objectName)
                if objectName.lower() in fluxStandardNames:
                    if expType.lower():
                        isNames.append(line)
#                        print('createLists: isNames = ',isNames)
#                        STOP
                else:
                    isntNames.append(line)
            else:
                print('createLists: exptype = <'+exptype+'>, object = <'+object+'>: expType = <'+expType+'>, objectName = <'+objectName+'>: adding <'+line+'> to isntNames')
                isntNames.append(line)
        print('createLists: writing isList <'+isList+'>')
        with open(isList,'w') as f:
            for name in isNames:
                f.write(addSuffixToFileName(name,suffix)+'\n')
        print('createLists: writing isntList <'+isntList+'>')
        with open(isntList,'w') as f:
            for name in isntNames:
                f.write(addSuffixToFileName(name,suffix)+'\n')

    lines = []
    path = inList[:inList.rfind('/')]
    with open(inList,'r') as f:
        lines = [os.path.join(path,line.strip('\n')) for line in f]

    for suffix in suffixes:
        listOutName = ''
        listOutNameObs = ''

        objectNames = []
        expTypes = []

        if (exptypes is None) and (objects is None):
            nonZerosOutName = os.path.join(path,'nonZeros.list')
            nonZerosOutName = addSuffixToFileName(nonZerosOutName, suffix)

            nonFlatsOutName = os.path.join(path,'nonFlats.list')
            nonFlatsOutName = addSuffixToFileName(nonFlatsOutName, suffix)

            # delete existing output lists
            for fileName in [nonZerosOutName, nonFlatsOutName]:
                silentRemove(fileName)

            for line in lines:
                hdulist = pyfits.open(line)

                fOutName = line
                fOutName = addSuffixToFileName(fOutName, suffix)

                expType = hdulist[0].header['EXPTYPE']
                print('separateFileList: line = ',line,': expType = ',expType)
                if expType.lower() not in [x['expType'].lower() for x in expTypes]:
                    listOutName = os.path.join(path,expType.lower()+'.list')
                    listOutName = addSuffixToFileName(listOutName,suffix)
                    print('separateFileList: writing to listOutName = <'+listOutName+'>')
                    expTypes.append({'expType':expType, 'listOutName':listOutName})
                    silentRemove(listOutName)
                else:
                    for x in expTypes:
                        if x['expType'].lower() == expType.lower():
                            listOutName = x['listOutName']
                print('separateFileList: listOutName = <'+listOutName+'>')

                objectName = hdulist[0].header['OBJECT']
                print('separateFileList: line = ',line,': objectName = ',objectName,', expType = ',expType)
                if expType == 'SCIENCE':
                    if objectName not in [x['objectName'] for x in objectNames]:
                        listOutNameObs = os.path.join(path,objectName.lower()+'.list')
                        listOutNameObs = addSuffixToFileName(listOutNameObs, suffix)
                        print('separateFileList: writing to listOutNameObs = <'+listOutNameObs+'>')
                        objectNames.append({'objectName':objectName,'listOutNameObs':listOutNameObs})
                        silentRemove(listOutNameObs)
                    else:
                        for x in objectNames:
                            if x['objectName'] == objectName:
                                listOutNameObs = x['listOutNameObs']

                if changeNames:
                    fOutName = os.path.join(fOutName[:fOutName.rfind('/')],
                                            expType+'_'+objectName+'_'+fOutName[fOutName.rfind('/')+1:])
                    if suffix == '':
                        copyfile(line, fOutName)

                # write output lists
                with open(listOutName,'a') as f:
                    f.write(fOutName+'\n')

                if expType.lower() == 'SCIENCE'.lower():
                    with open(listOutNameObs,'a') as f:
                        f.write(fOutName+'\n')

                if objectName.lower() != 'Bias'.lower():
                    with open(nonZerosOutName,'a') as f:
                        f.write(fOutName+'\n')

                if expType.lower() not in ['FLAT'.lower(),'BIAS'.lower()]:
                    with open(nonFlatsOutName,'a') as f:
                        f.write(fOutName+'\n')
        else:
            individualLists = []
            if len(exptypes) != len(objects):
                print('separateFileList: ERROR: lengths of <exptypes> and <objects> are not the same')
                STOP
            for iExptype in np.arange(0,len(exptypes),1):
                exptype = exptypes[iExptype].lower()
                print('exptype = ',exptype)
                for object in objects[iExptype]:
                    if  object.lower() == 'individual':
                        for line in lines:
                            hdulist = pyfits.open(line)
                            try:
                                expType = hdulist[0].header['EXPTYPE'].lower()
                            except:
                                expType = hdulist[0].header['OBJECT'].lower()
                            if ' ' in expType:
                                expType = expType[:expType.find(' ')].lower()
                            if expType.lower() == exptype.lower():
                                objectName = hdulist[0].header['OBJECT']
                                if ' ' in objectName:
                                    objectName = objectName[:objectName.find(' ')]
                                if objectName not in individualLists:
                                    individualLists.append(objectName)
                    else:
                        print('creating lists for exptype=',exptype,', object=',object,', lines=',lines,', suffix = ',suffix)
                        createLists(exptype,object,lines,suffix)
                        #if exptype == 'science':
                        #    STOP

                if len(individualLists) > 0:
                    for object in individualLists:
                        if exptype.lower() != 'FLUXSTDS'.lower():
                            createLists(exptype,object,lines,suffix)
# combine images in <ccdImages>, write output to <fitsOutName> if not None
#
# parameters:
# <ccdImages>: [ccdimage1, ccdimage2,...]
# <combinerMethod>: string 'average' or 'median'
# <clippingMethod>: (None, 'minmax', 'sigma', 'extrema')
# <clippingParameters>: minmax: {'min_clip':-0.3, - clip all pixels with a value below -0.3
#                                'max_clip':0.1} - clip all pixels with a value above 0.1
#                       sigma: {'niter':0, - 0: iterate until no more pixels are clipped
#                                            n: iterate n times
#                               'low_thresh':2, - mask pixels more than 2 standard deviations below
#                               'high_thresh':5, - or 5 standard deviations above
#                               'func':np.ma.median} - the median
#                       extrema: {'nlow':1, - mask the lowest pixel value
#                                 'nhigh':2} - and the highest two pixel values
# <scaling>:     boolean - normalize each image by its mean before combining
# <minVal>: float - replace all values less than <minVal> with <minVal>
# <fitsOutName>: string, default is don't write (None)
# <overwrite>:   boolean, default is True
#
# TODO: Add reprojection to common WCS (see https://ccdproc.readthedocs.io/en/latest/ccdproc/image_combination.html
#                                       and https://reproject.readthedocs.io/en/stable/ )
def combine(ccdImages,
            combinerMethod='median',
            clippingMethod='None',
            clippingParameters=None,
            scaling=False,
            minVal=None,
            fitsOutName=None,
            overwrite=True):
    ccdData = [CCDData.read(fname, unit="adu") for fname in ccdImages]
    combiner = ccdproc.Combiner(ccdData)
    if clippingMethod == 'minmax':
        print('combine: applying minmax clipping')
        combiner.minmax_clipping(min_clip=clippingParameters['min_clip'],
                                 max_clip=clippingParameters['max_clip'])

    elif clippingMethod == 'sigma':
        print('combine: applying sigma clipping')
        #print('combine: dir(clippingParameters = ',dir(clippingParameters))
        print('combine: clippingParameters.keys() = ',clippingParameters.keys())
        print("combine: clippingParameters['low_thresh'] = ",clippingParameters['low_thresh'])
        print("combine: clippingParameters['high_thresh'] = ",clippingParameters['high_thresh'])
        combiner.sigma_clipping(low_thresh=clippingParameters['low_thresh'],
                                high_thresh=clippingParameters['high_thresh'],
                                func=clippingParameters['func'])
        iIter = 1
        if (clippingParameters['niter'] == 0) or (iIter < clippingParameters['niter']):
            old_n_masked = 0  # dummy value to make loop execute at least once
            new_n_masked = combiner.data_arr.mask.sum()
            while (new_n_masked > old_n_masked):
                combiner.sigma_clipping(func=np.ma.median)
                old_n_masked = new_n_masked
                new_n_masked = combiner.data_arr.mask.sum()
                if iIter == clippingParameters['niter']:
                    break
                iIter += 1

    elif clippingMethod == 'extrema':
        print('combine: applying extrema clipping')
        combiner.clip_extrema(nlow=clippingParameters['nlow'],
                              nhigh=clippingParameters['nhigh'])

    if scaling:
        print('combine: applying image scaling')
        scaling_func = lambda arr: 1/np.ma.average(arr)
        combiner.scaling = scaling_func

    combinedImage = None
    if combinerMethod == 'average':
        combinedImage = combiner.average_combine()
    elif combinerMethod == 'median':
        combinedImage = combiner.median_combine()

    # replace all values less than <minVal> with <minVal>
    if minVal is not None:
        combinedImage.data = np.maximum(combinedImage.data, minVal)

    if fitsOutName is not None:
        print('combine: combinedImage.data.shape = ',combinedImage.data.shape)
        print('combine: combinedImage.header = ',combinedImage.header)
        writeFits(combinedImage, ccdImages[0], fitsOutName, ['COMBINED'], ['from %d images' % len(ccdImages)], overwrite=overwrite)

    return combinedImage

# subtract overscan and trim the images in <fitsFilesIn>, write result to <fitsFilesOut>
#
# parameters:
# <fitsFilesIn>:     [fitsFile1, fitsFile2, ..., fitsFileN]
# <overscanSection>: '[91:100, :]' : uses all rows of columns 90 through 99 as the overscan
# <trimSection>:     '[:91, :]': remove columns 91 until the end for all rows, default is None
# <fitsFilesOut>:    [fitsFile1, fitsFile2, ..., fitsFileN], default is None for don't write
# <overwrite>:       boolean, default is True
def subtractOverscan(fitsFilesIn, overscanSection, trimSection=None, fitsFilesOut=None, overwrite=True):
    dataOut = []
    for iFile in np.arange(0,len(fitsFilesIn),1):
        print('subtractOverscan: reading file '+fitsFilesIn[iFile])
        ccdData = CCDData.read(fitsFilesIn[iFile], unit="adu")
        ccdDataNoOverscan = ccdproc.subtract_overscan(ccdData, fits_section=overscanSection)
        if trimSection is None:
            if fitsFilesOut is not None:
                ccdDataNoOverscan.write(fitsFilesOut[iFile], overwrite=overwrite)
            dataOut.append(ccdDataNoOverscan)
        else:
            trimmed = ccdproc.trim_image(ccdDataNoOverscan, fits_section=trimSection)
            if fitsFilesOut is not None:
                writeFits(trimmed, fitsFilesIn[iFile], fitsFilesOut[iFile], ['OVERSCAN'], ['subtracted'], overwrite=overwrite)
            dataOut.append(trimmed)
    return dataOut

# subtract <masterBias> from the images in <fitsFilesIn>, write results to <fitsFilesOut>
#
# parameters:
# <fitsFilesIn>:     [fitsFile1, fitsFile2, ..., fitsFileN]
# <masterBias>:      master Bias to subtract
# <fitsFilesOut>:    [fitsFile1, fitsFile2, ..., fitsFileN], default is None for don't write
# <overwrite>:       boolean, default is True
def subtractBias(fitsFilesIn,
                 masterBias,
                 fitsFilesOut=None,
                 overwrite=True):
    dataOut = []
    masterBiasArr = CCDData.read(masterBias, unit="adu")
    for iFile in np.arange(0,len(fitsFilesIn),1):
        ccdData = CCDData.read(fitsFilesIn[iFile], unit="adu")
        print('subtractBias: inFile = ',fitsFilesIn[iFile],', masterBias = ',masterBias)
        ccdDataBiasSubtracted = ccdproc.subtract_bias(ccdData,
                                                      masterBiasArr)
        dataOut.append(ccdDataBiasSubtracted)
        if fitsFilesOut is not None:
            writeFits(ccdDataBiasSubtracted, fitsFilesIn[iFile], fitsFilesOut[iFile], ['BIAS'], ['subtracted %s' % masterBias], overwrite=overwrite)
    return dataOut

# find and remove cosmic rays from the images in <fitsFilesIn>, write results to <fitsFilesOut>
#
# parameters:
# <fitsFilesIn>:      [ccdimage1, ccdimage2,...]
# <cosmicMethod>:     string 'lacosmic' or 'median'
# <cosmicParameters>: lacosmic: {'sigclip':4.5, float, optional - Laplacian-to-noise limit for cosmic ray detection.
#                                               Lower values will flag more pixels as cosmic rays. Default: 4.5.
#                                'sigfrac':0.3, float, optional - Fractional detection limit for neighboring pixels.
#                                               For cosmic ray neighbor pixels, a Laplacian-to-noise detection limit
#                                               of sigfrac * sigclip will be used. Default: 0.3.
#                                'objlim':5.0, float, optional - Minimum contrast between Laplacian image and the
#                                              fine structure image. Increase this value if cores of bright stars are
#                                              flagged as cosmic rays. Default: 5.0.
#                                'pssl':0.0, float, optional - Previously subtracted sky level in ADU. We always need
#                                            to work in electrons for cosmic ray detection, so we need to know the sky
#                                            level that has been subtracted so we can add it back in. Default: 0.0.
#                                'gain':1.0, float, optional - Gain of the image (electrons / ADU). We always need to
#                                            work in electrons for cosmic ray detection. Default: 1.0
#                                'readnoise':6.5, float, optional - Read noise of the image (electrons). Used to
#                                                 generate the noise model of the image. Default: 6.5.
#                                'satlevel':65535.0, float, optional - Saturation level of the image (electrons). This
#                                                    value is used to detect saturated stars and pixels at or above
#                                                    this level are added to the mask. Default: 65535.0.
#                                'niter':4, int, optional - Number of iterations of the LA Cosmic algorithm to perform.
#                                           Default: 4.
#                                'sepmed':True, bool, optional - Use the separable median filter instead of the full
#                                               median filter. The separable median is not identical to the full median
#                                               filter, but they are approximately the same and the separable median
#                                               filter is significantly faster and still detects cosmic rays well.
#                                               Default: True
#                                'cleantype':'meanmask', str, optional - Set which clean algorithm is used:
#                                                        "median": An unmasked 5x5 median filter.
#                                                        "medmask": A masked 5x5 median filter.
#                                                        "meanmask": A masked 5x5 mean filter.
#                                                        "idw": A masked 5x5 inverse distance weighted interpolation.
#                                                        Default: "meanmask".
#                                'fsmode':'median', str, optional - Method to build the fine structure image:
#                                                   "median": Use the median filter in the standard LA Cosmic algorithm.
#                                                   "convolve": Convolve the image with the psf kernel to calculate the
#                                                               fine structure image.
#                                                   Default: "median".
#                                'psfmodel':'gauss', str, optional - Model to use to generate the psf kernel if
#                                                    fsmode == 'convolve' and psfk is None. The current choices are
#                                                    Gaussian and Moffat profiles:
#                                                    "gauss" and "moffat" produce circular PSF kernels.
#                                                    The "gaussx" and "gaussy" produce Gaussian kernels in the x and y
#                                                    directions respectively.
#                                                    Default: "gauss".
#                                'psffwhm':2.5, float, optional - Full Width Half Maximum of the PSF to use to generate
#                                               the kernel. Default: 2.5.
#                                'psfsize':7, int, optional
#                                             Size of the kernel to calculate. Returned kernel will have size psfsize
#                                             x psfsize. psfsize should be odd. Default: 7.
#                                'psfk':None, numpy.ndarray (with float dtype) or None, optional
#                                             PSF kernel array to use for the fine structure image if fsmode == 'convolve'.
#                                             If None and fsmode == 'convolve', we calculate the psf kernel using psfmodel.
#                                             Default: None.
#                                'psfbeta':4.765, float, optional - Moffat beta parameter. Only used if fsmode=='convolve'
#                                                 and psfmodel=='moffat'. Default: 4.765.
#                                'verbose':False bool, optional - Print to the screen or not. Default: False.
#
#                     median: {'thresh':5., float, optional - Threshold for detecting cosmic rays. Default is 5.
#                              'error_image': numpy.ndarray, float or None, optional.
#                                             Error level. If None, the task will use the standard deviation
#                                             of the data. If an ndarray, it should have the same shape as data.
#                                             Default is None.
#                              'mbox':11, - int, optional. Median box for detecting cosmic rays. Default is 11.
#                              'gbox':0, - int, optional. Box size to grow cosmic rays. If zero, no growing
#                                           will be done. Default is 0.
#                              'rbox':0} - int, optional. Median box for calculating replacement values.
#                                          If zero, no pixels will be replaced. Default is 0
# <fitsFilesOut>:    ['ccdOut1.fits','ccdOut2.fits',...] default is don't write (None)
# <overwrite>:       boolean, default is True
def cleanCosmic(fitsFilesIn,
                cosmicMethod='lacosmic',
                cosmicParameters=None,
                fitsFilesOut=None,
                overwrite=True):
    ctype = cosmicMethod.lower().strip()
    ctypes = ['lacosmic', 'median']
    if not ctype in ctypes:
        raise('>>> Cosmic ray type "%s" NOT available [%s]' % (ctype, ' | '.join(ctypes)))
    dataOut = []
    if cosmicParameters is not None:
        if ctype == ctypes[0]:
            goodKeys = ['sigclip','sigfrac','objlim','pssl','gain','readnoise','satlevel',
                        'niter','sepmed','cleantype','fsmode','psfmodel','psffwhm','psfk',
                        'psfbeta','verbose']
        else:
            goodKeys = ['thresh','mbox','gbox','rbox','error_image']
        for key in cosmicParameters.keys():
            if key not in goodKeys:
                raise('key "'+key+'" not in available keys ',goodKeys)
    for iFile in np.arange(0,len(fitsFilesIn),1):
        ccdData = CCDData.read(fitsFilesIn[iFile], unit="adu")
        #print('cleanCosmic: dir(ccdData) = ',dir(ccdData))
        ccdDataArr = ccdData.data
        if ctype == 'lacosmic':
            sigclip = 4.5
            sigfrac = 0.3
            objlim = 5.0
            pssl = 0.0
            gain = 1.0
            readnoise = 6.5
            satlevel = 65535.0
            niter = 4
            sepmed = True
            cleantype = 'meanmask'
            fsmode = 'median'
            psfmodel = 'gauss'
            psffwhm = 2.5
            psfk = None
            psfbeta = 4.765
            verbose = False
            if cosmicParameters is not None:
                if 'cleantype' in cosmicParameters.keys():
                    cleantype = cosmicParameters['cleantype']
                if 'sigclip' in cosmicParameters.keys():
                    sigclip = cosmicParameters['sigclip']
                if 'sigfrac' in cosmicParameters.keys():
                    sigfrac = cosmicParameters['sigfrac']
                if 'objlim' in cosmicParameters.keys():
                    objlim = cosmicParameters['objlim']
                if 'pssl' in cosmicParameters.keys():
                    pssl = cosmicParameters['pssl']
                if 'gain' in cosmicParameters.keys():
                    gain = cosmicParameters['gain']
                if 'readnoise' in cosmicParameters.keys():
                    readnoise = cosmicParameters['readnoise']
                if 'satlevel' in cosmicParameters.keys():
                    satlevel = cosmicParameters['satlevel']
                if 'niter' in cosmicParameters.keys():
                    niter = cosmicParameters['niter']
                if 'sepmed' in cosmicParameters.keys():
                    sepmed = cosmicParameters['sepmed']
                if 'fsmode' in cosmicParameters.keys():
                    fsmode = cosmicParameters['fsmode']
                if 'psfmodel' in cosmicParameters.keys():
                    psfmodel = cosmicParameters['psfmodel']
                if 'psffwhm' in cosmicParameters.keys():
                    psffwhm = cosmicParameters['psffwhm']
                if 'psfk' in cosmicParameters.keys():
                    psfk = cosmicParameters['psfk']
                if 'psfbeta' in cosmicParameters.keys():
                    psfbeta = cosmicParameters['psfbeta']
                if 'verbose' in cosmicParameters.keys():
                    verbose = cosmicParameters['verbose']
            ccdDataCleaned = ccdproc.cosmicray_lacosmic(ccdDataArr,
                                                        sigclip=sigclip,
                                                        cleantype=cleantype,
                                                        sigfrac=sigfrac,
                                                        objlim=objlim,
                                                        pssl=pssl,
                                                        gain = gain,
                                                        readnoise = readnoise,
                                                        satlevel = satlevel,
                                                        niter = niter,
                                                        sepmed = sepmed,
                                                        fsmode = fsmode,
                                                        psfmodel = psfmodel,
                                                        psffwhm = psffwhm,
                                                        psfk = psfk,
                                                        psfbeta = psfbeta,
                                                        verbose = verbose)
        elif ctype == 'median':
            thresh = 5.
            mbox = 11
            gbox = 0
            rbox = 0
            error_image = None
#            print('cleanCosmic: ccdData.uncertainty = ',ccdData.uncertainty)
#            error_image = ccdData.uncertainty.array
#            print('cleanCosmic: error_image = ',error_image)
            if cosmicParameters is not None:
                if 'thresh' in cosmicParameters.keys():
                    thresh = cosmicParameters['thresh']
                if 'mbox' in cosmicParameters.keys():
                    mbox = cosmicParameters['mbox']
                if 'gbox' in cosmicParameters.keys():
                    gbox = cosmicParameters['gbox']
                if 'rbox' in cosmicParameters.keys():
                    rbox = cosmicParameters['rbox']
                if 'error_image' in cosmicParameters.keys():
                    error_image = cosmicParameters['error_image']
            ccdDataCleaned = ccdproc.cosmicray_median(ccdDataArr,
                                                      mbox=mbox,
                                                      rbox=rbox,
                                                      gbox=gbox,
                                                      thresh=thresh,
                                                      error_image=error_image)
#        ccdDataCleaned.header['COSMIC'] = ctype.upper()
        dataOut.append(ccdDataCleaned)
        ccdData.data = ccdDataCleaned
        if fitsFilesOut is not None:
            writeFits(ccdData, fitsFilesIn[iFile], fitsFilesOut[iFile], ['COSMICS'], [cosmicMethod], overwrite=overwrite)
    return dataOut


# axis: 0 (columns) or 1 (rows), or 2 for 2D box
# width: odd number
def boxCarMedianSmooth(imageData, axis, width):
    print('boxCarMedianSmooth: imageData.shape = ',imageData.shape)
    print('boxCarMedianSmooth: len(imageData.shape) = ',len(imageData.shape))
    print('boxCarMedianSmooth: imageData = ',imageData)
    newDataArray = None
    if len(imageData.shape) > 1:
        newDataArray = np.zeros(shape=imageData.shape, dtype=type(imageData[0,0]))
        if axis == 0:
            for iRow in range(imageData.shape[0]):
                for iCol in range(imageData.shape[1]):
                    if iCol < int(width/2.0):
                        iColStart = 0
                        iColEnd = iCol+int(width/2.0)+1
                    elif iCol > imageData.shape[1]-int(width/2.0)-1:
                        iColStart = iCol - int(width/2.0)
                        iColEnd = imageData.shape[1]
                    else:
                        iColStart = iCol - int(width/2.0)
                        iColEnd = iCol + int(width/2.0) + 1
                    newDataArray[iRow,iCol] = np.median(imageData[iRow,iColStart:iColEnd])
    #                print 'iRow = ',iRow,', iCol = ',iCol,': iColStart = ',iColStart,', iColEnd = ',iColEnd,': imageData[iRow,iColStart:iColEnd] = ',imageData[iRow,iColStart:iColEnd],': median = ',newDataArray[iRow,iCol]
        elif axis == 1:
            for iCol in range(imageData.shape[1]):
                for iRow in range(imageData.shape[0]):
                    if iRow < int(width/2.0):
                        iRowStart = 0
                        iRowEnd = iRow+int(width/2.0)+1
                    elif iRow > imageData.shape[0]-int(width/2.0)-1:
                        iRowStart = iRow - int(width/2.0)
                        iRowEnd = imageData.shape[1]
                    else:
                        iRowStart = iRow - int(width/2.0)
                        iRowEnd = iRow + int(width/2.0) + 1
                    newDataArray[iRow,iCol] = np.median(imageData[iRowStart:iRowEnd,iCol])
    #                print 'iCol = ',iCol,', iRow = ',iRow,': iRowStart = ',iRowStart,', iRowEnd = ',iRowEnd,': imageData[iRowStart:iRowEnd,iCol] = ',imageData[iRowStart:iRowEnd,iCol],': median = ',newDataArray[iRow,iCol]
        elif axis == 2:
            for iCol in range(imageData.shape[1]):
                if iCol < int(width/2.0):
                    iColStart = 0
                    iColEnd = iCol+int(width/2.0)+1
                elif iCol > imageData.shape[0]-int(width/2.0)-1:
                    iColStart = iCol - int(width/2.0)
                    iColEnd = imageData.shape[1]
                else:
                    iColStart = iCol - int(width/2.0)
                    iColEnd = iCol + int(width/2.0) + 1
                for iRow in range(imageData.shape[0]):
                    if iRow < int(width/2.0):
                        iRowStart = 0
                        iRowEnd = iRow+int(width/2.0)+1
                    elif iRow > imageData.shape[0]-int(width/2.0)-1:
                        iRowStart = iRow - int(width/2.0)
                        iRowEnd = imageData.shape[1]
                    else:
                        iRowStart = iRow - int(width/2.0)
                        iRowEnd = iRow + int(width/2.0) + 1
                    newDataArray[iRow,iCol] = np.median(imageData[iRowStart:iRowEnd,iColStart:iColEnd])
    #                print 'iColStart = ',iColEnd,', iRow = ',iRow,': iRowStart = ',iRowStart,', iRowEnd = ',iRowEnd,': imageData[iRowStart:iRowEnd,iColStart:iColEnd] = ',imageData[iRowStart:iRowEnd,iColStart:iColEnd],': median = ',newDataArray[iRow,iCol]
        else:
            print('boxCarMedianSmooth: ERROR: axis(=',axis,') out of bounds [0,1,2]')
    else:
        newDataArray = np.zeros(shape=imageData.shape, dtype=type(imageData[0]))
        for iCol in range(imageData.shape[0]):
            if iCol < int(width/2.0):
                iColStart = 0
                iColEnd = iCol+int(width/2.0)+1
            elif iCol > imageData.shape[0]-int(width/2.0)-1:
                iColStart = iCol - int(width/2.0)
                iColEnd = imageData.shape[0]
            else:
                iColStart = iCol - int(width/2.0)
                iColEnd = iCol + int(width/2.0) + 1
            newDataArray[iCol] = np.median(imageData[iColStart:iColEnd])
    return newDataArray

def imDivide(inFileNames, quotientImageName, outFileNames, zeroVal=0.):
    quotientArr = CCDData.read(quotientImageName, unit="adu").data
    for i in np.arange(0,len(inFileNames),1):
        data = CCDData.read(inFileNames[i], unit="adu")
        data.data = data.data / quotientArr
        data.data[np.where(quotientArr == 0.)] = zeroVal
        writeFits(data, inFileNames[i], outFileNames[i], overwrite=True)

def makeMasterFlat(combinedFlatIn, boxSize, minSNR, outFileNameMasterFlat=None, outFileNameMasterFlatSmoothed=None):
    ccdDataFlat = CCDData.read(combinedFlatIn, unit="adu")
    print('makeMasterFlat: mean(ccdDataFlat) = ',np.mean(ccdDataFlat.data))
    print('makeMasterFlat: ccdDataFlat.data = ',ccdDataFlat.data.shape,': ',ccdDataFlat.data)
    smoothedFlatArr = ndimage.median_filter(ccdDataFlat.data, boxSize)
    print('makeMasterFlat: mean(smoothedFlatArr) = ',np.mean(smoothedFlatArr))
    print('makeMasterFlat: smoothedFlatArr = ',smoothedFlatArr.shape,': ',smoothedFlatArr)
    masterFlatArr = smoothedFlatArr / ccdDataFlat.data
#    masterFlatArr = masterFlatArr / np.amax(masterFlatArr)
    print('makeMasterFlat: 1. mean(masterFlatArr) = ',np.mean(masterFlatArr))
    print('makeMasterFlat: 1. masterFlatArr = ',masterFlatArr.shape,': ',masterFlatArr)
    whereSNRltMinSNR = np.where(np.sqrt(ccdDataFlat.data) < minSNR)
    masterFlatArr[whereSNRltMinSNR] = 1.0
    print('makeMasterFlat: 2. mean(masterFlatArr) = ',np.mean(masterFlatArr))
    if outFileNameMasterFlat is not None:
        ccdDataFlat.data = masterFlatArr
        print('makeMasterFlat: mean(master ccdDataFlat) = ',np.mean(ccdDataFlat.data))
        writeFits(ccdDataFlat, combinedFlatIn, outFileNameMasterFlat, ['BOXSIZE','MINSNR'], ['%d' % boxSize, '%d' % minSNR], overwrite=True)
    if outFileNameMasterFlat is not None:
        ccdDataFlat.data = smoothedFlatArr
        writeFits(ccdDataFlat, combinedFlatIn, outFileNameMasterFlatSmoothed[0:outFileNameMasterFlatSmoothed.rfind('.')]+'_notScaled.fits', ['BOXSIZE','MINSNR'], ['%d' % boxSize, '%d' % int(minSNR)], overwrite=True)
        ccdDataFlat.data = smoothedFlatArr / np.amax(smoothedFlatArr)
        print('makeMasterFlat: mean(smoothed ccdDataFlat) = ',np.mean(ccdDataFlat.data))
        writeFits(ccdDataFlat, combinedFlatIn, outFileNameMasterFlatSmoothed, ['BOXSIZE','MINSNR'], ['%d' % boxSize, '%d' % int(minSNR)], overwrite=True)
    return masterFlatArr, smoothedFlatArr

# Correct the image for flat fielding.
# The flat field image is normalized by its mean or a user-supplied value before flat correcting.
#
# parameters:
# fitsFilesIn: ['ccdimage1.fits', 'ccdimage2.fits',...]
# flat:        string, name of fits file containing the combined Flatfield to apply to the data.
# min_value : float or None, optional
#             Minimum value for flat field. The value can either be None and no minimum value is
#             applied to the flat or specified by a float which will replace all values in the
#             flat by the min_value. Default is None.
# norm_value : float or None, optional
#              If not None, normalize flat field by this argument rather than the mean of the image.
#              This allows fixing several different flat fields to have the same scale. If this
#              value is negative or 0, a ValueError is raised. Default is None.
# add_keyword : str, Keyword or dict-like, optional
#               Item(s) to add to metadata of result. Set to False or None to completely disable
#               logging. Default is to add a dictionary with a single item: The key is the name
#               of this function and the value is a string containing the arguments the function
#               was called with, except the value of this argument.
# fitsFilesOut: ['ccdimageOut1.fits', 'ccdimageOut2.fits',...]
#               list of output fits files, must have same length as fitsFilesOut
#               Default is None.
# overwrite: boolean - Default is True
def flatCorrect(fitsFilesIn,
                flat,
                min_value=None,
                norm_value=None,
                add_keyword=True,
                fitsFilesOut=None,
                overwrite=True):
    ccdDataFlat = CCDData.read(flat, unit="adu")
    if ccdDataFlat.shape[0] == 2:
        ccdDataFlat.data = ccdDataFlat.data[0,:,:]
    print('flatCorrect: ccdDataFlat = ',ccdDataFlat)
    print('flatCorrect: min(ccdDataFlat) = ',np.min(ccdDataFlat),', max(ccdDataFlat) = ',np.max(ccdDataFlat))
    print('flatCorrect: ccdDataFlat[np.where(np.isnan(ccdDataFlat))] = ',ccdDataFlat[np.where(np.isnan(ccdDataFlat))])
    print('flatCorrect: ccdDataFlat[np.where(np.isinf(ccdDataFlat))] = ',ccdDataFlat[np.where(np.isinf(ccdDataFlat))])
    dataOut = []

    for iFile in np.arange(0,len(fitsFilesIn),1):
        ccdData = CCDData.read(fitsFilesIn[iFile], unit="adu")
        if ccdData.shape[0] == 2:
            ccdData.data = ccdData.data[0,:,:]
        print('flatCorrect: iFile = ',iFile,': ccdData = ',ccdData)
        print('flatCorrect: min(ccdData) = ',np.min(ccdData),', max(ccdData) = ',np.max(ccdData))
        print('flatCorrect: ccdData[np.where(np.isnan(ccdData))] = ',ccdData[np.where(np.isnan(ccdData))])
        print('flatCorrect: ccdData[np.where(np.isinf(ccdData))] = ',ccdData[np.where(np.isinf(ccdData))])
        ccdDataFlattened = ccdproc.flat_correct(ccdData,
                                                flat=ccdDataFlat,
                                                min_value=min_value,
                                                norm_value=norm_value,
                                                add_keyword=add_keyword)
        print('flatCorrect: iFile = ',iFile,' ',fitsFilesIn[iFile],': ccdDataFlattened = ',ccdDataFlattened)
        dataOut.append(ccdDataFlattened)
        if fitsFilesOut is not None:
            print('iFile = ',iFile,', len(fitsFilesIn) = ',len(fitsFilesIn),', len(fitsFilesOut) = ',len(fitsFilesOut))
            #print('fitsFilesIn = ',fitsFilesIn)
            #print('fitsFilesOut = ',fitsFilesOut)
            writeFits(ccdDataFlattened, fitsFilesIn[iFile], fitsFilesOut[iFile], ['FLATCORR'], [flat], overwrite=overwrite)
    return dataOut

# calculate Chebyshev polynomial for normalized x values [-1.0,...,1.0] and give coefficients
#
# parameters:
# xNorm : python array of normalized x values in the range [-1.,1.]
# coeffs: python array of Chebyshev polynomial coefficients
def chebyshev(xNorm, coeffs):
    if False:
        y = []
        z = [1.0]
        z.append(xNorm[0])
        for i in np.arange(2,len(coeffs)):
            z.append(2.0 * xNorm[0] * z[i-1] - z[i-2])
        print('chebyshev: z = ',z)
        y = []
        for x in xNorm:
            z[1] = x
            for i in np.arange(2,len(coeffs)):
                z[i] = (2.0 * x * z[i-1]) - z[i-2]
            yTemp = (coeffs[0] * z[0]) + (coeffs[1] * z[1])
            for i in np.arange(2,len(coeffs),1):
                yTemp += coeffs[i] * z[i]
            y.append(yTemp)
        print('chebyshev: xNorm = ',len(xNorm),': ',xNorm)
        print('chebyshev: y = ',len(y),': ',y)
        return y

    yCheck = chebval(xNorm, coeffs)
    print('chebyshev: yCheck = ',yCheck)
    return yCheck

# calculate Legendre polynomial for normalized x values [-1.0,...,1.0] and give coefficients
#
# parameters:
# xNorm : python array of normalized x values in the range [-1.,1.]
# coeffs: python array of Legendre polynomial coefficients
def legendre(xNorm, coeffs):
    if False:
        y = []
        z = [1.0]
        z.append(xNorm[0])
        for i in np.arange(2,len(coeffs)):
            z.append(((((2.0*(i+1.0))-3.0) * xNorm[0] * z[i-1]) - ((i-1.0) * z[i-2]) / (i)))
        print('legendre: z = ',z)
        y = []
        for x in xNorm:
            z[1] = x
            for i in np.arange(2,len(coeffs)):
                z[i] = ((((2.0*(i+1.0))-3.0) * x * z[i-1]) - ((i-1.0) * z[i-2])) / (i)
            print('legendre: z = ',z)
            yTemp = (coeffs[0] * z[0]) + (coeffs[1] * z[1])
            for i in np.arange(2,len(coeffs),1):
                yTemp += coeffs[i] * z[i]
                print('legendre: i = ',i,': yTemp = ',yTemp)
            y.append(yTemp)
        print('legendre: xNorm = ',len(xNorm),': ',xNorm)
        print('legendre: y = ',len(y),': ',y)
        return y

    yCheck = legval(xNorm,coeffs)
#    print('legendre: yCheck = ',yCheck)
    return yCheck

# calculate linear spline for normalized x values [-1.0,...,1.0] and give coefficients
#
# parameters:
# xRange : python array of length 2 giving the range for the x values
# order  : order of the linear spline (TODO: calculate from number of coefficients)
# coeffs : python array of linear spline coefficients
def linearSpline(xRange, order, coeffs):
    y = []
    print('linearSpline: order = ',order,', len(coeffs) = ',len(coeffs))
    for x in np.arange(xRange[0],xRange[1]+1):
        s = (x - xRange[0]) / (xRange[1] - xRange[0]) * order
        j = int(s)
        a = (j + 1) - s
        b = s - j
        print('linearSpline: x = ',x,': s = ',s,', j = ',j,', a = ',a,', b = ',b)
        y.append(coeffs[j] * a + coeffs[j+1] * b)
    print('linear spline: y = ',y)
    return y

# calculate cubic spline for normalized x values [-1.0,...,1.0] and give coefficients
#
# parameters:
# xRange : python array of length 2 giving the range for the x values
# order  : order of the linear spline (TODO: calculate from number of coefficients)
# coeffs : python array of cubic spline coefficients
def cubicSpline(xRange, order, coeffs):
    print('cubicSpline: order = ',order,', len(coeffs) = ',len(coeffs))
    y = []
    for x in np.arange(xRange[0],xRange[1]+1):
        s = (x - xRange[0]) / (xRange[1] - xRange[0]) * order
        j = int(s)
        a = (j + 1) - s
        b = s - j
        z_0 = a**3
        z_1 = 1.0 + 3.0 * a * (1.0 + a * b)
        z_2 = 1.0 + 3.0 * b * (1.0 + a * b)
        z_3 = b**3
        print('cubicSpline: x = ',x,': s = ',s,', j = ',j,', a = ',a,', b = ',b,', z_0 = ',z_0,', z_1 = ',z_1,', z_2 = ',z_2,', z_3 = ',z_3)
        y.append((coeffs[j] * z_0) + (coeffs[j+1] * z_1) + (coeffs[2+j] * z_2) + (coeffs[3+j] * z_3))
    print('cubic spline: y = ',y)
    return y

# calculate the trace function from an IRAF database file for a given aperture number
#
# parameters:
# dbFile: string, name of the database file including path
# apNum : int, aperture number (starting with 0)
# xRange: [1,size]
#
# return: x=[0...size-1], y=[0...size-1]
def calcTrace(dbFile, apNum=0, xRange = None, apOffsetX = 0.):
    imFileName = dbFile[:dbFile.find('database/')]+dbFile[dbFile.rfind('/ap')+3:]+'.fits'
    print('calcTrace: imFileName = <'+imFileName+'>')
    imageData = CCDData.read(imFileName, unit="adu")
    print('calcTrace: imageData.shape=',imageData.shape)
    records = iu.get_records(dbFile)
#    print('records[',apNum,'] = ',records[apNum])
#    print(dir(records[apNum]))
#    print('records[',apNum,'].fields = ',records[apNum].fields)
#    apdata = records[apNum].get_fields()['begin']
#    print('apdata = ',apdata)
#    STOP
    xCenter, yCenter = [float(c) for c in records[apNum].get_fields()['center'].split(' ')]
    xCenter += apOffsetX
    curve = records[apNum].get_fields()['curve']
    print('calcTrace: curve = ',curve)
    funcs = ['none','chebyshev','legendre','cubicSpline','linearSpline']
    function = funcs[int(curve[0])]
    function = records[apNum].get_fields()['function']
    print('calcTrace: function = ',function)
    #STOP
    order = curve[1][0]
    print('calcTrace: order = ',order)
    if xRange is None:
        xRange = [curve[2][0],curve[3][0]]
    print('calcTrace: xRange = ',xRange)
    xNorm = []
    xArr = np.arange(xRange[0],xRange[1]+1)
#    xArr = np.arange(1,imageData.shape[0]+1,1)#np.arange(xRange[0],xRange[1]+1)
    for x in xArr:
        if (function == funcs[1]) or (function == funcs[2]):
            xNorm.append((2.0 * x - (xRange[1] + xRange[0])) / (xRange[1] - xRange[0]))

    coeffs = [c[0] for c in curve[4:]]

    y = None
    if function == funcs[1]:
        print('function "chebyshev" detected')
        y = chebyshev(xNorm, coeffs) + xCenter
    elif function == funcs[2]:
        print('function "legendre" detected')
        y = legendre(xNorm, coeffs) + xCenter
    elif function == funcs[3]:
        y = cubicSpline(xRange, order, coeffs) + xCenter
    elif function == funcs[4]:
        y = linearSpline(xRange, order, coeffs) + xCenter
    else:
        print('calcTrace: could not identify function <'+function+'>')
    #print('y-1 = ',y-1)
    #STOP

    if len(xArr) != imageData.shape[0]:
        print('calcTrace: len(xArr)(=',len(xArr),') != imageData.shape[0](=',imageData.shape[0],': interpolating trace to get full image size')
        xArrFull = np.arange(1,imageData.shape[0]+1,1)
        f = interp1d(np.array(xArr), np.array(y), bounds_error = False,fill_value='extrapolate')
        y = f(np.array(xArrFull))
        xArr = xArrFull
    return [xArr-1.0,y-1.0]

# mark centers in fits file by setting the int(center) to zero
#
# parameters:
# imFileIn: string, name of fits file in which to mark the apertures
# trace: [xArr,yArr]
# imFileOut: name of output file if output file shall be created
#
# return:
# image with centers of aperture set to zero
def markCenter(imFileIn, trace, imFileOut=None):
    image = CCDData.read(imFileIn, unit="adu")
#    print('markCenter: image.shape = ',image.shape)
#    print('markCenter: trace = ',len(trace),': ',trace)
#    print('markCenter: trace[0].shape = ',trace[0].shape)
    for i in np.arange(0,trace[0].shape[0],1):
        image.data[int(trace[0][i]), int(trace[1][i])] = 0.
    if imFileOut is not None:
        writeFits(image, imFileIn, imFileOut, overwrite=True)
    return image

# returns height, width, and position of the top left corner of the largest
#  rectangle with the given value in mat
def max_size(mat, value=0):
    it = iter(mat)
    hist = [(el==value) for el in next(it, [])]
    max_size_start, start_row = max_rectangle_size(hist), 0
    for i, row in enumerate(it):
        hist = [(1+h) if el == value else 0 for h, el in zip(hist, row)]
        mss = max_rectangle_size(hist)
        if area(mss) > area(max_size_start):
            max_size_start, start_row = mss, i+2-mss[0]
    return max_size_start[:2], (start_row, max_size_start[2])

# returns height, width, and start column of the largest rectangle that
#  fits entirely under the histogram
def max_rectangle_size(histogram):
    stack = []
    top = lambda: stack[-1]
    max_size_start = (0, 0, 0) # height, width, start of the largest rectangle
    pos = 0 # current position in the histogram
    for pos, height in enumerate(histogram):
        start = pos # position where rectangle starts
        while True:
            if not stack or height > top().height:
                stack.append(Info(start, height)) # push
            elif stack and height < top().height:
                max_size_start = max(
                    max_size_start,
                    (top().height, pos - top().start, top().start),
                    key=area)
                start, _ = stack.pop()
                continue
            break # height == top().height goes here

    pos += 1
    for start, height in stack:
        max_size_start = max(max_size_start, (height, pos - start, start),
            key=area)

    return max_size_start

def area(size): return size[0]*size[1]

def transpose(inFileName, outFileName):
    image = np.array(CCDData.read(inFileName, unit="adu"))
    print('transpose: image.shape = ',image.shape)
    imageTransposed = image.transpose()
    print('transpose: imageTransposed.shape = ',imageTransposed.shape)
    writeFits(imageTransposed,
              inFileName,
              outFileName,
              metaKeys=['NAXIS1','NAXIS2',],
              metaData=[getHeaderValue(inFileName,'NAXIS2'),getHeaderValue(inFileName,'NAXIS1')],
              overwrite=True)
    image = np.array(CCDData.read(outFileName, unit="adu"))
    print('image.shape')
    plt.imshow(image, cmap='gray')
    plt.show()

# NOTE that the horizontal trace needs to come from an image that was rotated (IRAF rotate) by
# 90 degrees and flipped along the long axis
# WHICH IS equivalent to imtranspose
def interpolateTraceIm(imFiles, dbFileVerticalTrace, dbFileHorizontalTrace, markCenters=False):
    # read imFile for dimensions
    if len(imFiles) > 0:
        image = CCDData.read(imFiles[0], unit="adu")
        print('interpolateTraceIm: image.shape = ',image.shape)
        print('interpolateTraceIm: image.data.shape = ',image.data.shape)
        if image.data.shape[0] == 2:
            image.data = image.data[0,:,:]
        print('interpolateTraceIm: image.shape = ',image.shape)
        print('interpolateTraceIm: image.data.shape = ',image.data.shape)

        records = iu.get_records(dbFileVerticalTrace)
        verticalTraces = []
        for i in np.arange(0,len(records),1):
            verticalTraces.append(calcTrace(dbFileVerticalTrace,i,[1,image.data.shape[0]]))
    #        x, y = verticalTraces[len(verticalTraces)-1]

        records = iu.get_records(dbFileHorizontalTrace)
        horizontalTraces = []
        for i in np.arange(0,len(records),1):
            horizontalTraces.append(calcTrace(dbFileHorizontalTrace,i,[1,image.data.shape[1]]))
    #        x, y = horizontalTraces[len(horizontalTraces)-1]

        xyOrig = np.ndarray(shape=(image.data.shape[0] * image.data.shape[1],2), dtype=np.float32)
        xyFit = np.ndarray(shape=(image.data.shape[0] * image.data.shape[1],2), dtype=np.float32)
        zOrig = np.ndarray(shape=(image.data.shape[0] * image.data.shape[1]), dtype=np.float32)
        for ix in np.arange(0,image.data.shape[0],1):
            for iy in np.arange(0,image.data.shape[1],1):
                xyOrig[(ix*image.data.shape[1]) + iy,0] = ix
                xyOrig[(ix*image.data.shape[1]) + iy,1] = iy
                xyFit[(ix*image.data.shape[1]) + iy,0] = ix + horizontalTraces[0][1][iy] - horizontalTraces[0][1][0]
                xyFit[(ix*image.data.shape[1]) + iy,1] = iy + verticalTraces[0][1][ix] - verticalTraces[0][1][0]

        for imFile in imFiles:
            image = CCDData.read(imFile, unit="adu")
            if image.data.shape[0] == 2:
                image.data = image.data[0,:,:]

            if markCenters:
                inFile = ''
                outFile = ''
                for i in np.arange(0,len(verticalTraces),1):
                    if i == 0:
                        inFile = imFile
                        outFile = imFile[:imFile.rfind('.')]+'_vCenter%dMarked.fits' % (i)
                    else:
                        inFile = outFile
                        outFile = imFile[:imFile.rfind('.')]+'_vCenter%dMarked.fits' % (i)
                    markCenter(inFile,verticalTraces[i],outFile)
                horFile = outFile
                for i in np.arange(0,len(horizontalTraces),1):
                    if i == 0:
                        inFile = horFile
                        outFile = inFile[:inFile.rfind('.')]+'_hCenter%dMarked.fits' % (i)
                    else:
                        inFile = outFile
                        outFile = horFile[:horFile.rfind('.')]+'_hCenter%dMarked.fits' % (i)
                    markCenter(inFile,[horizontalTraces[i][1],horizontalTraces[i][0]],outFile)

            zOrig = np.ndarray(shape=(image.data.shape[0] * image.data.shape[1]), dtype=np.float32)
            print('interpolateTraceIm: image.data.shape = ',image.data.shape)
            for ix in np.arange(0,image.data.shape[0],1):
                for iy in np.arange(0,image.data.shape[1],1):
    #                print('interpolateTraceIm: image.data[',ix,',',iy,'] = ',image.data[ix,iy])
                    zOrig[(ix*image.data.shape[1]) + iy] = image.data[ix,iy]
            zFit = griddata(xyOrig, zOrig, xyFit, method='linear')

            # re-order vector back to 2D array
            for ix in np.arange(0,image.data.shape[0],1):
                for iy in np.arange(0,image.data.shape[1],1):
                    image.data[ix, iy] = zFit[(ix*image.data.shape[1]) + iy]

            #trim image to only contain good data inside the original trace
            maxSizeArr = np.zeros(image.data.shape)
            maxSizeArr[np.where(np.isnan(image.data))] = 1
            tempArr = max_size(maxSizeArr)
            image.data = image.data[tempArr[1][0]:tempArr[1][0]+tempArr[0][0],tempArr[1][1]:tempArr[1][1]+tempArr[0][1]]
    #        image.mask = image.mask[tempArr[1][0]:tempArr[1][0]+tempArr[0][0],tempArr[1][1]:tempArr[1][1]+tempArr[0][1]]
    #        image.uncertainty = image.uncertainty[tempArr[1][0]:tempArr[1][0]+tempArr[0][0],tempArr[1][1]:tempArr[1][1]+tempArr[0][1]]

            image.header['NAXIS1'] = tempArr[0][1]
            image.header['NAXIS2'] = tempArr[0][0]

            writeFits(image, imFile, imFile[:imFile.rfind('.')]+'i.fits', ['STRAIGHT'], ['interpolated'], overwrite=True)
        return 1
    else:
        return 0

def makeSkyFlat(skyFileIn, skyFlatOut, rowMedianSmoothSize = 7):
    image = CCDData.read(skyFileIn, unit="adu")
    print('makeSkyFlat: image[np.where(np.isnan(image))] = ',image[np.where(np.isnan(image))])
    print('makeSkyFlat: image[np.where(np.isinf(image))] = ',image[np.where(np.isinf(image))])
    print('makeSkyFlat: min(image) = ',np.min(image),', max(image) = ',np.max(image))
    profileImage = np.ndarray(image.data.shape, dtype=np.float32)

    # normalize each column of the sky flat
    for col in np.arange(0,profileImage.shape[1],1):
#        us = UnivariateSpline(np.arange(0,profileImage.shape[0],1),image.data[:,col], s=2. * float(profileImage.shape[0]))
#        profileImage[:,col] = us(np.arange(0,profileImage.shape[0],1))

        profileImage[:,col] = image[:,col] / np.median(image[:,col])
    print('makeSkyFlat: 1) profileImage[np.where(np.isnan(profileImage))] = ',profileImage[np.where(np.isnan(profileImage))])
    print('makeSkyFlat: 1) profileImage[np.where(np.isinf(profileImage))] = ',profileImage[np.where(np.isinf(profileImage))])

    # smooth each row
    for row in np.arange(0,profileImage.shape[0],1):
        profileImage[row,:] = ndimage.median_filter(profileImage[row,:], rowMedianSmoothSize)
    print('makeSkyFlat: profileImage = ',profileImage.shape,': ',profileImage)
    print('makeSkyFlat: 2) profileImage[np.where(np.isnan(profileImage))] = ',profileImage[np.where(np.isnan(profileImage))])
    print('makeSkyFlat: 2) profileImage[np.where(np.isinf(profileImage))] = ',profileImage[np.where(np.isinf(profileImage))])

    # set values <= to 0.01
    profileImage[np.where(profileImage <= 0.00001)] = 0.001
    print('makeSkyFlat: profileImage > 0. = ',profileImage.shape,': ',profileImage)
    print('makeSkyFlat: profileImage[np.where(np.isnan(profileImage))] = ',profileImage[np.where(np.isnan(profileImage))])
    print('makeSkyFlat: profileImage[np.where(np.isinf(profileImage))] = ',profileImage[np.where(np.isinf(profileImage))])

    image.data = profileImage
    writeFits(image, skyFileIn, skyFlatOut, ['SKY_FLAT'], ['rowMedianSmoothSize = %d' % rowMedianSmoothSize], overwrite=True)

# simple sum extraction
# imageFileIn: string or ndarray: str: name of fits file to extract; ndarray:image to extract
# dispAxis: string: [row, column]
def extractSum(imageFileIn, dispAxis, fNameOut = None):
    if isinstance(imageFileIn,str):
        image = CCDData.read(imageFileIn, unit="adu")
    else:
        image = imageFileIn
    print('extractSum: image.shape = ',image.shape)

    axis = 0
    if dispAxis == 'row':
        axis=0
    elif dispAxis == 'column':
        axis=1
    else:
        raise('dispAxis(=',dispAxis,' not valid, must be either "row" or "column"')
    flux = np.sum(image, axis = axis)
    if not fNameOut is None:
        writeFits1D(flux, fNameOut, wavelength=np.arange(1,len(flux)+1,1), header=imageFileIn, CRVAL1=None, CRPIX1=None, CDELT1=None)
    return flux

def lambdaCal(oneDImageFileIn, specOutName, func, coeffs):
    spec = np.array(CCDData.read(oneDImageFileIn, unit="adu"))
    print('lambdaCal: spec.shape = ',spec.shape)
    xSpec = range(spec.shape[0])
    xSpecNorm = normalizeX(xSpec)
    wLenSpec = func(xSpecNorm, coeffs)
    with open(specOutName,'w') as f:
        for i in xSpec:
            f.write('%.5f %d' % (wLenSpec[i], spec[i]))

# @brief: calculate the emission line profile for aperture number 'apNumber'
#         as provided in database/ap<twoDImageFileIn>
# @param twoDImageFileIn: string: name of input fits file (ARC)
# @param apNumber: int: number of aperture defined in database/ap<twoDImageFileIn>, starting with 0
# @param halfWidth: int: half width of emission line
# @param dxFit: float: dx for output fitted profile
# @param plot: bool: plot debugging plots?
def calcLineProfile(twoDImageFileIn,
                    apNumber,
                    halfWidth,
                    dxFit=0.01,
                    plot=False,
                    apOffsetX = 0.,
                    markCenter=False):
    image = np.array(CCDData.read(twoDImageFileIn, unit="adu"))
    print('calcLineProfile: image.shape = ',image.shape)

    centerRowIdx = int(image.shape[0]/2)

    tempFile = os.path.join(twoDImageFileIn[0:twoDImageFileIn.rfind('/')],'tmp'+twoDImageFileIn[twoDImageFileIn.rfind('/')+1:])
    print('calcLineProfile: tempFile = <'+tempFile+'>')

    copyfile(twoDImageFileIn,tempFile)

#    iraf.noao()
#    iraf.noao.twodspec()
#    iraf.noao.twodspec.apextract()
#    iraf.noao.twodspec.apextract.aptrace(input=tempFile,
#                                         apertures="",
#                                         references="",
#                                         interactive="yes",
#                                         find="no",
#                                         recenter='no',
#                                         resize='no',
#                                         edit='yes',
#                                         trace='yes',
#                                         fittrace='yes',
#                                         line=centerRowIdx,
#                                         nsum=3,
#                                         step=1,
#                                         nlost=10,
#                                         function='legendre',
#                                         order=5,
#                                         sample='*',
#                                         naverage=1,
#                                         niterate=1,
#                                         low_reject=3.0,
#                                         high_reject=3.0,
#                                         grow=0.0)


    dbFile = os.path.join(tempFile[:tempFile.rfind('/')],'database')
    dbFile = os.path.join(dbFile,'ap'+tempFile[tempFile.rfind('/')+1:tempFile.rfind('.')])
    print('calcLineProfile: dbFile = <'+dbFile+'>')
    row,center = calcTrace(dbFile, apNum=apNumber, xRange = None, apOffsetX=apOffsetX)
    if markCenter:
        markCenter(tempFile, [row, center], tempFile)
    #print('calcLineProfile: row = ',len(row),': ',row)
    #print('calcLineProfile: center = ',len(center),': ',center)

    colNumber = int(center[int(len(row)/2)])

    centerRow = image[centerRowIdx,:]
    print('calcLineProfile: twoDImageFileIn = ',twoDImageFileIn,': colNumber = ',colNumber,', halfWidth = ',halfWidth)
    print('calcLineProfile: centerRow[colNumber-halfWidth:colNumber+halfWidth+1] = ',centerRow[colNumber-halfWidth:colNumber+halfWidth+1])
    print('calcLineProfile: np.amax(centerRow[colNumber-halfWidth:colNumber+halfWidth+1]) = ',np.amax(centerRow[colNumber-halfWidth:colNumber+halfWidth+1]))
    print('calcLineProfile: centerRow = ',centerRow)
    maxPos = np.where(centerRow == np.amax(centerRow[colNumber-halfWidth:colNumber+halfWidth+1]))[0][0]
    print('calcLineProfile: twoDImageFileIn = ',twoDImageFileIn,': maxPos = ',maxPos)

    center += maxPos - center[centerRowIdx]
    #print('calcLineProfile: center = ',center)

    if False:#'PNG'in twoDImageFileIn:
        plot = True
    else:
        plot = False
    if plot:
        plt.plot(row,center)
        plt.show()
    if markCenter:
        markCenter(tempFile, [row,center], imFileOut=tempFile[:-5]+'_centerMarked'+str(apNumber)+'.fits')

    if plot:
        plt.plot(centerRow)
        plt.plot([colNumber-halfWidth,colNumber+halfWidth],[0.,0.])
        plt.show()
        #STOP

    profileDataX = []
    profileDataY = []
    xProfInt = np.arange(-halfWidth,halfWidth+1,1)
    print('calcLineProfile: xProfInt = ',xProfInt)
    rectangles = []
    intensities = []
    print('calcLineProfile: image.shape = ',image.shape)
    for row in np.arange(0,image.shape[0],1):
        image[row,xProfInt+int(center[row])] = image[row,xProfInt+int(center[row])] / np.sum(image[row,xProfInt+int(center[row])])
        for x in xProfInt:
#            print('calcLineProfile: x = ',x,': center[',row,',] = ',center[row],', int(center[row]) = ',int(center[row]))
#            print('calcLineProfile: x + 0.5 - center[row] + int(center[row]) = ',x + 0.5 - center[row] + int(center[row]))
            profileDataX.append(x + 0.5 - center[row] + int(center[row]))
            profileDataY.append(image[row,x+int(center[row])])

#            rectangles.append((x, row))
            rectangles.append((x+int(center[row]), row))
            intensities.append(image[row,x+int(center[row])])
#        print('calcLineProfile: np.sum(image[row,xProfInt+int(center[row])]) = ',np.sum(image[row,xProfInt+int(center[row])]))
    sortedIndices = np.argsort(profileDataX)
    #print('calcLineProfile: sortedIndices = ',sortedIndices)
    profileDataX = np.array(profileDataX)[sortedIndices]
    profileDataY = np.array(profileDataY)[sortedIndices]
    if plot:
        plt.scatter(profileDataX,profileDataY)

    if xProfInt[0] < profileDataX[0]:
        xProfInt = xProfInt[1:]
    if xProfInt[len(xProfInt)-1] > profileDataX[len(profileDataX)-1]:
        xProfInt = xProfInt[:-1]

    # interpolate Cubic Spline
    #print(xProfInt.shape, profileDataX.shape, profileDataY.shape)
    #print('calcLineProfile: profileDataX = ',profileDataX)
    #print('calcLineProfile: profileDataY = ',profileDataY)
    #cs = CubicSpline(profileDataX,profileDataY,bc_type=((1,0),(1,0)))
    #print('calcLineProfile: profileDataX[0] = ',profileDataX[0])
    #print('profileDataX[profileDataX.shape[0]-1]+dxFit = ',profileDataX[profileDataX.shape[0]-1]+dxFit)
    x = np.arange(profileDataX[0],profileDataX[profileDataX.shape[0]-1]+dxFit, dxFit)
    #print('calcLineProfile: x = ',x)

    #plt.plot(x,cs(x),'b-', label='Cubic Spline')

    #interpolate lsq spline
    t = xProfInt#[-1,-0.5,0,0.5,1]
    k = 3
    #print('calcLineProfile: (profileDataX[0],)*(k+1) = ',(profileDataX[0],)*(k+1))
    #print('calcLineProfile: t=',t)
    #print('calcLineProfile: (profileDataX[-1],)*(k+1) = ',(profileDataX[-1],)*(k+1))
    t = np.r_[(profileDataX[0],)*(k+1),
              t,
              (profileDataX[-1],)*(k+1)]
    #print('calcLineProfile: t=',t)
    #print('calcLineProfile: profileDataX = ',profileDataX.shape,': ',profileDataX)
    #for ind in range(profileDataX.shape[0]):
    #    print('profileDataX[',ind,'] = ',profileDataX[ind])
    #print('calcLineProfile: profileDataY = ',profileDataY.shape,': ',profileDataY)
    #print(profileDataX.ndim)
    #print(np.any(profileDataX[1:] <= profileDataX[:-1]))
    spl = make_lsq_spline(profileDataX, profileDataY, t, k)
    yFit = spl(x)
    if plot:
        plt.plot(x, yFit, 'y-', lw=3, label='LSQ spline')

        plt.show()


#    print('calcLineProfile: rectangles = ',rectangles)
#    print('calcLineProfile: intensities = ',intensities)
    normal = plt.Normalize(np.array(intensities).min(), np.array(intensities).max())
#    print('calcLineProfile: normal = ',type(normal),': ',normal)
    colors = plt.cm.Greys(normal(intensities))
#    print('calcLineProfile: colors = ',type(colors),': ',colors)
    cmap=plt.cm.Greys
#    print('calcLineProfile: cmap = ',type(cmap),': ',cmap)
    c = cmap((np.array(colors) - np.amin(colors))/(np.amax(colors) - np.amin(colors)))
#    print('calcLineProfile: c = ',type(c),': ',c)

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for iRectangle in np.arange(0,len(rectangles),1):
            rect = patches.Rectangle(rectangles[iRectangle], 1., 1., color=colors[iRectangle])
            ax.add_patch(rect)
        cax, _ = cbar.make_axes(ax)
        cb2 = cbar.ColorbarBase(cax, cmap=cmap,norm=normal)
        ax.set_xlim(np.amin(center)-halfWidth,np.amax(center)+halfWidth)
        ax.set_ylim(0,image.shape[0])
        ax.plot(center, np.arange(0,image.shape[0],1))
        plt.show()

#    plt.imshow(image)
#    plt.show()

    # cut off last half pixel on both sides as they are sometimes bad
#    print('calcLineProfile: x = ',x)
#    print('calcLineProfile: 0.5/dxFit = ',0.5/dxFit)
    x = x[int(0.5/dxFit):x.shape[0]-int(0.5/dxFit)]
    yFit = yFit[int(0.5/dxFit):yFit.shape[0]-int(0.5/dxFit)]
#    print('calcLineProfile: x trimmed = ',x)

    # subtract background and re-normalize to an integral of 1
    yFit = yFit - np.amin(yFit)
    yFit = yFit / simps(yFit,dx=x[1]-x[0])

    # center to maximum
    maxPos = np.where(yFit == np.amax(yFit))
    print('calcLineProfile: maxPos = ',maxPos)
    print('calcLineProfile: x[maxPos] = ',x[maxPos[0]])
    x = x-x[maxPos[0]]

    return [x, yFit]

# @brief: return x and y inside xRange
def getInsideRange(x,y,xRange):
    xInsideRange = x[np.where(x >= xRange[0])]
    yInsideRange = y[np.where(x >= xRange[0])]
    yInsideRange = yInsideRange[np.where(xInsideRange <= xRange[1])]
    xInsideRange = xInsideRange[np.where(xInsideRange <= xRange[1])]
    return xInsideRange, yInsideRange

# @brief: Return interpolated y at xAt
def getYAt(x,y,xAt):
#    print('getYAt: x = ',x.shape,': ',x)
#    print('getYAt: y = ',y.shape,': ',y)
#    print('getYAt: xAt = ',xAt.shape,': ',xAt)
#    xArr = np.array(x)
#    yArr = np.array(y)
#    print('getYAt: len(x) = ',len(x),', len(y) = ',len(y))
    f = interp1d(np.array(x), np.array(y), bounds_error = False,fill_value=0.)
    yOut = f(np.array(xAt))
#    yOut = []
#    if xArr.shape[0] != yArr.shape[0]:
#        raise Exception('getYAt(x,y,xAt): ERROR: x and y must have same size')
#    for xXAt in xAt:
#        if xXAt < x[0]:
#            yOut.append(yArr[0])
#            print('getYAt: xXAt = ',xXAt,' < x[0]=',x[0],': appended ',yOut[len(yOut)-1],' to yOut')
#        elif xXAt > xArr[xArr.shape[0]-1]:
#            yOut.append(yArr[xArr.shape[0]-1])
#            print('getYAt: xXAt = ',xXAt,' > xArr[',xArr.shape[0]-1,']=',xArr[xArr.shape[0]-1],': appended ',yOut[len(yOut)-1],' to yOut')
#        else:
#            for iX in np.arange(1,xArr.shape[0],1):
#                if (xArr[iX-1] <= xXAt) and (xArr[iX] >= xXAt):
#                    yOut.append(yArr[iX-1] + ((yArr[iX]-yArr[iX-1]) * (xXAt - xArr[iX])/(xArr[iX]-xArr[iX-1])))
#                    print('getYAt: xXAt = ',xXAt,': appended ',yOut[len(yOut)-1],' to yOut')
#    yOut = np.array(yOut)
#    if np.array(xAt).shape[0] != yOut.shape[0]:
#        raise Exception('getYAt(x,y,xAt): ERROR: xAt.shape[0](=',xAt.shape[0],'] != yOut.shape[0](=',yOut.shape[0],']')
    return yOut

# @brief: cross-correlate 2D arrays static and moving
# @param static: 2D array not moving [x,y]
# @param moving: 2D moving array, must be smaller than static [x,y], x=0 in center
def xCor(static, moving):
    xCorChiSquares = []
    dxMoving = moving[0][1]-moving[0][0]
#    print('xCor: dx = ',dxMoving)
#    print('xCor: xMoving = ',moving[0])
#    print('xCor: xStatic = ',static[0])
    xMovingStart = 0. - moving[0][0]#(moving[0][moving[0].shape[0]-1]-moving[0][0])/2.
#    print('xCor: xMovingStart = ',xMovingStart)
#    print('xCor: xMoving at half = ',moving[0][int(moving[0].shape[0]/2)])
    xMovingEnd = static[0][static[0].shape[0]-1]-moving[0][moving[0].shape[0]-1]
#    print('xCor: xMovingEnd = ',xMovingEnd)
    xXCor = np.arange(xMovingStart,xMovingEnd,dxMoving)
#    print('xCor: xXCor = [',xXCor[0],',',xXCor[1],',...,',xXCor[xXCor.shape[0]-1],']')
    for iX in np.arange(0,xXCor.shape[0],1):
        xMovingPlot = moving[0]+xXCor[iX]
#        print('xCor: xMovingPlot = ',len(xMovingPlot),': ',xMovingPlot)
#        plt.plot(xMovingPlot,moving[1])
        xStaticPlot,yStaticPlot = getInsideRange(static[0],static[1], [xMovingPlot[0],xMovingPlot[xMovingPlot.shape[0]-1]])
        yStaticPlot = yStaticPlot - np.amin(yStaticPlot)
#        print('xCor: xStaticPlot = ',xStaticPlot)
#        print('xCor: yStaticPlot = ',yStaticPlot)
#        plt.plot(xStaticPlot,yStaticPlot/simps(yStaticPlot,xStaticPlot))
#        plt.xlim(xMovingPlot[0],xMovingPlot[xMovingPlot.shape[0]-1])
        yAt = getYAt(xMovingPlot,moving[1],xStaticPlot)
        yAt = yAt * np.amax(yStaticPlot) / np.amax(yAt)
        #yStaticPlot = yStaticPlot / np.amax(yStaticPlot)
#        plt.plot(xStaticPlot,yAt)
#        plt.show()
        xCorChiSquares.append(np.sum(np.square(yStaticPlot - yAt)) / yAt.shape[0])
#    print('xCor: xXCor = ',xXCor.shape,': ',xXCor)
#    print('xCor: static[0] = ',static[0].shape,': ',static[0])
    yAtXCor = getYAt(xXCor,xCorChiSquares,static[0])

#    print('xCor: xCorChiSquares = ',xCorChiSquares)
    if plot:
        plt.plot(xXCor,xCorChiSquares/np.amax(xCorChiSquares), label='xCor')
        plt.plot(static[0],static[1]/np.amax(static[1]), label='static')
#        print('xCor: static[1].shape = ',static[1].shape)
#        print('xCor: yAtXCor.shape = ',yAtXCor.shape)
#        print('xCor: yAtXCor = ',yAtXCor)
        yPlot = yAtXCor / static[1]
        yPlot = yPlot / np.amax(yPlot)
        plt.plot(static[0], yPlot, label='xCor/static')
        plt.legend()
        plt.show()
    return np.array(xXCor), np.array(xCorChiSquares)

def gauss(x,a,x0,sigma,yBackground=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground

def xCorFindMinimum(xCorX, xCorY):
    y = xCorY - np.amax(xCorY)
    a = np.amin(y)
    print('xCorFindMinimum: a = ',a)
    x0 = xCorX[np.where(y == np.amin(y))][0]
    print('xCorFindMinimum: x0 = ',x0)
    popt,pcov = curve_fit(gauss,xCorX,y,p0=[a,x0,1.,0.])
    return popt[1]


# @brief: fit Gaussians every <step> pixels and identify line positions
# @param spec: 1D np.array:
# @param xCorX: 1D np.array: x values (output from xCor)
# @param xCorY: 1D np.array: y values (output from xCor)
# @param sigma: sigma of Gaussians to fit
# @param peakHeight: float: minimum peak height in spec for find_peaks
# @param peakWidth: float: minimum peak width in spec for find_peaks
def findLines(spec,xCorX,xCorY,sigma,peakHeight=None, peakWidth=None, threshold=None, plot=False):
    maxPosDiff = 0.67
    maxSigDiff = sigma * 0.6
    print('findLines: peakHeight = ',peakHeight,', peakWidth = ',peakWidth,', threshold = ',threshold)
    peaks,properties = find_peaks(spec, height = peakHeight, width=peakWidth, threshold=threshold)#
    if plot:
        plt.plot(spec)
        plt.scatter(peaks,spec[peaks])
        plt.title('peaks')
        plt.show()

    """Check that all peaks are within xCorX range"""
    for i in np.arange(len(peaks)-1,-1,-1):
        if (peaks[i] < xCorX[0]) or (peaks[i] > xCorX[len(xCorX)-1]):
            peaks = np.delete(peaks,i)
    print('findLines: peaks = ',peaks)
    print('findLines: properties = ',properties)
    print('findLines: spec.shape = ',spec.shape)
    print('findLines: xCorX.shape = ',xCorX.shape)
    print('findLines: xCorX = ',xCorX)

    yNorm = xCorY / np.amax(xCorY)
#    peaks,properties = find_peaks(0.-yNorm)#, height = -0.05, width=[30,120], threshold=0.00001)
#    print('findLines: peaks = ',peaks)
    if plot:
        plt.plot(xCorX,yNorm)
    xCorPeaks = []
    for peak in peaks:
        print('xCorX - ',peak,' = ',xCorX-peak)
        xCorPeaks.append(np.where(np.absolute(xCorX - peak) < (xCorX[1]-xCorX[0])/2.)[0])
        print('findLines: peak = ',peak,': xCorPeaks[',len(xCorPeaks)-1,'] = ',xCorPeaks[len(xCorPeaks)-1])
    xCorPeaks = np.array(xCorPeaks)
    if plot:
        plt.scatter(xCorX[xCorPeaks],yNorm[xCorPeaks])
        plt.title('scaled peak profile centers')
        plt.show()

    xDiff = []
    yDiff = []
    xYFitParams = []
    if plot:
        plt.plot(xCorX,yNorm)
        plt.scatter(xCorX[xCorPeaks],yNorm[xCorPeaks])
    for i in np.arange(0,xCorPeaks.shape[0],1):
        print('findLines: xCorX.shape = ',xCorX.shape)
        print('findLines: xCorPeaks.shape = ',xCorPeaks.shape)
        print('findLines: xCorPeaks[',i,'] = ',xCorPeaks[i])
        print('findLines: sigma = ',sigma)
        indices = np.where((xCorX > (xCorX[xCorPeaks[i]]-sigma)) & (xCorX < (xCorX[xCorPeaks[i]]+sigma)))
#            print('findLines: indices = ',indices)
        xi = xCorX[indices]
#            print('findLines: xi = ',xi)
        yi = yNorm[indices]
#            print('findLines: yi = ',yi)
        xCenter = xi[int(xi.shape[0]/2)]
        print('findLines: xCenter = ',xCenter,', xCorX[xCorPeaks[',i,']] = ',xCorX[xCorPeaks[i]])
        try:
            popt,pcov = curve_fit(gauss,xi,yi,p0=[-np.amax(yi),xCenter,sigma,np.amax(yi)])
        except:
            continue
        print('findLines: popt = ',popt)
        if popt[2] < 0.:
            popt[2] = 0.-popt[2]
        maxAmp = -0.00001
        if ((popt[0] < maxAmp)
            and (np.absolute(xCenter - popt[1]) < maxPosDiff)
            and (np.absolute(sigma - popt[2]) < maxSigDiff)
           ):
            xDiff.append(xi[int(xi.shape[0]/2)])
            yFit = gauss(xi,*popt)
            if plot:
                plt.plot(xi,yi)
                plt.plot(xi,yFit)
                plt.scatter(xCorX[xCorPeaks[i]],yNorm[xCorPeaks[i]])
            yDiff.append(np.sum(((yi-yFit) / np.amax(yi))**2) / yi.shape[0])
            xYFitParams.append([xi,yi,yFit,popt,i])
        else:
            print('findLines: rejected fit for line at ',xCorX[xCorPeaks[i]],', fitted parameters = [a=',popt[0],', x0=',popt[1],', sigma=',popt[2],', background=',popt[3],']')
            if popt[0] >= maxAmp:
                print('findLines: amplitude >= ',maxAmp)
            if np.absolute(xCenter - popt[1]) >= maxPosDiff:
                print('findLines: np.absolute(xCenter - popt[1])(=',np.absolute(xCenter - popt[1]),') >= maxPosDiff=',maxPosDiff)
            if np.absolute(sigma - popt[2]) >= maxSigDiff:
                print('findLines: np.absolute(sigma - popt[2])(=',np.absolute(sigma - popt[2]),') >= maxSigDiff=',maxSigDiff)
            xDiff.append(xi[int(xi.shape[0]/2)])
            yFit = gauss(xi,*popt)
            if plot:
                plt.plot(xi,yi)
                plt.plot(xi,yFit)
                plt.scatter(xCorX[xCorPeaks[i]],yNorm[xCorPeaks[i]])
                plt.title('rejected line at %.2f' % (xCenter))
                plt.show()
            yDiff.append(np.sum(((yi-yFit) / np.amax(yi))**2) / yi.shape[0])
    if plot:
        plt.show()

        plt.plot(xDiff,yDiff)
        plt.show()
    sigmas = [par[3][2] for par in xYFitParams]
    print('findLines: good sigmas = ',sigmas)

    if plot:
        plt.plot(spec/np.amax(spec))
        plt.plot(xCorX,yNorm)
        plt.scatter(xCorX[xCorPeaks],yNorm[xCorPeaks])
        for iGoodFit in np.arange(0,len(xYFitParams),1):
            plt.plot(xYFitParams[iGoodFit][0],xYFitParams[iGoodFit][1])
            plt.plot(xYFitParams[iGoodFit][0],xYFitParams[iGoodFit][2])
            #plt.scatter(xCorX[xCorPeaks[xYFitParams[iGoodFit][3]]],yNorm[xCorPeaks[xYFitParams[iGoodFit][3]]])
        plt.show()

    return [par[3][1] for par in xYFitParams]

# NOTE: currently it is required for the x values to be accending
def normalizeX(x):
    xZero = x - x[0]
    return 2.0 * xZero / xZero[-1] - 1.0

# @brief : fit background and subtract from y
# @param x : 1D array of x-values
# @param y : 1D array of y-values
# @param deg : int: degree for Legendre Polynomial
# @param indicesToIgnore : 1D array of indices to ignore for fitting the background
def subtractBackground(x,y,deg,indicesToIgnore=None):
#    print('subtractBackground: x = ',x)
#    print('subtractBackground: y = ',y)
    xArr = np.array(x)
    yArr = np.array(y)
    xFit = []
    yToFit = []
    for iX in np.arange(0,xArr.shape[0],1):
        if (not indicesToIgnore) or (iX not in indicesToIgnore):
            xFit.append(x[iX])
            yToFit.append(y[iX])
    xFit = np.array(xFit)
    yToFit = np.array(yToFit)
#    print('subtractBackground: xFit = ',xFit.shape,': ',xFit)
#    print('subtractBackground: yToFit = ',yToFit.shape,': ',yToFit)

    #interpolate lsq spline
#    xKnots = xFit[1::3]
#    print('subtractBackground: xKnots = ',xKnots.shape,': ',xKnots)

#    t = xKnots#[-1,-0.5,0,0.5,1]
#    k = 3
#    t = np.r_[(xFit[0],)*(k+1),
#              t,
#              (xFit[-1],)*(k+1)]
#    print('subtractBackground: t = ',t.shape,': ',t)
#    spl = make_lsq_spline(xFit, yToFit, t, k)
    nx = normalizeX(xFit)
    coeffs = np.polynomial.legendre.legfit(nx, yToFit, deg)
    nx = normalizeX(xArr)
    yFit = np.polynomial.legendre.legval(nx, coeffs)
#    plt.plot(x,y,label='y(x)')
#    plt.plot(xFit,yToFit,label='yToFit(xFit)')
#    plt.plot(x,yFit,label='yFit(x)')
#    plt.legend()
#    plt.show()

    return yArr-yFit

def calcDispersion(lineList, xRange, degree=3, delimiter=' ', display=False):
    pixels = [xRange[0]]# we need to include the range limits to get the correct normlization
    wLens = [-1]        # we need to include the range limits to get the correct normlization

    if isinstance(lineList,str):
        with open(lineList,'r') as f:
            lines = f.readlines()
        print('calcDispersion: lines = ',lines)

        for line in lines:
            line = line.rstrip()
            line = line.split(delimiter)
            pixels.append(float(line[0]))
            wLens.append(float(line[1]))
    else:
        for line in lineList:
            pixels.append(line[0])
            wLens.append(line[1])
    pixels.append(xRange[1])# we need to include the range limits to get the correct normlization
    wLens.append(-1)        # we need to include the range limits to get the correct normlization

    pixels = np.array(pixels)
    wLens = np.array(wLens)

    #normalize x to [-1,1]
    nx = np.array(normalizeX(pixels))
    print('calcDispersion: pixels = ',pixels.shape,': ',pixels)
    print('calcDispersion: nx = ',nx.shape,': ',nx)

    #fit Legendre polynomial excluding the normalized range limits
    coeffs = np.polynomial.legendre.legfit(nx[1:nx.shape[0]-1], wLens[1:nx.shape[0]-1], degree)

    # calculate fitted values excluding the normalized range limits
    yFit = np.polynomial.legendre.legval(nx[1:nx.shape[0]-1], coeffs)
    differences = wLens[1:nx.shape[0]-1] - yFit
    errors = (differences) ** 2
    rms = np.sqrt(np.sum(errors) / (wLens.shape[0]-2))
    for i in range(len(pixels)-2):
        print('calcDispersion: ',pixels[i+1],wLens[i+1],differences[i],errors[i])
    print('calcDispersion: RMS = ',rms)

    # plot original values and fit
    if display:
        plt.scatter(pixels[1:nx.shape[0]-1],wLens[1:nx.shape[0]-1])
        plt.errorbar(pixels[1:nx.shape[0]-1],wLens[1:nx.shape[0]-1],yerr=errors)
        plt.plot(pixels[1:nx.shape[0]-1],yFit)
        plt.show()

    #p = L.fit(pixels, wLens, 3)
    #print('calcDispersion: p = ',p)
    return [coeffs,rms]

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
    print('sigmaReject: iter = 0: indices = ',indices.shape,': ',indices)
    for iter in np.arange(1,nIter,1):
        yGood, goodIndices = sigmaRej(yGood, lowReject, highReject, replace=replace, adjustSigLevels=adjustSigLevels, useMean=useMean, keepFirstAndLastX=keepFirstAndLastX)
        indices = indices[goodIndices]
        print('sigmaReject: iter = ',iter,': yGood = ',yGood.shape,': ',yGood)
        print('sigmaReject: iter = ',iter,': indices = ',indices.shape,': ',indices)
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

# NOTE: only takes the median of the lower half of the sorted values in order to allow the sky subtraction for extended objects
def subtractMedianSky(twoDImageFilesIn):
    for twoDImageFileIn in twoDImageFilesIn:
        image = np.array(CCDData.read(twoDImageFileIn, unit="adu"))
        objectSpec = []
        skyImage = np.zeros(image.shape)
        for col in range(image.shape[1]):
            skyImage[:,col] = np.ones(image.shape[0]) * np.median(np.sort(image[:,col])[0:int(image.shape[0]/2)])
        image = image - np.array(skyImage)
        writeFits(skyImage, twoDImageFileIn, twoDImageFileIn[:-5]+'MedianSky.fits', metaKeys=None, metaData=None, overwrite=True)
        writeFits(image, twoDImageFileIn, twoDImageFileIn[:-5]+'-MedianSky.fits', metaKeys=None, metaData=None, overwrite=True)

#@brief: extract sky-subtracted object spectrum
#@param twoDImageFileIn: string: name of 2d input fits file
#@param specOut: string: name of 1d output fits file
#@param yRange: [int start,int stop]: extracted area is where y in [start, stop] (including stop)
#@param skyAbove=None: [int start,int stop]: sky above area is where y in [start, stop] (including stop)
#@param skyBelow=None: [int start,int stop]: sky below area is where y in [start, stop] (including stop)
#@param extractionMethod='sum': 'sum', 'median', 'skyMedian', or filename:
#       if 'sum': sum up pixels where y in [start, stop] (including stop)
#       if 'median': extracted value = (stop-start+1) * median(y in [start, stop] (including stop)) - for all sky images
#       if 'skyMedian': extracted value = sum up pixels where y in [start, stop] (including stop) - ((stop-start+1) * median(y in [start, stop]) (including stop))
#                       skyAbove and skyBelow should be empty, should be default method if no PN spectrum is obvious
#       if filename: scale sky image with exposure time and subtract from twoDImageFileIn
#@param dispAxis='row': 'row' or 'column'
def extractObjectAndSubtractSky(twoDImageFileIn,
                                specOut,
                                yRange,
                                skyAbove=None,
                                skyBelow=None,
                                extractionMethod='sum',
                                dispAxis='row',
                                display=False):
    image = np.array(CCDData.read(twoDImageFileIn, unit="adu"))
    print('extractObjectAndSubtractSky: image.shape = ',image.shape)
    print('extractObjectAndSubtractSky: twoDImageFileIn = <'+twoDImageFileIn+'>')
    print('extractObjectAndSubtractSky: yRange = ',yRange)
    print('extractObjectAndSubtractSky: skyAbove = ',skyAbove)
    print('extractObjectAndSubtractSky: skyBelow = ',skyBelow)

    hdulist = pyfits.open(twoDImageFileIn)
    head = hdulist[0].header

    plt.rcParams["figure.figsize"] = [15., 7.0]
    plt.rcParams["figure.autolayout"] = True

    newImageData = None
    skyData = None
    if dispAxis == 'row':
#        if (skyAbove is not None) and (skyBelow is not None):
        imageTransposed = image.transpose()
        print('extractObjectAndSubtractSky: imageTransposed.shape = ',image.shape)
        if display:
            plt.imshow(image,vmin=0.,vmax=1.5*np.mean(image))
            print('image.shape = ',image.shape)
            if skyAbove is not None:
                plt.plot([0.,image.shape[1]],[skyAbove[0],skyAbove[0]],'r-')
                plt.plot([0.,image.shape[1]],[skyAbove[1],skyAbove[1]],'r-')
            if skyBelow is not None:
                plt.plot([0.,image.shape[1]],[skyBelow[0],skyBelow[0]],'b-')
                plt.plot([0.,image.shape[1]],[skyBelow[1],skyBelow[1]],'b-')
            plt.plot([0.,image.shape[1]],[yRange[0],yRange[0]],'y-')
            plt.plot([0.,image.shape[1]],[yRange[1],yRange[1]],'y-')
            mng = plt.get_current_fig_manager()
            print('dir(mng) = ',dir(mng))
            mng.full_screen_toggle()
            plt.title('original image '+twoDImageFileIn)
            plt.show()
        newImageData,skyData = subtractSky(imageTransposed,skyAbove,skyBelow,sigLow=3.0,sigHigh=3.0) if ((skyAbove is not None) and (skyBelow is not None)) else [imageTransposed,None]
        if display:
            plt.imshow(newImageData.transpose(),vmin=0.,vmax=1.5*np.mean(newImageData))
            if skyAbove is not None:
                plt.plot([0.,image.shape[1]],[skyAbove[0],skyAbove[0]],'r-')
                plt.plot([0.,image.shape[1]],[skyAbove[1],skyAbove[1]],'r-')
            if skyBelow is not None:
                plt.plot([0.,image.shape[1]],[skyBelow[0],skyBelow[0]],'b-')
                plt.plot([0.,image.shape[1]],[skyBelow[1],skyBelow[1]],'b-')
            plt.plot([0.,image.shape[1]],[yRange[0],yRange[0]],'y-')
            plt.plot([0.,image.shape[1]],[yRange[1],yRange[1]],'y-')
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            plt.title('sky subtracted image '+twoDImageFileIn)
            plt.show()
            if skyData is not None:
                plt.imshow(skyData.transpose())
            if skyAbove is not None:
                plt.plot([0.,image.shape[1]],[skyAbove[0],skyAbove[0]],'r-')
                plt.plot([0.,image.shape[1]],[skyAbove[1],skyAbove[1]],'r-')
            if skyBelow is not None:
                plt.plot([0.,image.shape[1]],[skyBelow[0],skyBelow[0]],'b-')
                plt.plot([0.,image.shape[1]],[skyBelow[1],skyBelow[1]],'b-')
            plt.plot([0.,image.shape[1]],[yRange[0],yRange[0]],'y-')
            plt.plot([0.,image.shape[1]],[yRange[1],yRange[1]],'y-')
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            plt.title('sky image '+twoDImageFileIn)
            plt.show()
    else:
        if (skyAbove is not None) and (skyBelow is not None):
            if display:
                plt.imshow(image)
                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                plt.show()
            newImageData,skyData = subtractSky(image,skyAbove,skyBelow,sigLow=3.0,sigHigh=3.0)
            if display:
                plt.imshow(newImageData)
                plt.show()
                plt.imshow(skyData)
                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                plt.show()

    profile = extractSum(image,'column' if dispAxis == 'row' else 'row')
    if display:
        plt.plot(profile)
        plt.title(twoDImageFileIn+' profile')
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()

    # subtract sky
    if dispAxis == 'row':
        if display:
            if skyAbove is not None:
                plt.plot([skyAbove[0],skyAbove[0]],[np.amin(profile),np.amax(profile)],'b-')
                plt.plot([skyAbove[1],skyAbove[1]],[np.amin(profile),np.amax(profile)],'b-')
            if skyBelow is not None:
                plt.plot([skyBelow[0],skyBelow[0]],[np.amin(profile),np.amax(profile)],'b-')
                plt.plot([skyBelow[1],skyBelow[1]],[np.amin(profile),np.amax(profile)],'b-')
            plt.plot([yRange[0],yRange[0]],[np.amin(profile),np.amax(profile)],'g-')
            plt.plot([yRange[1],yRange[1]],[np.amin(profile),np.amax(profile)],'g-')
        rowsSky = []
        if skyAbove is not None:
            for i in np.arange(skyAbove[0],skyAbove[1]+1,1):
                rowsSky.append(i)
        if skyBelow is not None:
            for i in np.arange(skyBelow[0],skyBelow[1]+1,1):
                rowsSky.append(i)
        if (skyAbove is not None) or (skyBelow is not None):
            rowsSky = np.array(rowsSky)
            for col in range(image.shape[1]):
                skyArr = []
                for i in rowsSky:
                    skyArr.append(image[i,col])

                f = interp1d(rowsSky, np.array(skyArr), bounds_error = False,fill_value='extrapolate')
                image[yRange[0]:yRange[1]+1,col] -= f(np.arange(yRange[0],yRange[1]+1))
        if display:
            plt.show()
            plt.imshow(image,vmin=0.,vmax=1.5*np.mean(image))
            plt.title(twoDImageFileIn+' sky subtracted')
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            plt.show()

        if extractionMethod == 'sum':
#            if (skyAbove is not None) and (skyBelow is not None):
            objectSpec = extractSum(newImageData.transpose()[yRange[0]:yRange[1]+1,:],dispAxis)
#            else:
#                objectSpec = extractSum(image[yRange[0]:yRange[1]+1,:],dispAxis)
        elif extractionMethod == 'median':
            objectSpec = []
            for col in range(image.shape[1]):
                objectSpec.append(np.median(image[yRange[0]:yRange[1]+1,col]) * (yRange[1]-yRange[0]+1))
            objectSpec = np.array(objectSpec)
        elif extractionMethod == 'skyMedian':
            print('extractObjectAndSubtractSky: ERROR: "skyMedian" as method has been abandoned - use the ...MedianSky.fits image for the sky instead')
            STOP
            objectSpec = []
            skyImage = np.zeros(image.shape)
            for col in range(image.shape[1]):
                skyImage[yRange[0]:yRange[1]+1,col] = np.ones(yRange[1]-yRange[0]+1) * np.median(image[yRange[0]:yRange[1]+1,col])
                objectSpec.append(np.sum(image[yRange[0]:yRange[1]+1,col] - skyImage[yRange[0]:yRange[1]+1,col]))
            objectSpec = np.array(objectSpec)
        else:
            skyImage = np.array(CCDData.read(extractionMethod, unit="adu"))
            result = np.array([minimize(chisqfunc, 1., (image[:,iCol], skyImage[:,iCol], np.sqrt(np.absolute(image[:,iCol])))).x[0] for iCol in range(image.shape[1])])
            print('extractObjectAndSubtractSky: result = ',result)
            xNorm = normalizeX(np.arange(0,image.shape[1],1.))
            #sfit(x, y, fittingFunction, solveFunction, order, nIterReject, nIterFit, lowReject, highReject, adjustSigLevels=False, useMean=False, display=False)
            coeffs, resultFit = sfit(xNorm,
                                     result,
                                     np.polynomial.legendre.legfit,
                                     np.polynomial.legendre.legval,
                                     order=7,
                                     nIterReject=2,
                                     nIterFit=2,
                                     lowReject=3.,
                                     highReject=1.8,
                                     adjustSigLevels=False,
                                     useMean=False,)
            if display:
                plt.plot(result)
                plt.plot(resultFit)
                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                plt.show()

            for iCol in range(image.shape[1]):
                image[:,iCol] = image[:,iCol] - (skyImage[:,iCol])# * resultFit[iCol])
            objectSpec = extractSum(image[yRange[0]:yRange[1]+1,:],dispAxis)
    else:
        if display:
            if skyAbove is not None:
                plt.plot([np.amin(profile),np.amax(profile)],[skyAbove[0],skyAbove[0]],'b-')
                plt.plot([np.amin(profile),np.amax(profile)],[skyAbove[1],skyAbove[1]],'b-')
            if skyBelow is not None:
                plt.plot([np.amin(profile),np.amax(profile)],[skyBelow[0],skyBelow[0]],'b-')
                plt.plot([np.amin(profile),np.amax(profile)],[skyBelow[1],skyBelow[1]],'b-')
            plt.plot([np.amin(profile),np.amax(profile)],[yRange[0],yRange[0]],'g-')
            plt.plot([np.amin(profile),np.amax(profile)],[yRange[1],yRange[1]],'g-')
        colsSky = []
        if skyAbove is not None:
            for i in np.arange(skyAbove[0],skyAbove[1]+1,1):
                colsSky.append(i)
        if skyBelow is not None:
            for i in np.arange(skyBelow[0],skyBelow[1]+1,1):
                colsSky.append(i)
        colsSky = np.array(colsSky)
        if (skyAbove is not None) or (skyBelow is not None):
            for row in range(image.shape[0]):
                skyArr = []
                for i in colsSky:
                    skyArr.append(image[row,i])

                f = interp1d(colsSky, np.array(skyArr), bounds_error = False,fill_value='extrapolate')
                image[row,yRange[0]:yRange[1]+1] -= f(np.arange(yRange[0],yRange[1]+1))
        if display:
            plt.show()
        if extractionMethod == 'sum':
            if (skyAbove is not None) and (skyBelow is not None):
                objectSpec = extractSum(newImageData[:,yRange[0]:yRange[1]+1],dispAxis)
            else:
                objectSpec = extractSum(image[:,yRange[0]:yRange[1]+1],dispAxis)
        elif extractionMethod == 'median':
            objectSpec = []
            for row in range(image.shape[0]):
                objectSpec.append(np.median(row,image[yRange[0]:yRange[1]+1]) * (yRange[1]-yRange[0]+1))
            objectSpec = np.array(objectSpec)
        elif extractionMethod == 'skyMedian':
            objectSpec = []
            for row in range(image.shape[0]):
                objectSpec.append(np.sum(row,image[yRange[0]:yRange[1]+1]) - (np.median(row,image[yRange[0]:yRange[1]+1]) * (yRange[1]-yRange[0]+1)))
            objectSpec = np.array(objectSpec)
        else:
            skyImage = np.array(CCDData.read(extractionMethod, unit="adu"))
            result = np.array([minimize(chisqfunc, 1., (image[iRow,:], skyImage[iRow,:], np.sqrt(np.absolute(image[iRow,:])))).x[0] for iRow in range(image.shape[0])])
            print('extractObjectAndSubtractSky: result = ',result)
            xNorm = normalizeX(np.arange(0,image.shape[0],1.))
            #sfit(x, y, fittingFunction, solveFunction, order, nIterReject, nIterFit, lowReject, highReject, adjustSigLevels=False, useMean=False, display=False)
            coeffs, resultFit = sfit(xNorm,
                                     result,
                                     np.polynomial.legendre.legfit,
                                     np.polynomial.legendre.legval,
                                     order=7,
                                     nIterReject=2,
                                     nIterFit=2,
                                     lowReject=3.,
                                     highReject=1.8,
                                     adjustSigLevels=True,
                                     useMean=False,
                                     )
            if display:
                plt.plot(result)
                plt.plot(resultFit)
                mng = plt.get_current_fig_manager()
                mng.full_screen_toggle()
                plt.show()

            for iRow in range(image.shape[0]):
                image[iRow,:] = image[iRow,:] - (skyImage[iRow,:] * resultFit[iRow])
            objectSpec = extractSum(image[:,yRange[0]:yRange[1]+1],dispAxis)
#    if display:
#        mng = plt.get_current_fig_manager()
#        mng.full_screen_toggle()
#        plt.show()

    if newImageData is not None:
        if dispAxis == 'row':
            writeFits(newImageData.transpose(), twoDImageFileIn, twoDImageFileIn[:-5]+'-sky.fits', metaKeys=None, metaData=None, overwrite=True)
            if skyData is not None:
                writeFits(skyData.transpose(), twoDImageFileIn, twoDImageFileIn[:-5]+'Sky.fits', metaKeys=None, metaData=None, overwrite=True)
        else:
            writeFits(newImageData, twoDImageFileIn, twoDImageFileIn[:-5]+'-sky.fits', metaKeys=None, metaData=None, overwrite=True)
            if skyData is not None:
                writeFits(skyData, twoDImageFileIn, twoDImageFileIn[:-5]+'Sky.fits', metaKeys=None, metaData=None, overwrite=True)
    else:
        if extractionMethod == 'skyMedian':
            writeFits(image-skyImage, twoDImageFileIn, twoDImageFileIn[:-5]+'-sky.fits', metaKeys=None, metaData=None, overwrite=True)
        else:
            writeFits(image, twoDImageFileIn, twoDImageFileIn[:-5]+'-sky.fits', metaKeys=None, metaData=None, overwrite=True)

    if display:
        plt.plot(objectSpec,'g-')
        plt.title(twoDImageFileIn[twoDImageFileIn.rfind('SCIENCE')+8:twoDImageFileIn.rfind('.')])
        plt.show()

    writeFits1D(objectSpec, specOut, wavelength=None, header=head, CRVAL1=1., CRPIX1=1, CDELT1=1.)

#@param xSpec: np.array1d
#@param ySpec: np.array1d
#@param lineProfile: [x:np.array1d,y:np.array1d]
def findGoodLines(xSpec,ySpec,lineProfile,outFileNameAllLines=None,outFileNameGoodLines=None,display=False,chiSquareLimit=0.25):
    xXCor, xCorChiSquares = xCor([xSpec,ySpec],lineProfile)
    linesX = findLines(ySpec,
                       xXCor,
                       xCorChiSquares,
                       3.,
                       peakHeight=np.amax(ySpec) * 0.0025,#/ (300000. / 14000.),
                       peakWidth=2.5,#3.,
                       threshold=300.,
                       plot=display,
                      )
    print('findGoodLines: linesX = ',linesX)
    print('lineProfile = ',lineProfile)
#    print('lineProfile[0][0] = ',lineProfile[0][0])
#    print('lineProfile[0][len(lineProfile[0])-1] = ',lineProfile[0][len(lineProfile[0])-1])
#    goodProfIndices = np.where(np.absolute( lineProfile[0] ) < np.min([np.absolute(lineProfile[0][0]), lineProfile[0][len(lineProfile[0])-1]]))
#    print('goodProfIndices = ',goodProfIndices)
#    lineProfile[0] = lineProfile[0][goodProfIndices]
#    lineProfile[1] = lineProfile[1][goodProfIndices]
#    print('lineProfile = ',lineProfile)

    if display:
        plt.plot(ySpec)
        plt.title('findGoodLines')
    chiSquares = []
#    chiSquareLimit = 0.25#0.16
    for line in linesX:
        print('findGoodLines: line = ',line)
        print('findGoodLines: lineProfile[0]+line = ',lineProfile[0]+line)
        print('findGoodLines: np.arange(int(np.amin(lineProfile[0]+line)),int(np.amax(lineProfile[0]+line)),1) = ',np.arange(int(np.amin(lineProfile[0]+line)),int(np.amax(lineProfile[0]+line)),1))
        print('findGoodLines: ySpec[np.arange(int(np.amin(lineProfile[0]+line)),int(np.amax(lineProfile[0]+line)),1)] = ',ySpec[np.arange(int(np.amin(lineProfile[0]+line)),int(np.amax(lineProfile[0]+line)),1)])
        print('findGoodLines: np.amax(lineProfile[1]) = ',np.amax(lineProfile[1]))
        xs = np.arange(int(np.amin(lineProfile[0]+line)),
                       int(np.amax(lineProfile[0]+line)),
                       1)
        print('findGoodLines: line = ',line,': xs = ',xs)
        background = np.amin(ySpec[xs])
        top = np.amax(ySpec[xs])
        if display:
            plt.plot(lineProfile[0]+line,
                     lineProfile[1]
                     * (top-background)
                     / np.amax(lineProfile[1])
                     + background,
                     label=str(line)+'_1')
        profYatX = getYAt(lineProfile[0]+line,lineProfile[1],xs)
        profYatX = profYatX / np.amax(profYatX)
        lineY = ySpec[xs]
        lineY = lineY - np.amin(lineY)
        lineY = lineY / np.amax(lineY)
        chiSquare = np.sum((profYatX - lineY) ** 2)
        print('findGoodLines: line = ',line,': chiSquare = ',chiSquare)
        #chiSquares.append(chiSquare)#np.sum((profYatX[int(len(lineY) / 4.):int(len(lineY) * 3. / 4.)] - lineY[int(len(lineY) / 4.):int(len(lineY) * 3. / 4.)]) ** 2))
#        print('findGoodLines: line = ',line,': chiSquare = ',chiSquares[len(chiSquares)-1])


#        xs = np.arange(int(np.amin(lineProfile[0][int(len(lineProfile[0])/4.):]+line)),
#                       int(np.amax(lineProfile[0][:int(len(lineProfile[0])*3./4.)]+line)),
#                       1)
#        print('findGoodLines: line = ',line,': xs = ',xs)
#        background = np.amin(ySpec[xs])
#        top = np.amax(ySpec[xs])
#        if display:
#            plt.plot(lineProfile[0]+line,
#                     lineProfile[1]
#                     * (top-background)
#                     / np.amax(lineProfile[1])
#                     + background,
#                     label=str(line)+'_2')
#        profYatX = getYAt(lineProfile[0]+line,lineProfile[1],xs)
#        profYatX = profYatX / np.amax(profYatX)
#        lineY = ySpec[xs]
#        lineY = lineY - np.amin(lineY)
#        lineY = lineY / np.amax(lineY)
#        chiSquare = np.sum((profYatX - lineY) ** 2)
#        print('findGoodLines: line = ',line,': chiSquare = ',chiSquare)
        chiSquares.append(chiSquare)#np.sum((profYatX[int(len(lineY) / 4.):int(len(lineY) * 3. / 4.)] - lineY[int(len(lineY) / 4.):int(len(lineY) * 3. / 4.)]) ** 2))
#        print('findGoodLines: line = ',line,': chiSquare = ',chiSquares[len(chiSquares)-1])
        if display and (chiSquare < chiSquareLimit):
            plt.scatter(line,0.)

    print('findGoodLines: chiSquareLimit = ',chiSquareLimit)
    if plot:
        plt.legend()
        plt.show()
#    if True:
        plt.scatter(linesX,chiSquares)
        plt.plot([0,ySpec.shape[0]],[chiSquareLimit,chiSquareLimit])
        plt.title('findGoodLines chiSquares')
        plt.show()
    if not outFileNameAllLines is None:
        with open(outFileNameAllLines,'w') as f:
            for i in range(len(linesX)):
                f.write('%.5f \n' % (linesX[i]))
    goodLines = []
    for i in range(len(linesX)):
        if chiSquares[i] < chiSquareLimit:
            goodLines.append(linesX[i])
            print('findGoodLines: line = ',linesX[i],' is okay')
        else:
            print('findGoodLines: line = ',linesX[i],' with chiSquare=',chiSquares[i],' rejected')
    if not outFileNameGoodLines is None:
        with open(outFileNameGoodLines,'w') as f:
            for line in goodLines:
                f.write('%.5f \n' % (line))
    return goodLines

def getNumberOfApertures(databaseFileNameIn):
    nAps = 0
    with open(databaseFileNameIn,'r') as f:
        lines = f.readlines()
    for line in lines:
        elems = line.strip('\n').strip(' ').strip('\t').split('\t')
        print('getNumberOfApertures: elems = ',elems)
        if elems[0].strip(' ').strip('\t') == 'aperture':
            nAps = int(elems[1])
    print('getNumberOfApertures: nAps = ',nAps)
    return nAps

def getApWidth(databaseFileNameIn):
    with open(databaseFileNameIn,'r') as f:
        lines = f.readlines()
    low = 0.
    high = 0.
    for line in lines:
        elems = line.strip('\n').strip(' ').strip('\t').split('\t')
        print('getApWidth: elems = ',elems)
        if elems[0].strip(' ').strip('\t') == 'low':
            low = float(elems[1].split(' ')[0])
        if elems[0].strip(' ').strip('\t') == 'high':
            high = float(elems[1].split(' ')[0])
    width = high-low
    print('getApWidth: width = ',width)
    return width

def getLineProfiles(arcFitsName2D,
                    dxFit=0.1,
                    display=False,
                    apOffsetX = 0.):
    lineProfiles = []
    print('getLineProfiles: arcFitsName2D = <'+arcFitsName2D+'>')
    tempFile = os.path.join(arcFitsName2D[0:arcFitsName2D.rfind('/')],'database','aptmp'+arcFitsName2D[arcFitsName2D.rfind('/')+1:-5])
    print('getLineProfiles: tempFile = <'+tempFile+'>')

    halfWidth = int(getApWidth(tempFile)/2.)

    for apNumber in np.arange(0,getNumberOfApertures(tempFile),1):
        lineProfiles.append(calcLineProfile(arcFitsName2D,
                                            apNumber,
                                            halfWidth,
                                            dxFit,
                                            apOffsetX=apOffsetX))

        if display:
            plt.plot(lineProfiles[len(lineProfiles)-1][0],lineProfiles[len(lineProfiles)-1][1],label='ap '+str(apNumber))
    if display:
        plt.title('getLineProfiles')
        plt.legend()
        plt.show()
    return lineProfiles

def getBestLineProfile(lineProfiles,outFileName=None,display=False):
    bestLineProfileIdx = 0
    maxValue = np.amax(lineProfiles[0][1])
    for lineProfileIdx in np.arange(1,len(lineProfiles),1):
        if np.amax(lineProfiles[lineProfileIdx][1]) > maxValue:
            bestLineProfileIdx = lineProfileIdx
    print('getBestLineProfile: best line profile found at index ',bestLineProfileIdx)
#    print('getBestLineProfile: lineProfiles = ',lineProfiles)
    print('getBestLineProfile: lineProfiles[bestLineProfileIdx] = ',lineProfiles[bestLineProfileIdx])
    if display:
        plt.plot(lineProfiles[bestLineProfileIdx][0],lineProfiles[bestLineProfileIdx][1])
        plt.show()
    if not outFileName is None:
        with open(outFileName,'w') as f:
            for i in np.arange(0,len(lineProfiles[bestLineProfileIdx][0]),1):
                f.write('%.5f %.5f\n' % (lineProfiles[bestLineProfileIdx][0][i],lineProfiles[bestLineProfileIdx][1][i]))
    return lineProfiles[bestLineProfileIdx]

#@brief: calculate distances for all lines in lineList to all lines in referenceLineList and identify the lines in lineList
#@param lineList: array like (1D)
#@param referenceLineList: array like (2D) [[pixel, wavelength],[...]]
#@return: array2D [[pixel, wavelength],[...]]
def crossCheckLines(lineList, referenceLineList):
    print('crossCheckLines: referenceLineList = ',len(referenceLineList),': ',referenceLineList)
    minDists = []
    minIdx = []
    for line in lineList:
        dists = []
        idx = []
        for iRef in range(len(referenceLineList)):
            dists.append(np.absolute(line-referenceLineList[iRef][0]))
            idx.append(iRef)
#        print('crossCheckLines: dists = ',dists)
#        print('crossCheckLines: idx = ',idx)
        dists = np.array(dists)
        idx = np.array(idx)
        minDists.append(np.amin(dists))
        where = np.where(dists == minDists[len(minDists)-1])
#        print('crossCheckLines: where = ',where)
        minIdx.append(idx[where[0][0]])
    print('crossCheckLines: minDists = ',minDists)
    print('crossCheckLines: minIdx = ',minIdx)

    meanMinDist = np.mean(minDists)/4.
    print('crossCheckLines: meanMinDist = ',meanMinDist)
    goodLines = []
    for i in range(len(minDists)):
        if minDists[i] < meanMinDist:
            print('crossCheckLines: referenceLineList[',minIdx[i],'] = ',referenceLineList[minIdx[i]],', referenceLineList[',minIdx[i],'][1] = ',referenceLineList[minIdx[i]][1])
            goodLines.append([lineList[i],referenceLineList[minIdx[i]][1]])
    if len(goodLines) < len(referenceLineList) / 4.:
        meanMinDist = np.median(minDists) * 1.3
        print('crossCheckLines: meanMinDist = ',meanMinDist)
        goodLines = []
        for i in range(len(minDists)):
            if minDists[i] < meanMinDist:
                print('crossCheckLines: referenceLineList[',minIdx[i],'] = ',referenceLineList[minIdx[i]],', referenceLineList[',minIdx[i],'][1] = ',referenceLineList[minIdx[i]][1])
                goodLines.append([lineList[i],referenceLineList[minIdx[i]][1]])
    return goodLines

# @brief: reidentify ARC lines and return new lineList
# @param specIn: array like (1D)
# @param lineListIn: array like (2D) [[position, wavelength], [position,wavelength],...]
# @param profileIn: array like (2D): [x,y]: x,y: array like (1D): integral normalized emission line profile
# @return: lineListOut: same as lineListIn with new positions
def reidentify(arcFitsName2D,
               arcFitsName2DForLineProfile,
               referenceApertureDefinitionFile,
               lineListIn,
               lineListOut=None,
               specOut=None,
               display=False,
               chiSquareLimit=0.25,
               degree=5,
               apOffsetX=0.):
    print('lineListIn = <'+lineListIn+'>')
#    STOP
    with open(referenceApertureDefinitionFile,'r') as f:
        lines = f.readlines()
    with open(referenceApertureDefinitionFile,'w') as f:
        for line in lines:
            f.write(line.replace(referenceApertureDefinitionFile[referenceApertureDefinitionFile.rfind('/')+3:],
                                 arcFitsName2DForLineProfile[arcFitsName2DForLineProfile.rfind('/')+1:arcFitsName2DForLineProfile.rfind('.')]))
    print('reidentify: updated new database file ',referenceApertureDefinitionFile)
#        STOP

    lineProfiles = getLineProfiles(arcFitsName2DForLineProfile, display=display, apOffsetX=apOffsetX)
    bestLineProfile = getBestLineProfile(lineProfiles,outFileName=None,display=display)
    specY = extractSum(arcFitsName2D,'row')
    writeFits1D(specY, arcFitsName2D[:arcFitsName2D.rfind('.')]+'Ec.fits', wavelength=None, header=arcFitsName2D, CRVAL1=1., CRPIX1=1., CDELT1=1.)
    if display:
        plt.plot(specY)
        plt.show()
    print('reidentify: specIn = ',specY.shape,': ',specY)
    print('reidentify: profileIn = ',bestLineProfile)

    goodLines = findGoodLines(np.arange(0,specY.shape[0],1),specY,bestLineProfile,outFileNameGoodLines=lineListOut,display=display,chiSquareLimit=chiSquareLimit)
    print('reidentify: found ',len(goodLines),' goodLines at ',goodLines)
    if not lineListOut is None:
        with open(lineListOut[:lineListOut.rfind('.')]+'_temp.dat','w') as f:
            for line in goodLines:
                f.write('%.5f\n' % line)

    try:
        gratingAngle = getHeaderValue(arc, 'GR-ANGLE', hduNum=0).strip()
    except:
        gratingAngle = '0.0'
    if '%' in lineListIn:
        lineList = readFileToArr(lineListIn % (int(gratingAngle[:gratingAngle.find('.')]),int(gratingAngle[gratingAngle.find('.')+1:])))
    else:
        lineList = readFileToArr(lineListIn)
    for line in lineList:
        print('reidentify: line = <'+line+'>')
        print('reidentify: line[:line.find(' ')] = <'+line[:line.find(' ')]+'>, line[line.find(' ')+1:] = <'+line[line.find(' ')+1:]+'>')
    refLineList = [[float(line[:line.find(' ')]),float(line[line.find(' ')+1:])] for line in lineList]
    print('refLineList = ',refLineList)
    lineListIdentified = crossCheckLines(goodLines,refLineList)
    print('reidentify: ',arcFitsName2D,': lineListIdentified = ',len(lineListIdentified),': ',lineListIdentified)
    if not lineListOut is None:
        with open(lineListOut,'w') as f:
            for line in lineListIdentified:
                f.write('%.5f %.5f\n' % (line[0],line[1]))

    xSpec = np.arange(0,specY.shape[0],1.)
    coeffs, rms = calcDispersion(lineListIdentified, xRange=[0,xSpec[xSpec.shape[0]-1]], degree=degree, display=display)
    if display:
        xSpecNorm = normalizeX(xSpec)
        wLenSpec = np.polynomial.legendre.legval(xSpecNorm, coeffs)
        plt.plot(wLenSpec,specY)
        minY = np.amin(specY)
        maxY = np.amax(specY)
        for line in lineListIdentified:
            plt.plot([line[1],line[1]],[minY,maxY])
        plt.title(arcFitsName2D[arcFitsName2D.rfind('/')+1:])
        plt.show()

    return [lineListIdentified, coeffs, [0,xSpec[xSpec.shape[0]-1]], rms]

def rebin(wavelength, spectrum, newWavelength, preserveFlux = True, outFileName = None, header = None):
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
        naNPos = naNPos[0]
    if naNPos.shape[0] > 0:
        specInterp.data[naNPos] = 0.

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

def shiftLineList(lineListIn,lineListOut,shift):
    with open(lineListIn,'r') as f:
        lines = f.readlines()
    print('lines = ',lines)
    lines = [line.strip() for line in lines]
    lines = [line.split(' ') for line in lines]
    print('lines = ',lines)
    lines = [[float(line[0])+shift,float(line[1])] for line in lines]
    with open(lineListOut,'w') as f:
        for line in lines:
            f.write('%.2f %.5f\n' % (line[0], line[1]))

def shiftApertureDefs(apDefFileIn,apDefFileOut,shift):
    print('apDefFileIn = ',apDefFileIn)
    with open(apDefFileIn,'r') as f:
        lines = f.readlines()
    with open(apDefFileOut,'w') as f:
        for line in lines:
            if line.find('center') >= 0:
                oldCenterTemp = line[line.find('center')+7:]
                print('shiftApertureDefs: oldCenterTemp = <'+oldCenterTemp+'>')
                oldCenter = float(oldCenterTemp[:oldCenterTemp.find(' ')])
                print('shiftApertureDefs: oldCenter = <',oldCenter,'>')
                newLine = line[:line.find('center')+7]+'%.5f' % (oldCenter+shift)
                f.write(newLine+oldCenterTemp[oldCenterTemp.find(' '):])
            else:
                f.write(line)

def getListOfFiles(fname):
    fList = readFileToArr(fname)
    print('fname = ',fname,': fList = ',fList)
    if len(fList) > 0:
        if fList[0].rfind('/') == -1:
            fList = [os.path.join(workPath, fileName) for fileName in fList]
    return fList

def resampleSpec(wLen, spec):
    dLam = np.min([np.absolute(wLen[1]-wLen[0]),np.absolute(wLen[wLen.shape[0]-1]-wLen[wLen.shape[0]-2])])
    resampled = np.arange(np.min([wLen[0], wLen[wLen.shape[0]-1]]), np.max([wLen[0], wLen[wLen.shape[0]-1]]), dLam)

    print('resampleSpec: wLen = ',wLen.shape,': ',wLen)
    print('resampleSpec: spec = ',spec.shape,': ',spec)
    print('resampleSpec: resampled = ',resampled.shape,': ',resampled)
    resampledSpec = rebin(wLen, spec, resampled, preserveFlux = False)
    return [resampled,resampledSpec]

def extractAndReidentifyARCs(arcListIn, refApDef, lineListIn, xCorSpecIn, display=False, chiSquareLimit=0.25,degree=5, apOffsetX=0.):
    print('refApDef = <'+refApDef)
    wavelengthsOrigOut = []
    wavelengthsResampledOut = []
    xCorSpec = getImageData(xCorSpecIn,0)#[12:-12]
    xCorSpecX = np.arange(0,len(xCorSpec),1)#[12:-12] - (len(xCorSpec) / 2.)
    for arc in arcListIn:
        try:
            gratingAngle = getHeaderValue(arc, 'GR-ANGLE', hduNum=0).strip()
        except:
            gratingAngle = '0.0'

        arcInterp = arc[:-5]+'i.fits'
        oneDSpecInterp = extractSum(arcInterp,'row')

        writeFits1D(oneDSpecInterp,
                    arcInterp[:-5]+'Ec.fits',
                    wavelength=None,
                    header=arcInterp,
                    CRVAL1=1,
                    CRPIX1=1,
                    CDELT1=1)
        oneDSpecInterpX = np.arange(0,len(oneDSpecInterp),1)
        if display:
            plt.plot(oneDSpecInterp)
            plt.show()
        print('extractAndReidentifyARCs: len(oneDSpecInterp) = ',len(oneDSpecInterp),', len(xCorSpec) = ',len(xCorSpec))
        corr = np.correlate(oneDSpecInterp,xCorSpec,'same')
        print('extractAndReidentifyARCs: corr = ',len(corr),': ',corr)
        if display:
            plt.plot(oneDSpecInterpX,oneDSpecInterp,label='spectrum')
            plt.plot(xCorSpecX,xCorSpec,label='reference spectrum')
            plt.legend()
            plt.show()
            plt.plot(corr)
            plt.show()
        maxPos = np.where(corr == np.max(corr))[0]
        print('extractAndReidentifyARCs: maxPos = ',maxPos)
        shift = maxPos - (len(xCorSpec) / 2.)
        print('extractAndReidentifyARCs: shift = ',shift)
        if display:
            plt.plot(oneDSpecInterpX,oneDSpecInterp,label='spectrum')
            plt.plot(xCorSpecX+shift,xCorSpec,label='reference spectrum')
            plt.legend()
            plt.show()
        if '%' in lineListIn:
            lineListTmp = lineListIn % (int(gratingAngle[:gratingAngle.find('.')]),int(gratingAngle[gratingAngle.find('.')+1:]))
        else:
            lineListTmp = lineListIn
        shiftLineList(lineListTmp,
                      lineListTmp+'tmp',
                      shift)
        print('extractAndReidentifyARCs: shifted line list = <'+lineListTmp+'tmp>')
        print('extractAndReidentifyARCs: referenceApertureDefinitionFile = <'+refApDef+'>')
        tempFile = arc[:arc.rfind('/')+1]+'database/aptmp'+arc[arc.rfind('/')+1:arc.rfind('.')]
        print('extractAndReidentifyARCs: tempFile = <'+tempFile+'>')
        if '%' in refApDef:
            refApDefTmp = refApDef % (int(gratingAngle[:gratingAngle.find('.')]),
                             int(gratingAngle[gratingAngle.find('.')+1:]))
        else:
            refApDefTmp = refApDef
        copyfile(refApDefTmp,tempFile)
        shiftApertureDefs(tempFile,tempFile,shift)
        lineListNew, coeffs, xRange, rms = reidentify(arcInterp,
                                                      arc,
                                                      tempFile,
                                                      lineListIn+'tmp',
                                                      lineListOut=arc[:arc.rfind('.')]+'_lines.dat',
                                                      specOut=arc[:-5]+'Ecd.fits',
                                                      display=display,
                                                      chiSquareLimit=chiSquareLimit,
                                                      degree=degree,
                                                      apOffsetX=apOffsetX)
        xSpec = np.arange(xRange[0],xRange[1]+1,1.)
        xSpecNorm = normalizeX(xSpec)
        wLenSpec = np.polynomial.legendre.legval(xSpecNorm, coeffs)
        wavelengthsOrigOut.append(wLenSpec)

        print('extractAndReidentifyARCs: np.absolute(wLenSpec[1]-wLenSpec[0]) = ',np.absolute(wLenSpec[1]-wLenSpec[0]))
        print('extractAndReidentifyARCs: np.absolute(wLenSpec[wLenSpec.shape[0]-1]-wLenSpec[wLenSpec.shape[0]-2]) = ',np.absolute(wLenSpec[wLenSpec.shape[0]-1]-wLenSpec[wLenSpec.shape[0]-2]))
        resampled,resampledSpec = resampleSpec(wLenSpec,oneDSpecInterp)
        wavelengthsResampledOut.append(resampled)

        if display:
            plt.plot(wLenSpec,oneDSpecInterp,label='original')
            plt.plot(resampled,resampledSpec,label='resampled')
            plt.legend()
            plt.title(arc)
            plt.show()

        writeFits1D(resampledSpec,
                    arcInterp[:-5]+'Ecd.fits',
                    wavelength=None,
                    header=arcInterp,
                    CRVAL1=resampled[0],
                    CRPIX1=1,
                    CDELT1=resampled[1]-resampled[0])

        pyfits.setval(arcInterp[:-5]+'Ecd.fits', 'CRVAL1', value=resampled[0])
        pyfits.setval(arcInterp[:-5]+'Ecd.fits', 'CRPIX1', value=1)
        pyfits.setval(arcInterp[:-5]+'Ecd.fits', 'CDELT1', value=resampled[1]-resampled[0])
        pyfits.setval(arcInterp[:-5]+'Ecd.fits', 'NLINES', value=len(lineListNew))
        for iLine in range(len(lineListNew)):
            pyfits.setval(arcInterp[:-5]+'Ecd.fits', 'LINE'+str(iLine), value='%.5f %.4f' % (lineListNew[iLine][0],lineListNew[iLine][1]))
        pyfits.setval(arcInterp[:-5]+'Ecd.fits', 'RMS', value=rms)

    return [wavelengthsOrigOut, wavelengthsResampledOut]

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
    closestBefore = getClosestInTime(specTime, arcTimes[whereLT])
    closestTemp = getClosestInTime(specTime, arcTimes[whereGT])
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
        wLenSpecAfter = wavelengthsOrigIn[closestArcs[1][0]]
        print('dispcor: wLenSpecAfter = ',wLenSpecAfter)
        factorAfter = closestArcs[1][1] / (closestArcs[0][1]+closestArcs[1][1])
        print('dispcor: factorAfter = ',factorAfter)

        wLenSpec = (wLenSpecBefore * factorBefore) + (wLenSpecAfter * factorAfter)
        print('dispCor: wLenSpec = ',wLenSpec)

        #read science header and append keywords
        headerSc['REFSPEC1'] = arcBefore[arcBefore.rfind('/')+1:]+' %.5f' % (factorBefore)
        headerSc['REFSPEC2'] = arcAfter[arcAfter.rfind('/')+1:]+' %.5f' % (factorAfter)

        #apply heliocentric radial velocity correction
        if doHelioCor:
            vrad = heliocor(observatoryLocation, headerSc, keywordRA, keywordDEC, keywordObsTime)
            headerSc['VHELIO'] = vrad
            wLenSpec = applyVRadCorrection(wLenSpec, vrad)
            print('dispCor: after heliocentric correction for vrad = ',vrad,': wLenSpec = ',wLenSpec)
        hdulist.close()

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
#        if 'dbs01541' in scienceListIn[iSpec]:
#            STOP

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

def applyVRadCorrection(wavelength, vRad):
    return np.array([wLen - (vRad * wLen / c0) for wLen in wavelength])

def readFluxStandardsList(fName='/Users/azuri/stella/referenceFiles/fluxStandards.txt'):
    with open(fName,'r') as f:
        lines = f.readlines()
    lines = [line.rstrip('\n') for line in lines]
    names = [line[line.rfind('/')+1:line.rfind('.')] for line in lines]
    dirs = [line[:line.rfind('/')] for line in lines]
    dirs = [d[d.rfind('/')+1:] for d in dirs]
#    stdDirs = []
#    for d in dirs:
#        if not d in stdDirs:
#            stdDirs.append(d)
    return [names, dirs, lines]

def readFluxStandardFile(fName):
    #with open(fName,'r') as f:
        #lines = f.readlines()
    lines = readFileToArr(fName)[1:]#[line.rstrip('\n').strip('\t').strip(' ') for line in lines]
    #print('readFluxStandardFile: lines = ',lines)
    lines = [line.replace('\t',' ') for line in lines]
    wavelengths = np.asarray([float(line[:line.find(' ')]) for line in lines])
    lines = [line[line.find(' '):].strip() for line in lines]
    fluxes = np.asarray([float(line[:line.find(' ')]) for line in lines])
    for i in range(len(lines)):
        print('readFluxStandardFile: lambda = ',wavelengths[i],', flux = ',fluxes[i])
    return [wavelengths, fluxes]

def calcResponse(fNameList,
                 arcList,
                wLenOrig,
                areas,
                stdStarNameEndsBefore = '_a',
                fluxStdandardList = '/Users/azuri/stella/referenceFiles/fluxStandards.txt',
                airmassExtCor='apoextinct.dat',
                display=False):
    fluxStandardNames, fluxStandardDirs, fluxStandardFileNames = readFluxStandardsList(fluxStdandardList)
    print('calcResponse: fluxStandardNames = ',fluxStandardNames)
    print('calcResponse: fluxStandardDirs = ',fluxStandardDirs)
    print('calcResponse: fluxStandardFileNames = ',fluxStandardFileNames)
    fluxStandardNames = np.asarray(fluxStandardNames)
    fluxStandardDirs = np.asarray(fluxStandardDirs)
    fluxStandardFileNames = np.asarray(fluxStandardFileNames)
    print('calcResponse: fluxStandardNames = ',fluxStandardNames)
    print('calcResponse: fluxStandardDirs = ',fluxStandardDirs)
    print('calcResponse: reading fNameList <'+fNameList+'>')
    with open(fNameList, 'r') as f:
        fNames = f.readlines()
    fNames = [n.rstrip('\n') for n in fNames]
    sensFuncs = []

    print('calcResponse: fNames = ',fNames)
    for fName in fNames:
        print("calcResponse: fName.rfind('SCIENCE_') = ",fName.rfind('SCIENCE_'))
        if fName.rfind('SCIENCE_') >= 0:
            stdName = fName[fName.rfind('SCIENCE_')+8:]
        else:
            stdName = fName[fName.rfind('FLUXSTDS_')+9:]
        stdName = stdName[:stdName.find(stdStarNameEndsBefore)].replace('_','').replace('-','')
        print('calcResponse: stdName = <'+stdName+'>')
#        wLenStd = getWavelengthArr(fName[:fName.rfind('.')]+'Ecd.fits',0)
#        wLenStd = getWavelengthArr(fName,0)

        indices = np.where(fluxStandardNames == stdName.lower())[0]
        print('calcResponse: indices = ',indices)
        print('calcResponse: Flux standard star '+stdName+' found in ',fluxStandardDirs[indices])
        print('calcResponse: Flux standard star '+stdName+' found in ',fluxStandardFileNames[indices])
        done = False
        for prior in fluxStdDirsByPriority:
            print('calcResponse: prior = <'+prior+'>')
            priorInFluxStandardDirs = False
            if prior.find('/') < 0:
                print("calcResponse: prior.rfind('/') = ",prior.find('/')," < 0")
                if prior in fluxStandardDirs[indices]:
                    print('calcResponse: prior( = <'+prior+'> in fluxStandardDirs[',indices,'] = ',fluxStandardDirs[indices])
                    priorInFluxStandardDirs = True
                    print('calcResponse: priorInFluxStandardDirs = True')
                else:
                    print('calcResponse: prior( = <'+prior+'> NOT in fluxStandardDirs[',indices,'] = ',fluxStandardDirs[indices])
            else:
                print("calcResponse: prior.rfind('/') = ",prior.rfind('/')," >= 0")
                if prior[prior.rfind('/')+1:] in fluxStandardDirs[indices]:
                    print("calcResponse: prior[prior.rfind('/')+1:] = <"+prior[prior.rfind('/')+1:]+"> in fluxStandardDirs[",indices,"] = ",fluxStandardDirs[indices])
                    priorInFluxStandardDirs = True
                    print('calcResponse: priorInFluxStandardDirs = True')
                else:
                    print("calcResponse: prior[prior.rfind('/')+1:] = <"+prior[prior.rfind('/')+1:]+"> NOT in fluxStandardDirs[",indices,"] = ",fluxStandardDirs[indices])

            print('calcResponse: priorInFluxStandardDirs = ',priorInFluxStandardDirs)
            if (not done) and (priorInFluxStandardDirs):
                print('calcResponse: dir = '+prior)
                img = CCDData.read(fName, unit=u.adu)
                print('calcResponse: dir(img.header) = ',dir(img.header))
                for key in img.header.keys():
                    print(key,': ',img.header[key])
                print('calcResponse: img.data.shape = ',img.data.shape)
                print('calcResponse: type(img.data[0]) = ',type(img.data[0]))
                print("calcResponse: type(img.header['EXPTIME']) = ",type(img.header['EXPTIME']))
                # put in units of ADU/s
                #img.data = img.data / float(img.header['EXPTIME'])
                #img.unit = u.adu / u.s

                # Trace & Extract the standard star spectrum. See the extract example demo for more details
#                tr = trace(img, display=True, nbins=25)
#                print('calcResponse: tr = ',tr)
#                ex_tbl = extract(img, tr, display=True, apwidth=8, skysep=3, skywidth=7)
                #scienceSpectra = []#extractSum(fn,'row',fn[:-5]+'Ec.fits') for fn in getListOfFiles(os.path.join(workPath,'SCIENCE_otzfif.list'))]
                extractedFileName = ''
                print('fName = ',fName)
                foundEx = False
                for i in range(areas.size()):
                    print('')
                    if areas.getData('fName',i) == fName:
                        extractedFileName = areas.getData('fName',i)[:-5]+'Ecd.fits'
                        print('calcResponse: reading extractedFileName = <'+extractedFileName+'>')
                        wapprox = getWavelengthArr(extractedFileName)
                        foundEx = True
                if not foundEx:
                    print('calcResponse: did not find fName = <'+fName+'> in areas')
                    STOP
                obj_flux = CCDData.read(extractedFileName, unit=u.adu)
                obj_flux = obj_flux.data / float(img.header['EXPTIME'])
                obj_flux = obj_flux * u.adu / u.s

                # this data comes from the APO DIS red channel, which has wavelength axis backwards
                # (despite not mentioning in the header)
                # wapprox = wLenOrig[getClosestArcs(fName,arcList)]#(np.arange(img.shape[1]) - img.shape[1]/2)[::-1] * img.header['DISPDW'] + img.header['DISPWC']

                wapprox = wapprox * u.angstrom

                # obj_flux = (flux_std - sky_std) * u.adu / u.s
                #obj_flux = ex_tbl['flux'] - ex_tbl['skyflux']
                print('calcResponse: type(obj_flux) = ',type(obj_flux))
                print('calcResponse: dir(obj_flux) = ',dir(obj_flux))
                print('calcResponse: obj_flux = ',obj_flux)
                print('calcResponse: obj_flux.data = ',obj_flux.data)

                if display:
                    plt.plot(wapprox, obj_flux.data)
#                    plt.errorbar(wapprox.value, obj_flux.data, alpha=0.25)#, yerr=ex_tbl['fluxerr'].data
                    plt.show()

                print('calcResponse: prior = <'+prior+'>')
                print("calcResponse: prior+'/'+stdName.lower()+'.dat' = <"+prior+'/'+stdName.lower()+'.dat'+'>')
                stdstar=onedstd(prior+'/'+stdName.lower()+'.dat')
                print('calcResponse: stdstar = ',stdstar)

                #obj_flux = ex_tbl['flux'] - ex_tbl['skyflux']
#                print('calcResponse: obj_flux.quantity = ',obj_flux.quantity)
                obj_spectrum = Spectrum1D(spectral_axis=wapprox,
                                          flux=obj_flux)#.quantity,)
                                          #uncertainty=StdDevUncertainty(ex_tbl['fluxerr']))

                print('calcResponse: obj_spectrum = ',obj_spectrum,', stdstar = ',stdstar)
                sensfunc_lin = standard_sensfunc(obj_spectrum, stdstar, display=True, mode='linear')
                print('calcResponse: sensfunc_lin = ',sensfunc_lin)
                # the actual sensitivity function(s), which in theory include some crude information about
                # the flat fielding (response) - though the reference spectrum is very coarse.
                if display:
                    plt.plot(sensfunc_lin['wave'], sensfunc_lin['S'])
                    plt.show()

                # now apply the sensfunc back to the std star to demonstrate
                # NOTE: this only works b/c wavelength is exactly the same. Normally use `apply_sensfunc`
                if display:
                    plt.plot(wapprox, obj_flux * sensfunc_lin['S'])
                    plt.scatter(stdstar['wave'], stdstar['flux'], c='C1')
                    plt.xlim(5500,7500)
                    plt.ylim(0, 0.3e-12)
                    plt.show()

                # now let's demo the Airmass correction
                Xfile = obs_extinction(airmassExtCor)

                try:
                    AIRVAL = float(img.header['AIRMASS'])
                except:
                    AIRVAL = float(img.header['SECZ'])
                print('calcResponse: AIRMASS = ',AIRVAL)
                Atest = airmass_cor(obj_spectrum, AIRVAL, Xfile)

                if display:
                    plt.plot(obj_spectrum.wavelength, obj_spectrum.flux)
                    plt.plot(Atest.wavelength, Atest.flux)
                    plt.show()

                sensFuncs.append(sensfunc_lin)

                # Now demo how to apply a sensfuc to a new spectrum (just happens to be the same spectrum here...)
                #Stest = apply_sensfunc(obj_spectrum, sensfunc_lin)

                #print('calcResponse: Stest = ',Stest)

                #plt.scatter(stdstar['wave'], stdstar['flux'], c='C1')

                #plt.plot(Stest.wavelength, Stest.flux, c='C2')
                #plt.xlim(5500,7500)
                #plt.ylim(0, 0.3e-12)
                #plt.show()

                done = True

            #wavelengths, fluxes = readFluxStandardFile(fluxStandardFileNames[ind])
 #           if (wavelengths[0] < wLenStd[0]) and (wavelengths[len(wavelengths)-1] > wLenStd[len(wLenStd)-1]):
 #               goodFiles.append([fName,fluxStandardFileNames[ind]])
    return sensFuncs



def applySensFuncs(objectSpectraIn, objectSpectraOut, sensFuncs, airmassExtCor='apoextinct.dat'):
    for iSpec in range(len(objectSpectraIn)):
        print('applySensFuncs: reading spectrum file <'+objectSpectraIn[iSpec]+'>')
        img = CCDData.read(objectSpectraIn[iSpec], unit=u.adu)
        # put in units of ADU/s
        img.data = img.data / float(img.header['EXPTIME'])
        img.unit = u.adu / u.s
#        print('applySensFuncs: img.header = ',img.header)
#        print('applySensFuncs: dir(img.header) = ',dir(img.header))
#        print('applySensFuncs: img.header.keys = ',img.header.keys)

        wLen = getWavelengthArr(objectSpectraIn[iSpec],0) * u.angstrom

        obj_spectrum = Spectrum1D(spectral_axis=wLen, flux=img.data * img.unit)#,
 #                                 uncertainty=StdDevUncertainty(ex_tbl['fluxerr']))
        Xfile = obs_extinction(airmassExtCor)

        try:
            AIRVAL = float(getHeaderValue(objectSpectraIn[iSpec],'AIRMASS'))#img.header['AIRMASS'])
        except:
            AIRVAL = float(img.header['SECZ'])
#        print('applySensFuncs: AIRMASS = ',AIRVAL)
#        print('applySensFuncs: obj_spectrum = ',obj_spectrum)
        obj_spectrum = airmass_cor(obj_spectrum, AIRVAL, Xfile)
#        print('applySensFuncs: obj_spectrum = ',obj_spectrum)
#        print('applySensFuncs: sensFuncs[0] = ',sensFuncs[0])
        objectSpectrumFluxCalibrated = apply_sensfunc(obj_spectrum, sensFuncs[0])
#        print('applySensFuncs: objectSpectrumFluxCalibrated.data = ',objectSpectrumFluxCalibrated.data)
#        print('applySensFuncs: dir(objectSpectrumFluxCalibrated.data) = ',dir(objectSpectrumFluxCalibrated.data))
#        print('applySensFuncs: img.header.keys = ',img.header.keys)
        crval = obj_spectrum.wavelength[0]
        cdelt = obj_spectrum.wavelength[1]-obj_spectrum.wavelength[0]
        crpix = 1
#        print('applySensFuncs: crval = ',crval)
#        print('applySensFuncs: dir(crval) = ',dir(crval))
#        print('applySensFuncs: crval.value = ',crval.value)
#        print('applySensFuncs: dir(crval.value) = ',dir(crval.value))
#        print('applySensFuncs: cdelt = ',cdelt)
#        print('applySensFuncs: dir(cdelt) = ',dir(cdelt))
#        print('applySensFuncs: crpix = ',crpix)
#        print('applySensFuncs: objectSpectrumFluxCalibrated.data = ',objectSpectrumFluxCalibrated.data)
        print('applySensFuncs: writing file ',objectSpectraOut[iSpec])
        writeFits1D(objectSpectrumFluxCalibrated.data,
                    objectSpectraOut[iSpec],
                    wavelength=None,
                    header=img.header,
                    CRVAL1=crval.value,
                    CRPIX1=crpix,
                    CDELT1=cdelt.value,
                   )

if False:#def fluxCalibrate(obsSpecFName, standardSpecFName):
    spec = getImageData(obsSpecFName,0)
    wLen = getWavelengthArr(obsSpecFName,0)
    objectSpectrum = Spectrum1D( flux=np.array(spec) * u.erg / (u.cm * u.cm) / u.s / u.AA,
                                 spectral_axis = np.array(wLen) * u.AA)

    FluxCal = fluxcal.FluxCalibration()
    print('fluxCalibrate: dir(FluxCal) = ',dir(FluxCal))
    try:
        airmass = float(getHeaderValue(obsSpecFName,'AIRMASS'))
    except:
        airmass = float(getHeaderValue(obsSpecFName, 'SECZ'))
    print('fluxCalibrate:type(airmass) = ',type(airmass))
    FluxCal(objectSpectrum, airmass)

    import matplotlib.pyplot as plt
    from specreduce.calibration_data import AtmosphericExtinction, SUPPORTED_EXTINCTION_MODELS

    fig, ax = plt.subplots(2, 1, sharex=True)
    for model in SUPPORTED_EXTINCTION_MODELS:
        ext = AtmosphericExtinction(model=model)
        ax[0].plot(ext.spectral_axis, ext.extinction_mag, label=model)
        ax[1].plot(ext.spectral_axis, ext.transmission)
    ax[0].legend(fancybox=True, shadow=True)
    ax[1].set_xlabel("Wavelength ($\AA$)")
    ax[0].set_ylabel("Extinction (mag)")
    ax[1].set_ylabel("Transmission")
    plt.tight_layout()
    fig.show()

    from specreduce.calibration_data import AtmosphericTransmission
    fig, ax = plt.subplots()
    ext_default = AtmosphericTransmission()
    ext_custom = AtmosphericTransmission(data_file="/Users/azuri/entwicklung/python/data_reduction/specreduce/docs/atm_transmission_secz1.5_1.6mm.dat")
    ax.plot(ext_default.spectral_axis, ext_default.transmission, label=r"sec $z$ = 1; 1 mm H$_{2}$O", linewidth=1)
    ax.plot(ext_custom.spectral_axis, ext_custom.transmission, label=r"sec $z$ = 1.5; 1.6 mm H$_{2}$O", linewidth=1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, fancybox=True, shadow=True)
    ax.set_xlabel("Wavelength (microns)")
    ax.set_ylabel("Transmission")
    fig.show()

    from specreduce.calibration_data import load_MAST_calspec
    spec = load_MAST_calspec("agk_81d266_stisnic_007.fits")

    fig, ax = plt.subplots()
    ax.step(spec.spectral_axis, spec.flux, where="mid")
    ax.set_yscale('log')
    ax.set_xlabel(f'Wavelength ({spec.spectral_axis.unit})')
    ax.set_ylabel(f"Flux ({spec.flux.unit})")
    ax.set_title("AGK+81 266")
    fig.show()

    from specreduce.calibration_data import load_MAST_calspec, load_onedstds
    s1 = load_MAST_calspec("ltt9491_002.fits", remote=True)
    s2 = load_onedstds("snfactory", "LTT9491.dat")
    s3 = load_onedstds("eso", "ctiostan/ltt9491.dat")
    fig, ax = plt.subplots()
    ax.step(s1.spectral_axis, s1.flux, label="MAST", where="mid")
    ax.step(s2.spectral_axis, s2.flux, label="SNFactory", where="mid")
    ax.step(s3.spectral_axis, s3.flux, label="ESO", where="mid")
    ax.set_yscale('log')
    ax.set_xlabel(f"Wavelength ({s1.spectral_axis.unit})")
    ax.set_ylabel(f"Flux ({s1.flux.unit})")
    ax.set_title("LTT 9491")
    ax.legend()
    fig.show()


    from specreduce.calibration_data import get_reference_file_path
    kpno_extinction_file = get_reference_file_path("extinction/kpnoextinct.dat")
    print('fluxCalibrate:kpno_extinction_file = <'+kpno_extinction_file+'>')

    print('fluxCalibrate:dir(calibration_data) = ',dir(calibration_data))
    a=calibration_data.load_MAST_calspec(kpno_extinction_file)
    print('fluxCalibrate:a = ',a)
    b=calibration_data.load_onedstds()
    print('fluxCalibrate:b = ',b)


#    extinctionCurve = FluxCal.obs_extinction('kpno')
#    print('fluxCalibrate:extinctionCurve = ',extinctionCurve)

def continuum(spectrumFileNameIn, spectrumFileNameOut, fittingFunction, evalFunction, order, nIterReject, nIterFit, lowReject, highReject, type='difference', adjustSigLevels=False, useMean=False, display = False):
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
        markEmissionLines(xRange, yRange)
        plt.show()

    writeFits1D(spec,
                spectrumFileNameOut,
                wavelength=None,
                header=spectrumFileNameIn,
                CRVAL1=getHeaderValue(spectrumFileNameIn,'CRVAL1'),
                CRPIX1=getHeaderValue(spectrumFileNameIn,'CRPIX1'),
                CDELT1=getHeaderValue(spectrumFileNameIn,'CDELT1'),
               )

def scombine(fileListName, spectrumFileNameOut, method='median', lowReject=None, highReject=None, adjustSigLevels=False, useMean=False, display=False):
    fileNames = readFileToArr(fileListName)
    fileNames = [fileName if '/' in fileName else os.path.join(fileListName[:fileListName.rfind('/')],fileName) for fileName in fileNames]
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
        #print('scombine: notNaN = ',notNaN)
        if (len(goodPix[pix]) % 2 == 0) or (method == 'mean'):
            combinedSpectrum[pix] = np.mean(np.array(goodPix[pix][notNaN]))
        else:
            combinedSpectrum[pix] = np.median(np.array(goodPix[pix][notNaN]))
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

def removeFilesFromListWithAngleNotEqualTo(inputFileList,outputFileList,value):
    path = inputFileList[:inputFileList.rfind('/')]
    with open(inputFileList,'r') as f:
        lines = f.readlines()
    print('inputFileList contains ',len(lines),' files')
    lines = [line.strip() for line in lines]
    print('lines = ',lines)
    with open(outputFileList,'w') as f:
        for line in lines:
            if getHeaderValue(os.path.join(path,line), 'GR-ANGLE', hduNum=0) == value:
                f.write(line+'\n')
            else:
                print('removing file '+line+': '+getHeaderValue(os.path.join(path,line), 'OBJECT', hduNum=0)+': '+getHeaderValue(os.path.join(path,line), 'EXPTYPE', hduNum=0))

def fixStdWidth(stdFileName):
    with open(stdFileName,'r') as f:
        lines = f.readlines()
    with open(stdFileName+'_fixed','w') as f:
        f.write(lines[0])
        wLen = []
        flux = []
        for line in lines[1:]:
            wLen.append(float(line[:line.find('\t')]))
            flux.append(float(line[line.find('\t'):line.rfind('\t')].strip()))
        for i in range(len(wLen)):
            if i == 0:
                f.write('%.1f\t %.2f\t %.1f\n' % (wLen[i],flux[i],wLen[1]-wLen[0]))
            elif i == len(wLen)-1:
                f.write('%.1f\t %.2f\t %.1f\n' % (wLen[i],flux[i], (wLen[i]-wLen[i-1]) ))
            else:
                f.write('%.1f\t %.2f\t %.1f\n' % (wLen[i],flux[i], ((wLen[i]-wLen[i-1]) / 2.) + ((wLen[i+1]-wLen[i]) / 2.) ))

def plotSpec(fitsFileName):
    wLen = getWavelengthArr(fitsFileName)
    spec = getImageData(fitsFileName,0)
    plt.plot(wLen,spec)
    plt.show()

def separateByObsDate(fitslist):
    import shutil
    for fitsfile in fitslist:
        if os.path.isfile(fitsfile):
            path = fitsfile[:fitsfile.rfind('/')]
            #print('checking fitsfile <'+fitsfile+'>')
            header = getHeader(fitsfile,0)
            #print('header = ',header)
            obsDate = getDate(header['DATE-OBS'])
            print('fitsfile = ',fitsfile,': obsDate = ',obsDate,' str(obsDate) = <'+str(obsDate)+'>')
            newPath = os.path.join(path,str(obsDate))
            if not os.path.exists(newPath):
                os.mkdir(newPath)
            shutil.move(fitsfile,os.path.join(newPath,fitsfile[fitsfile.rfind('/')+1:]))

def fixDBSHeaders(filelist):
    if type(filelist) == type('str'):
        with open(filelist,'r') as f:
            filelist = f.readlines()
        filelist = [f.strip() for f in filelist]
    print('filelist = ',filelist)
    fluxStandardNames, fluxStandardDirs, fluxStandardFileNames = readFluxStandardsList()
    for fitsFile in filelist:
        if getHeaderValue(fitsFile,'IMAGETYP').lower() == 'object':
            if getHeaderValue(fitsFile,'OBJECT').lower() == 'arc':
                print('IMAGETYP = ',getHeaderValue(fitsFile,'IMAGETYP'))
                setHeaderValue(fitsFile,'IMAGETYP','ARC')
                print('new IMAGETYP = ',getHeaderValue(fitsFile,'IMAGETYP'))
        if getHeaderValue(fitsFile,'IMAGETYP').lower() == 'flat':
            if getHeaderValue(fitsFile,'OBJECT').lower() == 'sky':
                print('IMAGETYP = ',getHeaderValue(fitsFile,'IMAGETYP'))
                setHeaderValue(fitsFile,'IMAGETYP','SKYFLAT')
                print('new IMAGETYP = ',getHeaderValue(fitsFile,'IMAGETYP'))
        if getHeaderValue(fitsFile,'IMAGETYP').lower() == 'flat':
            if getHeaderValue(fitsFile,'OBJECT').lower() == 'internal_flat':
                print('OBJECT = ',getHeaderValue(fitsFile,'OBJECT'))
                setHeaderValue(fitsFile,'OBJECT','DOMEFLAT')
                print('new OBJECT = ',getHeaderValue(fitsFile,'OBJECT'))
        if getHeaderValue(fitsFile,'IMAGETYP').lower() == 'object':
            if getHeaderValue(fitsFile,'OBJECT').lower() in fluxStandardNames:
                print('IMAGETYP = ',getHeaderValue(fitsFile,'IMAGETYP'))
                setHeaderValue(fitsFile,'IMAGETYP','FLUXSTDS')
                print('new IMAGETYP = ',getHeaderValue(fitsFile,'IMAGETYP'))
        if getHeaderValue(fitsFile,'IMAGETYP').lower() == 'domeflat':
            if getHeaderValue(fitsFile,'OBJECT').lower() == 'internal_flat':
                setHeaderValue(fitsFile,'IMAGETYP','FLAT')
                setHeaderValue(fitsFile,'OBJECT','DOMEFLAT')
        if getHeaderValue(fitsFile,'IMAGETYP').lower() == 'skyflat':
            if getHeaderValue(fitsFile,'OBJECT').lower() == 'sky':
                setHeaderValue(fitsFile,'IMAGETYP','FLAT')
                setHeaderValue(fitsFile,'OBJECT','SKYFLAT')


def open_image(imagename):
    from astropy.io import fits
    from astropy.wcs import WCS
    import montage_wrapper as montage
    hdu = fits.open(imagename)
    hdu = montage.reproject_hdu(hdu[0], north_aligned=True)
    image = hdu.data
    nans = np.isnan(image)
    image[nans] = 0
    header = hdu.header
    wcs = WCS(header)
    return image, header, wcs

def createFindingChartFromFits(fitsFileName,
                               widthInArcSeconds,
                               pnMajDiamInArcSeconds,
                               objectName,
                               ra,
                               dec,
                               outDirName,
                               display=False):
    from myUtils import angularDistanceFromXY#getArcsecDistance
    from astropy.nddata import Cutout2D
    from astropy import units as u
    import astropy.visualization as vis
    from astropy.wcs import WCS
    import montage_wrapper as montage
#    from myUtils import getRaDecFromXY

#    image, header, wcs = open_image(fitsFileName)
#    print('image.shape = ',image.shape)

#    fig,ax = plt.subplots(1,1, subplot_kw=dict(projection=wcs))
#    ra = ax.coords[0]
#    dec = ax.coords[1]
#    ra.set_major_formatter('hh:mm:ss.s')
#    dec.set_major_formatter('dd:mm:ss')

    #setting image scale
#    interval = vis.PercentileInterval(99.9)
#    vmin,vmax = interval.get_limits(image)
#    norm = vis.ImageNormalize(vmin=vmin, vmax=vmax, stretch=vis.LogStretch(1000))
#    ax.imshow(image, cmap =plt.cm.Reds, norm = norm, origin = 'lower')
#    ax.set_ylabel('Dec.')
#    ax.set_xlabel('RA')
#    plt.show()

    print('fitsFileName = ',fitsFileName[fitsFileName.rfind('/')+1:])
    fitsData = getImageData(fitsFileName,0)
#    plt.imshow(fitsData)
#    plt.show()
    fitsDataShape = fitsData.shape
    print('fitsDataShape = ',fitsDataShape)
#    x0 = 0
#    y0 = 0
#    x1 = fitsDataShape[0]-1
#    y1 = fitsDataShape[1]-1
#    dist = angularDistanceFromXY(fitsFileName, 0, 0, x1, y1)
#    print('x0 = ',x0,', y0 = ',y0,', x1 = ',x1,', y1 = ',y1,': dist = ',dist)

    x0 = 0
    y0 = 0
    x1 = 0
    y1 = fitsDataShape[1]-1
    dist = angularDistanceFromXY(fitsFileName, 0, 0, x1, y1)
    print('x0 = ',x0,', y0 = ',y0,', x1 = ',x1,', y1 = ',y1,': dist = ',dist)

    x0 = 0
    y0 = 0
    x1 = fitsDataShape[0]-1
    y1 = 0
    dist = angularDistanceFromXY(fitsFileName, 0, 0, x1, y1)
    print('x0 = ',x0,', y0 = ',y0,', x1 = ',x1,', y1 = ',y1,': dist = ',dist)

    center = [fitsDataShape[0]/2.,fitsDataShape[1]/2.]
    newWidthInPixels = int(fitsDataShape[0] * widthInArcSeconds / dist)
    print('new width in pixels = ',newWidthInPixels)
    pixPerArcSec = newWidthInPixels / widthInArcSeconds
    print('pixPerArcSec = ',pixPerArcSec)
    circleRadius = pnMajDiamInArcSeconds * pixPerArcSec / 2.
    print('circleRadius = ',circleRadius,' pixels')
    if np.min([fitsDataShape[0],fitsDataShape[1]]) < newWidthInPixels:
        cutout = np.full((newWidthInPixels,newWidthInPixels),np.min(fitsData))
        x0 = np.max([int(newWidthInPixels/2 - fitsDataShape[0]/2),0])
        x1 = np.min([int(newWidthInPixels/2 + fitsDataShape[0]/2),newWidthInPixels])
        y0 = np.max([int(newWidthInPixels/2 - fitsDataShape[1]/2),0])
        y1 = np.min([int(newWidthInPixels/2 + fitsDataShape[1]/2),newWidthInPixels])
        print('x0 = ',x0,', x1 = ',x1,', y0 = ',y0,', y1 = ',y1)
        print('x1-x0 = ',x1-x0)
        print('y1-y0 = ',y1-y0)
        dx0 = np.max([0,int(fitsData.shape[0]/2-newWidthInPixels/2)])
        dx1 = np.min([newWidthInPixels,fitsData.shape[0]])
        dy0 = np.max([0,int(fitsData.shape[1]/2-newWidthInPixels/2)])
        dy1 = np.min([newWidthInPixels,fitsData.shape[1]])
        print('dx0 = ',dx0,', dx1 = ',dx1,', dy0 = ',dy0,', dy1 = ',dy1)
        print('dx1-dx0 = ',dx1-dx0)
        print('dy1-dy0 = ',dy1-dy0)
        cutout[x0:x1,y0:y1] = fitsData[dx0:dx1,dy0:dy1]
    else:
        position = (center[0], center[1])
        size = newWidthInPixels * u.pixel
        print('position = ',position,', size = ',size)
        cutout = Cutout2D(fitsData,position,size)
    plt.gray()
    if fitsFileName.rfind('_shs') < 0:
        tmp = fitsFileName[:fitsFileName.rfind('_Ha')]
        outFileName = os.path.join(outDirName,fitsFileName[fitsFileName.rfind('/')+1:fitsFileName.rfind('_Ha')]+'_findingChart.png')
    else:
        tmp = fitsFileName[:fitsFileName.rfind('_shs')]
        outFileName = os.path.join(outDirName,fitsFileName[fitsFileName.rfind('/')+1:fitsFileName.rfind('_shs')]+'_findingChart.png')
    idPNMain = tmp[tmp.rfind('_')+1:]
    #ra, dec = getRaDecFromXY(fitsFileName,fitsDataShape[0]/2,fitsDataShape[1]/2)
    print('idPNMain = ',idPNMain)
    minVal = np.min(cutout.data)
    maxVal = np.max(cutout.data)
    vmax = minVal+(maxVal-minVal)/7.
    print('minVal = ',minVal,', maxVal = ',maxVal,', vmax = ',vmax)
    interval = vis.PercentileInterval(99.95)
    vmin,vmax = interval.get_limits(cutout.data)
    norm = vis.ImageNormalize(vmin=vmin, vmax=vmax, stretch=vis.LogStretch(1000))

    plt.imshow(cutout.data, cmap=plt.cm.binary, origin='lower', norm=norm)#, vmin = minVal, vmax = vmax)#, cmap='gray', vmin=0, vmax=255)norm=norm)#
    #plt.plot([newWidthInPixels/2,newWidthInPixels/2+newWidthInPixels/20],[newWidthInPixels/2,newWidthInPixels/2-newWidthInPixels/20],'r-')
    #plt.plot([newWidthInPixels/2,newWidthInPixels/2],[newWidthInPixels/2,newWidthInPixels/2-newWidthInPixels/20],'r-')
    #plt.plot([newWidthInPixels/2,newWidthInPixels/2+newWidthInPixels/20],[newWidthInPixels/2,newWidthInPixels/2],'r-')
    circle1 = plt.Circle((int(newWidthInPixels/2), int(newWidthInPixels/2)), circleRadius, color='r', fill=False)
    plt.gca().add_patch(circle1)
    plt.title('HASH ID '+idPNMain+' '+objectName)
#    cutout.plot_on_original(color='white')
    plt.xticks([])
    plt.yticks([])
    plt.xlabel('RA %s +/- %d arcsec' % (ra, widthInArcSeconds/2))
    plt.ylabel('DEC %s +/- %d arcsec' % (dec, widthInArcSeconds/2))
    if display:
        plt.show()
    else:
        plt.savefig(outFileName)
        plt.close()
#    STOP
    return outFileName


def createFindingChartsFromFits(dirName,widthInArcSeconds,pnMainWithMajDiamFileName):
    import csvFree,csvData
    pnMain = csvFree.readCSVFile(pnMainWithMajDiamFileName)
    fileList = os.listdir(dirName)
    for item in fileList:
        if item.endswith('.fits'):
            tmp = item[:item.rfind('_shs')]
            idPNMain = tmp[tmp.rfind('_')+1:]
            found = pnMain.find('idPNMain',idPNMain)[0]
            majDiam = float(pnMain.getData('MajDiam',found))
            name = pnMain.getData('Name',found)
            print('name = <'+name+'>: MajDiam = ',majDiam)
            createFindingChartFromFits(os.path.join(dirName,item),widthInArcSeconds,majDiam,name)
