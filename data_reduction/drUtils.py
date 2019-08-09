#from astropy import units as u
#from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
from astropy import units as u
import astropy.modeling.tests.irafutil as iu
from astropy.nddata import CCDData
import ccdproc# import Combiner, subtract_overscan
from collections import namedtuple
import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.legendre import legval
#import hammer
#import numpy as np
import os
from scipy.interpolate import interp1d
#from scipy import interpolate
from scipy.interpolate import griddata
from scipy.interpolate import UnivariateSpline
from scipy.signal import medfilt
from shutil import copyfile
#from sklearn import linear_model

# TODO: Make it possible to pass in CCDData instead of fits file names

Info = namedtuple('Info', 'start height')

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

def writeFits(ccdData, output_file, overwrite=True):
    hdulist = ccdData.to_hdu()
    hdulist.writeto(output_file, overwrite=overwrite)

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
def separateFileList(inList, suffixes, exptypes=None, objects=None, changeNames=False):
    def createLists(exptype, object, lines, suffix):
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
            hdulist = pyfits.open(line)

            fOutName = line
#            fOutName = addSuffixToFileName(fOutName, suffix)

            expType = hdulist[0].header['EXPTYPE']
            objectName = hdulist[0].header['OBJECT']
            if changeNames:
                fOutName = os.path.join(fOutName[:fOutName.rfind('/')],
                                        expType+'_'+objectName+'_'+fOutName[fOutName.rfind('/')+1:])
                if suffix == '':
                    copyfile(line, fOutName)
                line = fOutName
            if expType == exptype:
                if object == '*':
                    print('exptype = <'+exptype+'>, object = <'+object+'>: expType = <'+expType+'>, objectName = <'+objectName+'>: adding <'+line+'> to isNames')
                    isNames.append(line)
                else:
                    if objectName == object:
                        print('exptype = <'+exptype+'>, object = <'+object+'>: expType = <'+expType+'>, objectName = <'+objectName+'>: adding <'+line+'> to isNames')
                        isNames.append(line)
                    else:
                        print('exptype = <'+exptype+'>, object = <'+object+'>: expType = <'+expType+'>, objectName = <'+objectName+'>: adding <'+line+'> to isntNames')
                        isntNames.append(line)
            else:
                print('exptype = <'+exptype+'>, object = <'+object+'>: expType = <'+expType+'>, objectName = <'+objectName+'>: adding <'+line+'> to isntNames')
                isntNames.append(line)
        with open(isList,'w') as f:
            for name in isNames:
                f.write(addSuffixToFileName(name,suffix)+'\n')
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
                print('line = ',line,': expType = ',expType)
                if expType not in [x['expType'] for x in expTypes]:
                    listOutName = os.path.join(path,expType.lower()+'.list')
                    listOutName = addSuffixToFileName(listOutName,suffix)
                    print('writing to listOutName = <'+listOutName+'>')
                    expTypes.append({'expType':expType, 'listOutName':listOutName})
                    silentRemove(listOutName)
                else:
                    for x in expTypes:
                        if x['expType'] == expType:
                            listOutName = x['listOutName']
                print('listOutName = <'+listOutName+'>')

                objectName = hdulist[0].header['OBJECT']
                print('line = ',line,': objectName = ',objectName,', expType = ',expType)
                if expType == 'SCIENCE':
                    if objectName not in [x['objectName'] for x in objectNames]:
                        listOutNameObs = os.path.join(path,objectName.lower()+'.list')
                        listOutNameObs = addSuffixToFileName(listOutNameObs, suffix)
                        print('writing to listOutNameObs = <'+listOutNameObs+'>')
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

                if expType == 'SCIENCE':
                    with open(listOutNameObs,'a') as f:
                        f.write(fOutName+'\n')

                if objectName != 'Bias':
                    with open(nonZerosOutName,'a') as f:
                        f.write(fOutName+'\n')

                if expType not in ['FLAT','BIAS']:
                    with open(nonFlatsOutName,'a') as f:
                        f.write(fOutName+'\n')
        else:
            individualLists = []
            if len(exptypes) != len(objects):
                print('separateFileList: ERROR: lengths of <exptypes> and <objects> are not the same')
                STOP
            for iExptype in np.arange(0,len(exptypes),1):
                exptype = exptypes[iExptype]
                for object in objects[iExptype]:
                    if  object == 'individual':
                        for line in lines:
                            hdulist = pyfits.open(line)
                            expType = hdulist[0].header['EXPTYPE']
                            if expType == exptype:
                                objectName = hdulist[0].header['OBJECT']
                                if objectName not in individualLists:
                                    individualLists.append(objectName)
                    else:
                        createLists(exptype,object,lines,suffix)

                if len(individualLists) > 0:
                    for object in individualLists:
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
        print('applying minmax clipping')
        combiner.minmax_clipping(min_clip=clippingParameters['min_clip'],
                                 max_clip=clippingParameters['max_clip'])

    elif clippingMethod == 'sigma':
        print('applying sigma clipping')
        #print('dir(clippingParameters = ',dir(clippingParameters))
        print('clippingParameters.keys() = ',clippingParameters.keys())
        print("clippingParameters['low_thresh'] = ",clippingParameters['low_thresh'])
        print("clippingParameters['high_thresh'] = ",clippingParameters['high_thresh'])
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
        print('applying extrema clipping')
        combiner.clip_extrema(nlow=clippingParameters['nlow'],
                              nhigh=clippingParameters['nhigh'])

    if scaling:
        print('applying image scaling')
        scaling_func = lambda arr: 1/np.ma.average(arr)
        combiner.scaling = scaling_func

    combinedImage = None
    if combinerMethod == 'average':
        combinedImage = combiner.average_combine()
    elif combinerMethod == 'median':
        combinedImage = combiner.median_combine()

    # replace all values less than <minVal> with <minVal>
    if minVal is not None:
        combinedImage.data = np.minimum(combinedImage.data, minVal)

    if fitsOutName is not None:
        print('combinedImage.data.shape = ',combinedImage.data.shape)
        print('combinedImage.header = ',combinedImage.header)
        writeFits(combinedImage, fitsOutName, overwrite=overwrite)

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
        ccdData = CCDData.read(fitsFilesIn[iFile], unit="adu")
        ccdDataNoOverscan = ccdproc.subtract_overscan(ccdData, fits_section=overscanSection)
        if trimSection is None:
            if fitsFilesOut is not None:
                ccdDataNoOverscan.write(fitsFilesOut[iFile], overwrite=overwrite)
            dataOut.append(ccdDataNoOverscan)
        else:
            trimmed = ccdproc.trim_image(ccdDataNoOverscan, fits_section=trimSection)
            if fitsFilesOut is not None:
                writeFits(trimmed, fitsFilesOut[iFile], overwrite=overwrite)
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
        print('inFile = ',fitsFilesIn[iFile],', masterBias = ',masterBias)
        ccdDataBiasSubtracted = ccdproc.subtract_bias(ccdData,
                                                      masterBiasArr)
        dataOut.append(ccdDataBiasSubtracted)
        if fitsFilesOut is not None:
            writeFits(ccdDataBiasSubtracted, fitsFilesOut[iFile], overwrite=overwrite)
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
#                                                    fsmode == ‘convolve’ and psfk is None. The current choices are
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
            ccdDataCleaned = ccdproc.cosmicray_lacosmic(ccdData,
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
            ccdDataCleaned = ccdproc.cosmicray_median(ccdData,
                                                      mbox=mbox,
                                                      rbox=rbox,
                                                      gbox=gbox,
                                                      thresh=thresh,
                                                      error_image=error_image)
        ccdDataCleaned.header['COSMIC'] = ctype.upper()
        dataOut.append(ccdDataCleaned)
        if fitsFilesOut is not None:
            writeFits(ccdDataCleaned, fitsFilesOut[iFile], overwrite=overwrite)
    return dataOut

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
    print('flatCorrect: ccdDataFlat = ',ccdDataFlat)
    print('min(ccdDataFlat) = ',np.min(ccdDataFlat),', max(ccdDataFlat) = ',np.max(ccdDataFlat))
    print('ccdDataFlat[np.where(np.isnan(ccdDataFlat))] = ',ccdDataFlat[np.where(np.isnan(ccdDataFlat))])
    print('ccdDataFlat[np.where(np.isinf(ccdDataFlat))] = ',ccdDataFlat[np.where(np.isinf(ccdDataFlat))])
    dataOut = []

    for iFile in np.arange(0,len(fitsFilesIn),1):
        ccdData = CCDData.read(fitsFilesIn[iFile], unit="adu")
        print('flatCorrect: iFile = ',iFile,': ccdData = ',ccdData)
        print('min(ccdData) = ',np.min(ccdData),', max(ccdData) = ',np.max(ccdData))
        print('ccdData[np.where(np.isnan(ccdData))] = ',ccdData[np.where(np.isnan(ccdData))])
        print('ccdData[np.where(np.isinf(ccdData))] = ',ccdData[np.where(np.isinf(ccdData))])
        ccdDataFlattened = ccdproc.flat_correct(ccdData,
                                                flat=ccdDataFlat,
                                                min_value=min_value,
                                                norm_value=norm_value,
                                                add_keyword=add_keyword)
        print('flatCorrect: iFile = ',iFile,' ',fitsFilesIn[iFile],': ccdDataFlattened = ',ccdDataFlattened)
        dataOut.append(ccdDataFlattened)
        if fitsFilesOut is not None:
            writeFits(ccdDataFlattened, fitsFilesOut[iFile], overwrite=overwrite)
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
        print('z = ',z)
        y = []
        for x in xNorm:
            z[1] = x
            for i in np.arange(2,len(coeffs)):
                z[i] = (2.0 * x * z[i-1]) - z[i-2]
            yTemp = (coeffs[0] * z[0]) + (coeffs[1] * z[1])
            for i in np.arange(2,len(coeffs),1):
                yTemp += coeffs[i] * z[i]
            y.append(yTemp)
        print('xNorm = ',len(xNorm),': ',xNorm)
        print('y = ',len(y),': ',y)
        return y

    yCheck = chebval(xNorm, coeffs)
#    print('yCheck = ',yCheck)
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
        print('z = ',z)
        y = []
        for x in xNorm:
            z[1] = x
            for i in np.arange(2,len(coeffs)):
                z[i] = ((((2.0*(i+1.0))-3.0) * x * z[i-1]) - ((i-1.0) * z[i-2])) / (i)
            print('z = ',z)
            yTemp = (coeffs[0] * z[0]) + (coeffs[1] * z[1])
            for i in np.arange(2,len(coeffs),1):
                yTemp += coeffs[i] * z[i]
                print('i = ',i,': yTemp = ',yTemp)
            y.append(yTemp)
        print('xNorm = ',len(xNorm),': ',xNorm)
        print('y = ',len(y),': ',y)
        return y

    yCheck = legval(xNorm,coeffs)
#    print('yCheck = ',yCheck)
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
        print('x = ',x,': s = ',s,', j = ',j,', a = ',a,', b = ',b)
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
        print('x = ',x,': s = ',s,', j = ',j,', a = ',a,', b = ',b,', z_0 = ',z_0,', z_1 = ',z_1,', z_2 = ',z_2,', z_3 = ',z_3)
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
def calcTrace(dbFile, apNum=0, xRange = None):
    records = iu.get_records(dbFile)
    xCenter, yCenter = [float(c) for c in records[apNum].get_fields()['center'].split(' ')]
    curve = records[apNum].get_fields()['curve']
    funcs = ['none','chebyshev','legendre','cubicSpline','linearSpline']
    function = funcs[int(curve[0])]
    order = curve[1][0]
    if xRange is None:
        xRange = [curve[2][0],curve[3][0]]
    xNorm = []
    xArr = np.arange(xRange[0],xRange[1]+1)
    for x in xArr:
        if (function == funcs[1]) or (function == funcs[2]):
            xNorm.append((2.0 * x - (xRange[1] + xRange[0])) / (xRange[1] - xRange[0]))

    coeffs = [c[0] for c in curve[4:]]

    y = None
    if function == funcs[1]:
        y = chebyshev(xNorm, coeffs) + xCenter
    elif function == funcs[2]:
        y = legendre(xNorm, coeffs) + xCenter
    elif function == funcs[3]:
        y = cubicSpline(xRange, order, coeffs) + xCenter
    elif function == funcs[4]:
        y = linearSpline(xRange, order, coeffs) + xCenter
    else:
        print('calcTrace: could not identify function <'+function+'>')

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
    print('image = ',image)
    print('trace = ',len(trace),': ',trace)
    print('markCenter: trace[0].shape = ',trace[0].shape)
    for i in np.arange(0,trace[0].shape[0],1):
        image.data[int(trace[0][i]), int(trace[1][i])] = 0.
    if imFileOut is not None:
        writeFits(image, imFileOut, overwrite=True)
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

# NOTE that the horizontal trace needs to come from an image that was rotated by
# 90 degrees and flipped along the long axis
def interpolateTraceIm(imFiles, dbFileVerticalTrace, dbFileHorizontalTrace, markCenters=False):
    # read imFile for dimensions

    image = CCDData.read(imFiles[0], unit="adu")

    records = iu.get_records(dbFileVerticalTrace)
    verticalTraces = []
    for i in np.arange(0,len(records),1):
        verticalTraces.append(calcTrace(dbFileVerticalTrace,i,[1,image.shape[0]]))
#        x, y = verticalTraces[len(verticalTraces)-1]

    records = iu.get_records(dbFileHorizontalTrace)
    horizontalTraces = []
    for i in np.arange(0,len(records),1):
        horizontalTraces.append(calcTrace(dbFileHorizontalTrace,i,[1,image.shape[1]]))
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
        for ix in np.arange(0,image.data.shape[0],1):
            for iy in np.arange(0,image.data.shape[1],1):
                zOrig[(ix*image.data.shape[1]) + iy] = image.data[ix,iy]
        zFit = griddata(xyOrig, zOrig, xyFit, method='cubic')

        # re-order vector back to 2D array
        for ix in np.arange(0,image.data.shape[0],1):
            for iy in np.arange(0,image.data.shape[1],1):
                image.data[ix, iy] = zFit[(ix*image.data.shape[1]) + iy]

        #trim image to only contain good data inside the original trace
        maxSizeArr = np.zeros(image.data.shape)
        maxSizeArr[np.where(np.isnan(image.data))] = 1
        tempArr = max_size(maxSizeArr)
        image.data = image.data[tempArr[1][0]:tempArr[1][0]+tempArr[0][0],tempArr[1][1]:tempArr[1][1]+tempArr[0][1]]
        image.mask = image.mask[tempArr[1][0]:tempArr[1][0]+tempArr[0][0],tempArr[1][1]:tempArr[1][1]+tempArr[0][1]]
        image.uncertainty = image.uncertainty[tempArr[1][0]:tempArr[1][0]+tempArr[0][0],tempArr[1][1]:tempArr[1][1]+tempArr[0][1]]

        image.header['NAXIS1'] = tempArr[0][1]
        image.header['NAXIS2'] = tempArr[0][0]

        writeFits(image, imFile[:imFile.rfind('.')]+'i.fits', overwrite=True)
    return 1

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
        profileImage[row,:] = medfilt(profileImage[row,:], rowMedianSmoothSize)
    print('makeSkyFlat: profileImage = ',profileImage.shape,': ',profileImage)
    print('makeSkyFlat: 2) profileImage[np.where(np.isnan(profileImage))] = ',profileImage[np.where(np.isnan(profileImage))])
    print('makeSkyFlat: 2) profileImage[np.where(np.isinf(profileImage))] = ',profileImage[np.where(np.isinf(profileImage))])

    # set values <= to 0.01
    profileImage[np.where(profileImage <= 0.00001)] = 0.001
    print('makeSkyFlat: profileImage > 0. = ',profileImage.shape,': ',profileImage)
    print('makeSkyFlat: profileImage[np.where(np.isnan(profileImage))] = ',profileImage[np.where(np.isnan(profileImage))])
    print('makeSkyFlat: profileImage[np.where(np.isinf(profileImage))] = ',profileImage[np.where(np.isinf(profileImage))])

    image.data = profileImage
    writeFits(image, skyFlatOut, overwrite=True)
