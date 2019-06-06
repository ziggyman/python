#from astropy import units as u
#from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
from astropy import units as u
import astropy.modeling.tests.irafutil as iu
from astropy.nddata import CCDData
import ccdproc# import Combiner, subtract_overscan
import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.legendre import legval
#import hammer
#import numpy as np
import os
#from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.interpolate import griddata
from shutil import copyfile
#from sklearn import linear_model

# TODO: Make it possible to pass in CCDData instead of fits file names

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
# <changeNames>: if True copy each <file> in <inList> to <EXPTYPE_OBJECT_><file>
def separateFileList(inList, suffixes, changeNames=False):
    lines = []
    path = inList[:inList.rfind('/')]
    with open(inList,'r') as f:
        lines = [os.path.join(path,line.strip('\n')) for line in f]

    for suffix in suffixes:
        listOutName = ''
        listOutNameObs = ''

        objectNames = []
        expTypes = []

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

    if fitsOutName is not None:
        combinedImage.write(fitsOutName, overwrite=overwrite)

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
                trimmed.write(fitsFilesOut[iFile], overwrite=overwrite)
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
            ccdDataBiasSubtracted.write(fitsFilesOut[iFile], overwrite=overwrite)
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
#                                               psfsize : int, optional
#                                               Size of the kernel to calculate. Returned kernel will have size psfsize
#                                               x psfsize. psfsize should be odd. Default: 7.
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
            ccdDataCleaned.write(fitsFilesOut[iFile], overwrite=overwrite)
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
    dataOut = []

    for iFile in np.arange(0,len(fitsFilesIn),1):
        ccdData = CCDData.read(fitsFilesIn[iFile], unit="adu")
        ccdDataFlattened = ccdproc.flat_correct(ccdData,
                                                flat=ccdDataFlat,
                                                min_value=min_value,
                                                norm_value=norm_value,
                                                add_keyword=add_keyword)
        dataOut.append(ccdDataFlattened)
        if fitsFilesOut is not None:
            ccdDataFlattened.write(fitsFilesOut[iFile], overwrite=overwrite)
    return dataOut

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

def linearSpline(x, xRange, order, coeffs):
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

def cubicSpline(x, xRange, order, coeffs):
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

#xRange: [1,size]
#return x=[0...size-1], y=[0...size-1]
def calcTrace(dbFile, apNum=0, xRange = None):
    records = iu.get_records(dbFile)
#    print('calcTrace: dbFile = <'+dbFile+'>: len(records) = ',len(records))
    xCenter, yCenter = [float(c) for c in records[apNum].get_fields()['center'].split(' ')]
#    print('xCenter = ',xCenter,', yCenter = ',yCenter)
    curve = records[apNum].get_fields()['curve']
    nCoeffs = len(curve)
#    print('calcTrace: nCoeffs = ',nCoeffs)
    funcs = ['none','chebyshev','legendre','cubicSpline','linearSpline']
    function = funcs[int(curve[0])]
#    print('calcTrace: function = ',function)
    order = curve[1][0]
#    print('calcTrace: order = ',order)
    if xRange is None:
        xRange = [curve[2][0],curve[3][0]]
#    print('calcTrace: xRange = ',xRange)
    xNorm = []
    xArr = np.arange(xRange[0],xRange[1]+1)
    for x in xArr:
        if (function == funcs[1]) or (function == funcs[2]):
            xNorm.append((2.0 * x - (xRange[1] + xRange[0])) / (xRange[1] - xRange[0]))

#    print('calcTrace: xNorm = ',xNorm)
    coeffs = [c[0] for c in curve[4:]]
#    print('calcTrace: coeffs = ',coeffs)

    y = None
    if function == funcs[1]:
        y = chebyshev(xNorm, coeffs) + xCenter
    elif function == funcs[2]:
        y = legendre(xNorm, coeffs) + xCenter
    elif function == funcs[3]:
        y = cubicSpline(np.arange(xRange[0],xRange[1]+1), xRange, order, coeffs) + xCenter
    elif function == funcs[4]:
        y = linearSpline(np.arange(xRange[0],xRange[1]+1), xRange, order, coeffs) + xCenter
    else:
        print('calcTrace: could not identify function <'+function+'>')
#    print('calcTrace: function = <'+function+'>: y = ',y)

    return [xArr-1.0,y-1.0]

# imFile: string
# trace: [xArr,yArr]
def markCenter(imFileIn, trace, imFileOut=None):
    image = CCDData.read(imFileIn, unit="adu")
    print('image = ',image)
    print('trace = ',len(trace),': ',trace)
    print('markCenter: trace[0].shape = ',trace[0].shape)
    for i in np.arange(0,trace[0].shape[0],1):
#        print('markCenter: trace[0][',i,'] = ',trace[0][i])
#        tempIm = image[int(trace[0][i])]
        print('markCenter: int(trace[0][',i,']) = ',int(trace[0][i]),', int(trace[1][',i,']) = ',int(trace[1][i]))
#        print('markCenter: trace[1][',i,'] = ',trace[1][i])
#        print('markCenter: image[trace[0][',i,'],trace[1][',i,']] = ',image[int(trace[0][i]),int(trace[1][i])])
#        print('markCenter: image[',int(trace[0][i]),', ',int(trace[1][i]),'] = ',image[int(trace[0][i]),int(trace[1][i])])
#        print('markCenter: setting [',trace[0][i],', ',trace[1][i],'] to 0')
        image.data[int(trace[0][i]), int(trace[1][i])] = 0.
#        print('markCenter: image[',int(trace[0][i]),', ',int(trace[1][i]),'] = ',image[int(trace[0][i]),int(trace[1][i])])
    if imFileOut is not None:
        image.write(imFileOut, overwrite=True)
    return image

#def interpolatePixel(image, x, y):

# NOTE that the horizontal trace needs to come from an image that was rotated by
# 90 degrees and flipped along the long axis
def interpolateTraceIm(imFile, dbFileVerticalTrace, dbFileHorizontalTrace):
    # read imFile for dimensions
    image = CCDData.read(imFile, unit="adu")
    print('image.shape = ',image.shape)

    records = iu.get_records(dbFileVerticalTrace)
    print('len(records) = ',len(records))
    verticalTraces = []
    for i in np.arange(0,len(records),1):
        verticalTraces.append(calcTrace(dbFileVerticalTrace,i,[1,image.shape[0]]))
        x, y = verticalTraces[len(verticalTraces)-1]
        print('vertical trace ',len(verticalTraces)-1,': x.shape = ',x.shape,', y.shape = ',y.shape)
        for i in np.arange(0,x.shape[0],1):
            print('x[',i,'] = ',x[i],', y[',i,'] = ',y[i])

    records = iu.get_records(dbFileHorizontalTrace)
    print('len(records) = ',len(records))
    horizontalTraces = []
    for i in np.arange(0,len(records),1):
        horizontalTraces.append(calcTrace(dbFileHorizontalTrace,i,[1,image.shape[1]]))
        x, y = horizontalTraces[len(horizontalTraces)-1]
        print('horizontal trace ',len(horizontalTraces)-1,': x.shape = ',x.shape,', y.shape = ',y.shape)
        for i in np.arange(0,x.shape[0],1):
            print('x[',i,'] = ',x[i],', y[',i,'] = ',y[i])

    image = CCDData.read(imFile, unit="adu")
    print('image.shape = ',image.shape)
    inFile = ''
    outFile = ''
    for i in np.arange(0,len(verticalTraces),1):
        if i == 0:
            inFile = imFile
            outFile = imFile[:imFile.rfind('.')]+'_vCenter%dMarked.fits' % (i)
        else:
            inFile = outFile
            outFile = imFile[:imFile.rfind('.')]+'_vCenter%dMarked.fits' % (i)
        print('marking vertical trace ',i)
        markCenter(inFile,verticalTraces[i],outFile)
    horFile = outFile
    for i in np.arange(0,len(horizontalTraces),1):
        if i == 0:
            inFile = horFile
            outFile = inFile[:inFile.rfind('.')]+'_hCenter%dMarked.fits' % (i)
        else:
            inFile = outFile
            outFile = horFile[:horFile.rfind('.')]+'_hCenter%dMarked.fits' % (i)
        print('marking horizontal trace ',i)
        markCenter(inFile,[horizontalTraces[i][1],horizontalTraces[i][0]],outFile)

    print('len(verticalTraces) = ',len(verticalTraces))
    print('len(verticalTraces[0]) = ',len(verticalTraces[0]))
    print('verticalTraces[0][0].shape = ',verticalTraces[0][0].shape)
    print('horizontalTraces[0][1][:] - horizontalTraces[0][1][0] = ',horizontalTraces[0][1][:] - horizontalTraces[0][1][0])
    print('verticalTraces[0][1][:] - verticalTraces[0][1][0] = ',verticalTraces[0][1][:] - verticalTraces[0][1][0])
    coordsFit = np.ndarray(shape=(image.shape[0],image.shape[1],2), dtype=np.float32)
    for i in np.arange(0,image.shape[0],1):
        coordsFit[i,:,0] = float(i) + horizontalTraces[0][1][:] - horizontalTraces[0][1][0]
    for i in np.arange(0,image.shape[1],1):
        coordsFit[:,i,1] = float(i) + verticalTraces[0][1][:] - verticalTraces[0][1][0]

    xx = np.arange(0,verticalTraces[0][1].shape[0])
    print('xx = ',xx.shape,': ',xx)
    yy = np.arange(0,horizontalTraces[0][1].shape[0],1)
    print('yy = ',yy.shape,': ',yy)
    f = interpolate.interp2d(yy,xx, image, kind='linear')
    print('interpolated function = ',f)
    print('horizontalTraces[0][1] = ',horizontalTraces[0][1].shape,': ',horizontalTraces[0][1])
    print('verticalTraces[0][1] = ',verticalTraces[0][1].shape,': ',verticalTraces[0][1])

    xyOrig = np.ndarray(shape=(image.data.shape[0] * image.data.shape[1],2), dtype=np.float32)
    xyFit = np.ndarray(shape=(image.data.shape[0] * image.data.shape[1],2), dtype=np.float32)
    zOrig = np.ndarray(shape=(image.data.shape[0] * image.data.shape[1]), dtype=np.float32)
    print('image.data.shape = ',image.data.shape)
    nPoints = 0
    for ix in np.arange(0,image.data.shape[0],1):
        for iy in np.arange(0,image.data.shape[1],1):
            xyOrig[(ix*image.data.shape[1]) + iy,0] = ix
            xyOrig[(ix*image.data.shape[1]) + iy,1] = iy
            zOrig[(ix*image.data.shape[1]) + iy] = image.data[ix,iy]
            xyFit[(ix*image.data.shape[1]) + iy,0] = ix + horizontalTraces[0][1][iy] - horizontalTraces[0][1][0]
            xyFit[(ix*image.data.shape[1]) + iy,1] = iy + verticalTraces[0][1][ix] - verticalTraces[0][1][0]
            print('xOrig[',(ix*image.data.shape[1]) + iy,'] = ',ix,', yOrig[',(ix*image.data.shape[1]) + iy,'] = ',iy,': z = ',
                   zOrig[(ix*image.data.shape[1]) + iy],', xFit = ',xyFit[(ix*image.data.shape[1]) + iy,0],', yFit = ',xyFit[(ix*image.data.shape[1]) + iy,1])
            nPoints += 1

    print('nPoints set = ',nPoints,': xyOrig.shape = ',xyOrig.shape,', zOrig.shape = ',zOrig.shape,', xyFit.shape = ',xyFit.shape)
    zFit = griddata(xyOrig, zOrig, xyFit, method='nearest')

    for ix in np.arange(0,image.data.shape[0],1):
        for iy in np.arange(0,image.data.shape[1],1):
            image.data[ix, iy] = zFit[(ix*image.data.shape[1]) + iy]
#    fIm = f(horizontalTraces[0][1], verticalTraces[0][1])
#    print('fIm.shape = ',fIm.shape)
    image.write(imFile[:imFile.rfind('.')]+'_interpolated.fits', overwrite=True)

    if False:
    #    coords[0,0,0] = 0.
    #    coords[0,0,1] = 0.
        for i in np.arange(1,image.shape[0],1):
            for j in np.arange(1,image.shape[1],1):
    #            print('i = ',i,', j = ',j,': coords[',i-1,'][',j-1,'][0] = ',coords[i-1][j-1][0])
    #            print('i = ',i,', j = ',j,': verticalTraces[0][1][',i,' = ',verticalTraces[0][1][i])
    #            print('i = ',i,', j = ',j,': verticalTraces[0][1][',i-1,' = ',verticalTraces[0][1][i-1])
    #            print('i = ',i,', j = ',j,': horizontalTraces[0][1][',i,' = ',horizontalTraces[0][1][i])
    #            print('i = ',i,', j = ',j,': horizontalTraces[0][1][',i-1,' = ',horizontalTraces[0][1][i-1])
    #            coords[i,j,0] = coords[i-1][j-1][0] + verticalTraces[0][1][i] - verticalTraces[0][1][i-1]
    #            coords[i,j,1] = coords[i-1][j-1][1] + horizontalTraces[0][1][i] - horizontalTraces[0][1][i-1]
                print('[',i,', ',j,'] = [',coords[i,j,0],', ',coords[i,j,1],']')
    #            STOP
        if verticalTraces[0][1][0] < verticalTraces[0][1][verticalTraces[0][1].shape[0]-1]:
            print('y[0] = ',verticalTraces[0][1][0],' < y[',verticalTraces[0][1].shape[0]-1,'] = ',verticalTraces[0][1][verticalTraces[0][1].shape[0]-1])

        else:
            print('y[0] = ',verticalTraces[0][1][0],' > y[',verticalTraces[0][1].shape[0]-1,'] = ',verticalTraces[0][1][verticalTraces[0][1].shape[0]-1])
    return 1
