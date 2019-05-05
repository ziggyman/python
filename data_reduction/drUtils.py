#from astropy import units as u
#from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
from astropy import units as u
from astropy.nddata import CCDData
from ccdproc import Combiner
import numpy as np
#import hammer
#import numpy as np
import os
from shutil import copyfile
#from sklearn import linear_model

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
# <scaling>: boolean - normalize each image by its mean before combining
# <fitsOutName>: string, default is don't write (None)
# <overwrite>:       boolean, default is True
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
    combiner = Combiner(ccdData)
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
        ccdDataNoOverscan = subtract_overscan(ccdData, fits_section=overscanSection)
        if trimSection is None:
            if fitsFilesOut is not None:
                ccdDataNoOverscan.write(fitsFilesOut[iFile], overwrite=overwrite)
            dataOut.append(trimmed)
        else:
            trimmed = trim_image(ccdDataNoOverscan, fits_section=trimSection)
            if fitsFilesOut is not None:
                trimmed.write(fitsFilesOut[iFile], overwrite=overwrite)
            dataOut.append(trimmed)

