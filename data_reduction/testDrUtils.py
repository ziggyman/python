from astropy.nddata import CCDData
from drUtils import addSuffixToFileName, combine, separateFileList, silentRemove
from drUtils import subtractOverscan, subtractBias, cleanCosmic, flatCorrect,interpolateTraceIm
import numpy as np
import os
from shutil import copyfile
from myUtils import readFileToArr

overscanSection = '[1983:,:]'
trimSection = '[17:1982,38:97]'

def test_separateFileList(inList='/Volumes/work/azuri/spectra/saao/saao_may2019/20190506/allFits.list'):
    suffixes = ['','ot','otz','otzx','otzxf']

#    copyfile(inList, inList+'bak')
#    silentRemove(inList[:inList.rfind('/')+1]+'*.list')
#    copyfile(inList+'bak', inList)
    separateFileList(inList, suffixes, True)

def test_combine():
    #combine(ccdImages,
    #        combinerMethod='median',
    #        clippingMethod='None',
    #        clippingParameters=None,
    #        scaling=False,
    #        fitsOutName=None)
    inFiles = readFileToArr('/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/bias.list')
    combinedImage = combine(inFiles, fitsOutName='/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/combinedBias.fits')
    print('combinedImage = ',type(combinedImage))
    print('dir(combinedImage) = ',dir(combinedImage))
    print('combinedImage.shape = ',combinedImage.shape)
    print('combinedImage.ndim = ',combinedImage.ndim)
    print('combinedImage.size = ',combinedImage.size)
    print('mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='minmax',
                    clippingParameters={'min_clip':0.,
                                        'max_clip':700.},
                    scaling=True,
                    fitsOutName='/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average minmax: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='extrema',
                    clippingParameters={'nlow':2,
                                        'nhigh':2},
                    scaling=True,
                    fitsOutName='/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average extrema: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':0,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=True,
                    fitsOutName='/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average sigma 0: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':1,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=True,
                    fitsOutName='/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average sigma 1: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':3,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=True,
                    fitsOutName='/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average sigma 3: mean(combinedImage) = ',np.mean(combinedImage))

#subtractOverscan(fitsFilesIn, overscanSection, trimSection=None, fitsFilesOut=None, overwrite=True):
def test_subtractOverscan():
    inFiles = readFileToArr('/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/bias.list')
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=None,
                               fitsFilesOut=None,
                               overwrite=True)
    print('len(inFiles) = ',len(inFiles),', len(outFiles) = ',len(outArrs))
    meanDiffs = []
    for iFile in np.arange(0, len(inFiles), 1):
        meanA = np.mean(CCDData.read(inFiles[iFile], unit="adu"))
        meanB = np.mean(outArrs[iFile])
        meanDiff = meanA-meanB
        meanDiffs.append(meanDiff)
        print('iFile = ',iFile,': meanA - meanB = ',meanDiff)

    outFiles = [addSuffixToFileName(fileName, 'o') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=None,
                               fitsFilesOut=outFiles,
                               overwrite=True)
    for iFile in np.arange(0, len(inFiles), 1):
        meanA = np.mean(CCDData.read(inFiles[iFile], unit="adu"))
        meanB = np.mean(CCDData.read(outFiles[iFile], unit="adu"))
        meanDiff = meanA-meanB
        print('iFile = ',iFile,': meanA - meanB = ',meanDiff,
              ': difference to previous calculation = ',meanDiff - meanDiffs[iFile])

    outFiles = [addSuffixToFileName(fileName, 'ot') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=outFiles,
                               overwrite=True)
    for iFile in np.arange(0, len(inFiles), 1):
        meanA = np.mean(CCDData.read(inFiles[iFile], unit="adu"))
        meanB = np.mean(CCDData.read(outFiles[iFile], unit="adu"))
        meanDiff = meanA-meanB
        print('iFile = ',iFile,': meanA - meanB = ',meanDiff,
              ': difference to previous calculation = ',meanDiff - meanDiffs[iFile])

    try:
        outArrs = subtractOverscan(inFiles,
                                   overscanSection,
                                   trimSection=trimSection,
                                   fitsFilesOut=outFiles,
                                   overwrite=False)
    except:
        pass

def test_subtractBias():
    inFiles = readFileToArr('/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/bias.list')
    otFiles = [addSuffixToFileName(fileName, 'ot') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFiles,
                               overwrite=True)
    print('len(inFiles) = ',len(inFiles),', len(otFiles) = ',len(outArrs))
    meanDiffs = []
    for iFile in np.arange(0, len(inFiles), 1):
        meanA = np.mean(CCDData.read(inFiles[iFile], unit="adu"))
        meanB = np.mean(outArrs[iFile])
        meanDiff = meanA-meanB
        meanDiffs.append(meanDiff)
        print('iFile = ',iFile,': meanA - meanB = ',meanDiff)

    # create master bias
    masterBias = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/combinedBias_ot.fits'
    combinedImage = combine(otFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':0,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=True,
                    fitsOutName=masterBias)
    print('average sigma 0: mean(combinedImage) = ',np.mean(combinedImage))

    otzFiles = [addSuffixToFileName(fileName, 'otz') for fileName in inFiles]
    outArrs = subtractBias(otFiles,
                           masterBias,
                           fitsFilesOut=otzFiles,
                           overwrite=True)
    for iFile in np.arange(0, len(inFiles), 1):
        meanA = np.mean(CCDData.read(otFiles[iFile], unit="adu"))
        meanB = np.mean(CCDData.read(otzFiles[iFile], unit="adu"))
        meanDiff = meanA-meanB
        print('iFile = ',iFile,': meanA - meanB = ',meanDiff,
              ': difference to previous calculation = ',meanDiff - meanDiffs[iFile])

def test_cleanCosmic():
    inFiles = readFileToArr('/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/bias.list')
    if inFiles[0].rfind('/') == -1:
        inFiles = [os.path.join(path, fileName) for fileName in inFiles]
    otFiles = [addSuffixToFileName(fileName, 'ot') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFiles,
                               overwrite=True)

    # create master bias
    masterBias = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/combinedBias_ot.fits'
    combinedImage = combine(otFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':0,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=True,
                    fitsOutName=masterBias)
    print('average sigma 0: mean(combinedImage) = ',np.mean(combinedImage))

    inFiles = readFileToArr('/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/hd82106.list')
    if inFiles[0].rfind('/') == -1:
        inFiles = [os.path.join(path, fileName) for fileName in inFiles]
    otFiles = [addSuffixToFileName(fileName, 'ot') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFiles,
                               overwrite=True)

    otzFiles = [addSuffixToFileName(fileName, 'otz') for fileName in inFiles]
    outArrs = subtractBias(otFiles,
                           masterBias,
                           fitsFilesOut=otzFiles,
                           overwrite=True)

    otzxFiles = [addSuffixToFileName(fileName, 'otzx_std') for fileName in inFiles]
    outArrs = cleanCosmic(inFiles,
                          fitsFilesOut=otzxFiles,
                          overwrite=True)

    cosmicParameters = {'sigclip':4., 'cleantype':'medmask'}
    otzxFiles = [addSuffixToFileName(fileName, 'otzx_lacos_pars') for fileName in inFiles]
    outArrs = cleanCosmic(inFiles,
                          cosmicParameters = cosmicParameters,
                          fitsFilesOut=otzxFiles,
                          overwrite=True)

    cosmicParameters = {'sigclip':4., 'cleantype':'medmask', 'testkey':0}
    try:
        outArrs = cleanCosmic(inFiles,
                              cosmicParameters = cosmicParameters,
                              fitsFilesOut=otzxFiles,
                              overwrite=True)
    except:
        pass

    otzxFiles = [addSuffixToFileName(fileName, 'otzx_median') for fileName in inFiles]
    outArrs = cleanCosmic(inFiles,
                          cosmicMethod='median',
                          fitsFilesOut=otzxFiles,
                          overwrite=True)

def test_flatCorrect(objectFiles = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/hd82106.list'):
    path = objectFiles[0:objectFiles.rfind('/')+1]
    inFiles = readFileToArr(os.path.join(path,'bias.list'))
    if inFiles[0].rfind('/') == -1:
        inFiles = [os.path.join(path, fileName) for fileName in inFiles]
    otFiles = [addSuffixToFileName(fileName, 'ot') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFiles,
                               overwrite=True)

    # create master bias
    masterBias = os.path.join(path,'combinedBias_ot.fits')
    combinedImage = combine(otFiles,
                            combinerMethod='average',
                            clippingMethod='sigma',
                            clippingParameters={'niter':0,
                                                'low_thresh':-3.,
                                                'high_thresh':3.,
                                                'func':np.ma.median},
                            scaling=True,
                            fitsOutName=masterBias)
    print('average sigma 0: mean(combinedImage) = ',np.mean(combinedImage))

    inFiles = readFileToArr(os.path.join(path,'flat.list'))
    otFiles = [addSuffixToFileName(fileName, 'ot') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFiles,
                               overwrite=True)

    otzFiles = [addSuffixToFileName(fileName, 'otz') for fileName in inFiles]
    outArrs = subtractBias(otFiles,
                           masterBias,
                           fitsFilesOut=otzFiles,
                           overwrite=True)

    combinedFlat = os.path.join(path,'combinedFlat.fits')
    flat = combine(otzFiles,
                   combinerMethod='average',
                   clippingMethod='sigma',
                   clippingParameters={'niter':0,
                                       'low_thresh':-3.,
                                       'high_thresh':3.,
                                       'func':np.ma.median},
                   scaling=True,
                   fitsOutName=combinedFlat)

    flatOut = addSuffixToFileName(combinedFlat, 'flattened')
    flat = flatCorrect([combinedFlat],
                       combinedFlat,
                       fitsFilesOut=[flatOut])


    inFiles = readFileToArr(objectFiles)
    if inFiles[0].rfind('/') == -1:
        inFiles = [os.path.join(path, fileName) for fileName in inFiles]
    otFiles = [addSuffixToFileName(fileName, 'ot') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFiles,
                               overwrite=True)

    otzFiles = [addSuffixToFileName(fileName, 'otz') for fileName in inFiles]
    outArrs = subtractBias(otFiles,
                           masterBias,
                           fitsFilesOut=otzFiles,
                           overwrite=True)

    otzxFiles = [addSuffixToFileName(fileName, 'otzx') for fileName in inFiles]
    outArrs = cleanCosmic(otzFiles,
                          fitsFilesOut=otzxFiles,
                          overwrite=True)

    otzxfFiles = [addSuffixToFileName(fileName, 'otzxf') for fileName in inFiles]
    outArrs = flatCorrect(otzxFiles,
                          combinedFlat,
                          fitsFilesOut=otzxfFiles)

    return otzxfFiles

def test_combineSkyFlats(inFiles = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/skyflats.list'):
    inFiles = test_flatCorrect(inFiles)

    combinedFlat = inFiles[0][:inFiles[0].rfind('/')+1]+'combinedSkyFlat.fits'
    flat = combine(inFiles,
                   combinerMethod='average',
                   clippingMethod='sigma',
                   clippingParameters={'niter':0,
                                       'low_thresh':-3.,
                                       'high_thresh':3.,
                                       'func':np.ma.median},
                   scaling=True,
                   fitsOutName=combinedFlat)
    return combinedFlat

path = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/'
#test_separateFileList(os.path.join(path,'allFits.list'))
#test_combine()
#test_subtractOverscan()
#test_subtractBias()
#test_cleanCosmic()
#test_flatCorrect(os.path.join(path,'allFits.list'))
#test_combineSkyFlats()
interpolateTraceIm('/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/vertical_trace.fits',#FLAT_Domeflat_a1061047_otz.fits',
            '/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/database/apvertical_trace',
            '/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/database/aphorizontal_tracer90')
