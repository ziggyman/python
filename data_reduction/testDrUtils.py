from astropy.nddata import CCDData
from drUtils import addSuffixToFileName, combine, separateFileList, silentRemove
from drUtils import subtractOverscan, subtractBias, cleanCosmic, flatCorrect,interpolateTraceIm
from drUtils import makeSkyFlat, makeMasterFlat, imDivide, extractSum, calcLineProfile,xCor
from drUtils import findLines, getYAt, calcDispersion, normalizeX, extract,readFileToArr
from drUtils import reidentify,getLineProfiles,getBestLineProfile,findGoodLines,rebin
from drUtils import writeFits1D,scombine,continuum
import matplotlib.pyplot as plt
import numpy as np
import os
from shutil import copyfile
#from myUtils import readFileToArr

overscanSection = '[1983:,:]'
trimSection = '[17:1982,38:97]'
#testPath = '/Users/azuri/spupnik/data/20190501/'
testPath = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190904/'
testPath = '/Users/azuri/spectra/saao/saao_sep2019/20190905/'
#testPath = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190501/'
#testPath = '/Volumes/work/azuri/spectra/saao/saao_may2020/20200517/'


def test_separateFileList(inList='/Volumes/work/azuri/spectra/saao/saao_may2019/20190506/allFits.list'):
    suffixes = ['','ot','otz','otzf','otzfi','otzfif','otzx','otzxf','otzxfi','otzxfif']

#    copyfile(inList, inList+'bak')
#    silentRemove(inList[:inList.rfind('/')+1]+'*.list')
#    copyfile(inList+'bak', inList)
    exptypes = ['BIAS','FLAT','ARC','SCIENCE']
    objects = [['*'],['*','Domeflat','Skyflat'],['*'],['*','individual']]
    separateFileList(inList, suffixes, exptypes, objects, True)

def test_combine():
    #combine(ccdImages,
    #        combinerMethod='median',
    #        clippingMethod='None',
    #        clippingParameters=None,
    #        scaling=False,
    #        fitsOutName=None)
    inFiles = readFileToArr(os.path.join(testPath,'BIAS.list'))
    print('test_combine: inFiles = ',inFiles)
    combinedImage = combine(inFiles, fitsOutName=os.path.join(testPath,'combinedBias.fits'))
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
                    scaling=False,
                    fitsOutName=os.path.join(testPath,'combinedBias_av.fits'))
    print('average minmax: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='extrema',
                    clippingParameters={'nlow':2,
                                        'nhigh':2},
                    scaling=False,
                    fitsOutName=os.path.join(testPath,'combinedBias_av.fits'))
    print('average extrema: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':0,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=False,
                    fitsOutName=os.path.join(testPath,'combinedBias_av.fits'))
    print('average sigma 0: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':1,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=False,
                    fitsOutName=os.path.join(testPath,'combinedBias_av.fits'))
    print('average sigma 1: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':3,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=False,
                    fitsOutName=os.path.join(testPath,'combinedBias_av.fits'))
    print('average sigma 3: mean(combinedImage) = ',np.mean(combinedImage))

#subtractOverscan(fitsFilesIn, overscanSection, trimSection=None, fitsFilesOut=None, overwrite=True):
def test_subtractOverscan():
    inFiles = readFileToArr(os.path.join(testPath,'BIAS.list'))
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
    inFiles = readFileToArr(os.path.join(testPath,'BIAS.list'))
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
    masterBias = os.path.join(testPath,'combinedBias_ot.fits')
    combinedImage = combine(otFiles,
                            combinerMethod='median',
                            clippingMethod='sigma',
                            clippingParameters={'niter':0,
                                                'low_thresh':-3.,
                                                'high_thresh':3.,
                                                'func':np.ma.median},
                            scaling=False,
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
    inFiles = readFileToArr(os.path.join(testPath,'BIAS.list'))
    if inFiles[0].rfind('/') == -1:
        inFiles = [os.path.join(path, fileName) for fileName in inFiles]
    otFiles = [addSuffixToFileName(fileName, 'ot') for fileName in inFiles]
    outArrs = subtractOverscan(inFiles,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFiles,
                               overwrite=True)

    # create master bias
    masterBias = os.path.join(testPath,'combinedBias_ot.fits')
    combinedImage = combine(otFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':0,
                                        'low_thresh':-3.,
                                        'high_thresh':3.,
                                        'func':np.ma.median},
                    scaling=False,
                    fitsOutName=masterBias)
    print('average sigma 0: mean(combinedImage) = ',np.mean(combinedImage))

    inFiles = readFileToArr(os.path.join(testPath,'SCIENCE.list'))
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

def test_flatCorrect(objectFiles = os.path.join(testPath,'SCIENCE.list')):
    path = objectFiles[0:objectFiles.rfind('/')+1]
    inFilesBias = readFileToArr(os.path.join(path,'BIAS.list'))
    if inFilesBias[0].rfind('/') == -1:
        inFilesBias = [os.path.join(path, fileName) for fileName in inFilesBias]
    otFilesBias = [addSuffixToFileName(fileName, 'ot') for fileName in inFilesBias]
    otArrsBias = subtractOverscan(inFilesBias,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFilesBias,
                               overwrite=True)

    # create master bias
    masterBias = os.path.join(path,'combinedBias_ot.fits')
    combinedBias = combine(otFilesBias,
                            combinerMethod='median',
                            clippingMethod='sigma',
                            clippingParameters={'niter':0,
                                                'low_thresh':-3.,
                                                'high_thresh':3.,
                                                'func':np.ma.median},
                            scaling=False,
                            fitsOutName=masterBias)
    print('average sigma 0: mean(combinedBias) = ',np.mean(combinedBias))

    inFilesDomeflats = readFileToArr(os.path.join(path,'FLATDomeflat.list'))
    otFilesDomeflats = [addSuffixToFileName(fileName, 'ot') for fileName in inFilesDomeflats]
    otArrsDomeflats = subtractOverscan(inFilesDomeflats,
                                       overscanSection,
                                       trimSection=trimSection,
                                       fitsFilesOut=otFilesDomeflats,
                                       overwrite=True)

    otzFilesDomeflats = [addSuffixToFileName(fileName, 'otz') for fileName in inFilesDomeflats]
    otzArrsDomeflats = subtractBias(otFilesDomeflats,
                                    masterBias,
                                    fitsFilesOut=otzFilesDomeflats,
                                    overwrite=True)

    combinedFlat = os.path.join(path,'combinedFlat.fits')
    print('creating combinedFlat <'+combinedFlat+'>')
    print('otzFilesDomeflats = ',otzFilesDomeflats)
    flat = combine(otzFilesDomeflats,
                   combinerMethod='median',
                   clippingMethod='sigma',
                   clippingParameters={'niter':2,
                                       'low_thresh':-3.,
                                       'high_thresh':3.,
                                       'func':np.ma.median},
                   scaling=False,
                   minVal=0.0001,
                   fitsOutName=combinedFlat)

    masterFlat = os.path.join(path, 'masterDomeFlat.fits')
    smoothedFlat = os.path.join(path, 'smoothedDomeFlat.fits')
    makeMasterFlat(combinedFlat, 9, 80., outFileNameMasterFlat=masterFlat, outFileNameMasterFlatSmoothed=smoothedFlat)

    flatOut = addSuffixToFileName(combinedFlat, 'flattened')
    flat = flatCorrect([combinedFlat],
                       masterFlat,
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
    otzArrs = subtractBias(otFiles,
                           masterBias,
                           fitsFilesOut=otzFiles,
                           overwrite=True)

#    otzxFiles = [addSuffixToFileName(fileName, 'otzx') for fileName in inFiles]
#    otzxArrs = cleanCosmic(otzFiles,
#                           fitsFilesOut=otzxFiles,
#                           overwrite=True)

    otzfFiles = [addSuffixToFileName(fileName, 'otzf') for fileName in inFiles]
    otzfArrs = flatCorrect(otzFiles,
                           masterFlat,
                           norm_value = 1.,
                           fitsFilesOut=otzfFiles)

    return otzfFiles

def test_combineSkyFlats(inFiles = os.path.join(testPath,'FLATSkyflat.list')):
    inFiles = test_flatCorrect(inFiles)

    combinedFlat = inFiles[0][:inFiles[0].rfind('/')+1]+'combinedSkyFlat.fits'
    flat = combine(inFiles,
                   combinerMethod='median',
                   clippingMethod='sigma',
                   clippingParameters={'niter':2,
                                       'low_thresh':-3.,
                                       'high_thresh':3.,
                                       'func':np.ma.median},
                   scaling=True,
                   minVal=0.0001,
                   fitsOutName=combinedFlat)
    return combinedFlat

def testInterpolateTraceIm():
    interpolateTraceIm([os.path.join(testPath,'combinedSkyFlat.fits')],#FLAT_Domeflat_a1061047_otz.fits
                       os.path.join(testPath,'database/apvertical_trace'),
                       os.path.join(testPath,'database/aphorizontal_tracer90flipl'))

def testMakeSkyFlat():
    makeSkyFlat(os.path.join(testPath,'combinedSkyFlati.fits'),
                os.path.join(testPath,'combinedSkyFlati_flattened.fits'),
                7)

def test_skyFlatCorrect(objectFiles = os.path.join(testPath,'SCIENCE.list')):
    path = objectFiles[0:objectFiles.rfind('/')+1]
    inFilesBias = readFileToArr(os.path.join(path,'BIAS.list'))
    if inFilesBias[0].rfind('/') == -1:
        inFilesBias = [os.path.join(path, fileName) for fileName in inFilesBias]
    otFilesBias = [addSuffixToFileName(fileName, 'ot') for fileName in inFilesBias]
    otArrsBias = subtractOverscan(inFilesBias,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=otFilesBias,
                               overwrite=True)

    # create master bias
    masterBias = os.path.join(path,'combinedBias_ot.fits')
    combinedBias = combine(otFilesBias,
                            combinerMethod='median',
                            clippingMethod='sigma',
                            clippingParameters={'niter':0,
                                                'low_thresh':-3.,
                                                'high_thresh':3.,
                                                'func':np.ma.median},
                            scaling=False,
                            fitsOutName=masterBias)
    print('average sigma 0: mean(combinedBias) = ',np.mean(combinedBias))

    inFilesDomeflats = readFileToArr(os.path.join(path,'FLATDomeflat.list'))
    otFilesDomeflats = [addSuffixToFileName(fileName, 'ot') for fileName in inFilesDomeflats]
    otArrsDomeflats = subtractOverscan(inFilesDomeflats,
                                       overscanSection,
                                       trimSection=trimSection,
                                       fitsFilesOut=otFilesDomeflats,
                                       overwrite=True)

    otzFilesDomeflats = [addSuffixToFileName(fileName, 'otz') for fileName in inFilesDomeflats]
    otzArrsDomeflats = subtractBias(otFilesDomeflats,
                                    masterBias,
                                    fitsFilesOut=otzFilesDomeflats,
                                    overwrite=True)

    combinedFlat = os.path.join(path,'combinedFlat.fits')
    flat = combine(otzFilesDomeflats,
                   combinerMethod='median',
                   clippingMethod='sigma',
                   clippingParameters={'niter':2,
                                       'low_thresh':-3.,
                                       'high_thresh':3.,
                                       'func':np.ma.median},
                   scaling=False,
                   minVal=0.0001,
                   fitsOutName=combinedFlat)

    masterFlat = os.path.join(path, 'masterDomeFlat.fits')
    smoothedFlat = os.path.join(path, 'smoothedDomeFlat.fits')
    makeMasterFlat(combinedFlat, 9, 80., outFileNameMasterFlat=masterFlat, outFileNameMasterFlatSmoothed=smoothedFlat)

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
    otzArrs = subtractBias(otFiles,
                           masterBias,
                           fitsFilesOut=otzFiles,
                           overwrite=True)

#    otzxFiles = [addSuffixToFileName(fileName, 'otzx') for fileName in inFiles]
#    otzxArrs = cleanCosmic(otzFiles,
#                           fitsFilesOut=otzxFiles,
#                           overwrite=True)

    otzfFiles = [addSuffixToFileName(fileName, 'otzf') for fileName in inFiles]
    otzfArrs = flatCorrect(otzFiles,
                           masterFlat,
                           fitsFilesOut=otzfFiles)

#    otzffFiles = [addSuffixToFileName(fileName, 'otzff') for fileName in inFiles]
#    otzffArrs = imDivide(otzFiles, smoothedFlat, otzffFiles)


    test_combineSkyFlats()
    otzfArrs = flatCorrect([os.path.join(testPath,'combinedSkyFlat.fits')],
                           masterFlat,
                           fitsFilesOut=[os.path.join(testPath,'combinedSkyFlat_f.fits')])
#    otzffArrs = imDivide([os.path.join(testPath,'combinedSkyFlat_f.fits')], smoothedFlat, [os.path.join(testPath,'combinedSkyFlat_ff.fits')])

    interpolateTraceIm([os.path.join(testPath,'combinedSkyFlat_f.fits')],#FLAT_Domeflat_a1061047_otz.fits
                       os.path.join(testPath,'database/apvertical_trace'),
                       os.path.join(testPath,'database/aphorizontal_tracer90flipl'))

    interpolateTraceIm(otzfFiles,
                       os.path.join(testPath,'database/apvertical_trace'),
                       os.path.join(testPath,'database/aphorizontal_tracer90flipl'))

    makeSkyFlat(os.path.join(testPath,'combinedSkyFlat_fi.fits'),
                os.path.join(testPath,'combinedSkyFlat_fi_flattened.fits'),
                7)

    otzfiFiles = [addSuffixToFileName(fileName, 'otzfi') for fileName in inFiles]
    otzfifFiles = [addSuffixToFileName(fileName, 'otzfif') for fileName in inFiles]
    otzfifArrs = flatCorrect(otzfiFiles,
                              os.path.join(testPath,'combinedSkyFlat_fi_flattened.fits'),
                              fitsFilesOut=otzfifFiles)


    return otzfifFiles

def testExtractSum(oneDImageIn):
    extractSum(oneDImageIn)

def main():
    if False:
        test_separateFileList(os.path.join(testPath,'allFits.list'))
    if False:
        test_combine()
        test_subtractOverscan()
        test_subtractBias()
    if False:
        test_cleanCosmic()
    if False:
        test_flatCorrect(os.path.join(testPath,'allFits.list'))
        test_combineSkyFlats()
        testInterpolateTraceIm()
        testMakeSkyFlat()
        test_skyFlatCorrect()
    if False:
        test_skyFlatCorrect(os.path.join(testPath,'ARC.list'))
        inFiles = readFileToArr(os.path.join(testPath,'SCIENCEWKK98_201_otzfif.list'))
        wkk = combine(inFiles,
                       combinerMethod='average',
                       clippingMethod='sigma',
                       clippingParameters={'niter':2,
                                           'low_thresh':-3.,
                                           'high_thresh':3.,
                                           'func':np.ma.median},
                       scaling=True,
                       minVal=0.0001,
                       fitsOutName=os.path.join(testPath,'SCIENCEWKK98_201_otzxfif_combined.fits'))

    if False:
        cosmicParameters = {'sigclip':4.,
                            'cleantype':'meanmask',
                            'gain':0.219,
                            'readnoise':3.4,
                            'sepmed':False,
                            'fsmode':'convolve',
                            'psfmodel':'moffat',
                            'verbose':True,
                            'niter':1,
                            'psfbeta':1.5}
        inFiles = readFileToArr('/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/do_cosmics.list')
        outFiles = [addSuffixToFileName(fileName, 'x') for fileName in inFiles]
        outArrs = cleanCosmic(inFiles,
                              cosmicParameters = cosmicParameters,
                              fitsFilesOut=outFiles,
                              overwrite=True)

    if True:
        arc = 'ARC_MPA_J1354-6337_a1171127_otzf.fits'
        arcInterp = arc[:-5]+'i.fits'
    if False:
        lineProfiles = getLineProfiles(os.path.join(testPath,arc))
        bestLineProfile = getBestLineProfile(lineProfiles)

    if True:
        oneDSpec = extractSum(os.path.join(testPath,arc),'row')
        oneDSpecInterp = extractSum(os.path.join(testPath,arcInterp),'row')
        plt.plot(oneDSpec,label='raw')
        plt.plot(oneDSpecInterp,label='interpolated')
        plt.legend()
        plt.show()

        xSpec = np.arange(0,oneDSpecInterp.shape[0],1.)
        print('xSpec = ',xSpec)
    if False:
    #    xCor([xSpec,oneDSpecInterp],bestLineProfile)
        findGoodLines(xSpec,
                      oneDSpecInterp,
                      bestLineProfile,
                      outFileNameAllLines='/Volumes/work/azuri/spectra/saao/saao_refSpec_lines.dat',
                      outFileNameGoodLines='/Volumes/work/azuri/spectra/saao/saao_refSpec_lines_good.dat')
    if True:
        coeffs,rms = calcDispersion('/Volumes/work/azuri/spectra/saao/saao_refspec_lines_identified_good.dat', xRange=[0,xSpec[xSpec.shape[0]-1]], degree=3)
        xSpecNorm = normalizeX(xSpec)
        wLenSpec = np.polynomial.legendre.legval(xSpecNorm, coeffs)
        plt.plot(wLenSpec,oneDSpecInterp)
        plt.show()
    if False:
        extractList = [['SCIENCE_CTI173557.97-271213.963_a1061524_otzfif.fits', [20,29], [30,35], [16,18]]]
        extractList.append(['SCIENCE_Pre_22_a1061516_otzfif.fits', [25,32], [18,24],[33,41]])
        twoDSpec = os.path.join(testPath,extractList[1][0])
        specOut = twoDSpec[:-5]+'Ec.dat'

        extract(twoDSpec, specOut, extractList[1][1], extractList[1][2], extractList[1][3], 'row')

def test_scombine():
    inputList = '/Users/azuri/spectra/saao/saao_sep2019/20190905/MPA_J1924+0038_EcdF.list'
    outputFileName = '/Users/azuri/spectra/saao/saao_sep2019/20190905/MPA_J1924+0038_EcdF_combinedTest.fits'
    method = 'median'
    lowReject = 2.
    highReject = 2.
    adjustSigLevels = False
    useMean = False
    scombine(inputList,
             outputFileName[:-5]+'_'+method+'_lowRej'+str(lowReject)+'_highRej'+str(highReject)+'_adj='+str(adjustSigLevels),
             method=method,
             lowReject=lowReject,
             highReject=highReject,
             adjustSigLevels=adjustSigLevels,
             useMean=useMean,
             display=True)
    method = 'mean'
    scombine(inputList,
             outputFileName[:-5]+'_'+method+'_lowRej'+str(lowReject)+'_highRej'+str(highReject)+'_adj='+str(adjustSigLevels),
             method=method,
             lowReject=lowReject,
             highReject=highReject,
             adjustSigLevels=adjustSigLevels,
             useMean=useMean,
             display=True)

def test_continuum():
    spectrumFileNameIn = '/Users/azuri/spectra/saao/saao_sep2019/20190906/SCIENCE_NGC7293-NE1_a1171207_otzxfifEcdF.fits'
    spectrumFileNameOut = '/Users/azuri/spectra/saao/saao_sep2019/20190906/SCIENCE_NGC7293-NE1_a1171207_otzxfifEcdF_continuumTest.fits'
    fittingFunction = np.polynomial.legendre.legfit
    evalFunction = np.polynomial.legendre.legval
    order = 9
    nIterReject = 3
    nIterFit = 3
    lowReject = 2.
    highReject = 2.
    useMean = False
    continuum(spectrumFileNameIn, spectrumFileNameOut, fittingFunction, evalFunction, order, nIterReject, nIterFit, lowReject, highReject, type='difference', adjustSigLevels=False, useMean=useMean, display=True)

if __name__ == "__main__":
    #main()
    #test_cleanCosmic()
    test_scombine()
    test_continuum()
    if False:
#        arc = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190904/ARC_Feige110_a1171059_otzf.fits'
        arc = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190905/ARC_MPA_J1354-6337_a1171127_otzf.fits'
        arcInterp = arc[:-5]+'i.fits'
#        lineProfiles = getLineProfiles(arc)
#        bestLineProfile = getBestLineProfile(lineProfiles)

        #arc =
        oneDSpec = extractSum(arc,'row')
        oneDSpecInterp = extractSum(arcInterp,'row')
        plt.plot(oneDSpec,label='raw')
        plt.plot(oneDSpecInterp,label='interpolated')
        plt.legend()
        plt.show()

#        xSpec = np.arange(0,oneDSpecInterp.shape[0],1.)
#        print('xSpec = ',xSpec)
    #    xCor([xSpec,oneDSpecInterp],bestLineProfile)
#        findGoodLines(xSpec,
#                      oneDSpecInterp,
#                      bestLineProfile,
#                      outFileNameAllLines=arc[:-5]+'_lines.dat',
#                      outFileNameGoodLines=arc[:-5]+'_lines_good.dat')

#        arcSpec = extractSum(arcInterp,'row')
#        print('arcSpec = ',arcSpec.shape,': ',arcSpec)
#        profile = bestLineProfile
#        print('profile = ',profile)
        refApDef = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190905/database/aptmpARC_MPA_J1354-6337_a1171127_otzf'
        lineList = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_gr7_angle16_3_lines_identified_good.dat'
        lineListNew, coeffs, xRange, rms = reidentify(arcInterp,
                                                      arc,
                                                      refApDef,
                                                      lineList,
                                                      lineListOut=arc[:arc.rfind('.')]+'_lines.dat',
                                                      specOut=arc[:-5]+'Ecd.fits')
        xSpec = np.arange(xRange[0],xRange[1]+1,1.)
        xSpecNorm = normalizeX(xSpec)
        wLenSpec = np.polynomial.legendre.legval(xSpecNorm, coeffs)

        dLam = np.min([np.absolute(wLenSpec[1]-wLenSpec[0]),np.absolute(wLenSpec[wLenSpec.shape[0]-1]-wLenSpec[wLenSpec.shape[0]-2])])
        resampled = np.arange(np.min([wLenSpec[0], wLenSpec[wLenSpec.shape[0]-1]]), np.max([wLenSpec[0], wLenSpec[wLenSpec.shape[0]-1]]), dLam)
        resampledSpec = rebin(wLenSpec, oneDSpecInterp, resampled, preserveFlux = False)
        print('wLenSpec = ',wLenSpec.shape,': ',wLenSpec)
        print('resampled wavelength = ',resampled.shape,': ',resampled)
        plt.plot(wLenSpec,oneDSpecInterp,label='original')
        plt.plot(resampled,resampledSpec,label='resampled')
        plt.legend()
        plt.show()


        writeFits1D(resampledSpec, arcInterp[:-5]+'Ec.fits', wavelength=None, header=arcInterp, CRVAL1=resampled[0], CRPIX1=1, CDELT1=dLam)
