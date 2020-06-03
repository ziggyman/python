from astropy.nddata import CCDData
from drUtils import addSuffixToFileName, combine, separateFileList, silentRemove
from drUtils import subtractOverscan, subtractBias, cleanCosmic, flatCorrect,interpolateTraceIm
from drUtils import makeSkyFlat, makeMasterFlat, imDivide, extractSum, calcLineProfile,xCor
from drUtils import findLines, getYAt, calcDispersion, normalizeX, extract,readFileToArr
import matplotlib.pyplot as plt
import numpy as np
import os
from shutil import copyfile
#from myUtils import readFileToArr

overscanSection = '[1983:,:]'
trimSection = '[17:1982,38:97]'
#testPath = '/Users/azuri/spupnik/data/20190501/'
path = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190904/'
allInputFitsList = os.path.join(path,'allFits.list')

suffixes = ['','ot','otz','otzf','otzfi','otzfif','otzx','otzxf','otzxfi','otzxfif']

exptypes = ['BIAS','FLAT','ARC','SCIENCE']
objects = [['*'],['*','Domeflat','Skyflat'],['*'],['*','individual']]
separateFileList(allInputFitsList, suffixes, exptypes, objects, True)

# subtract overscan and trim images
for exptype in exptypes:
    runOnFilesList = os.path.join(path,exptype+'.list')
    outFiles = [addSuffixToFileName(fileName, 'ot') for fileName in runOnFilesList]
    outArrs = subtractOverscan(runOnFilesList,
                               overscanSection,
                               trimSection=trimSection,
                               fitsFilesOut=outFiles,
                               overwrite=True)

# create master bias
otFiles = readFileToArr(os.path.join(path,'BIAS_ot.list'))
masterBias = os.path.join(path,'combinedBias_ot.fits')
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

# subtract master Bias from all non-Biases
for exptype in ['FLAT','ARC','SCIENCE']:
    runOnFilesList = os.path.join(path,exptype+'.list')
    inFiles = [addSuffixToFileName(fileName, 'ot') for fileName in runOnFilesList]
    outFiles = [addSuffixToFileName(fileName, 'otz') for fileName in runOnFilesList]
    outArrs = subtractBias(otFiles,
                           masterBias,
                           fitsFilesOut=outFiles,
                           overwrite=True)

# create master DOMEFLAT
inFiles = readFileToArr(os.path.join(path,'FLATDomeflat_otz.list'))
otzFilesDomeflats = [addSuffixToFileName(fileName, 'otz') for fileName in inFiles]
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

# flatten images using master DOMEFLAT
for exptype in ['FLAT','ARC','SCIENCE']
    otzfFiles = [addSuffixToFileName(fileName, 'otzf') for fileName in inFiles]
    otzfArrs = flatCorrect(otzFiles,
                           masterFlat,
                           norm_value = 1.,
                           fitsFilesOut=otzfFiles)
