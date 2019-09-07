from astropy.nddata import CCDData
from drUtils import addSuffixToFileName, combine, separateFileList, silentRemove
from drUtils import subtractOverscan, subtractBias, cleanCosmic, flatCorrect,interpolateTraceIm
from drUtils import makeSkyFlat, makeMasterFlat, imDivide
import numpy as np
import os
from shutil import copyfile

overscanSection = '[1983:,:]'
trimSection = '[17:1982,38:97]'
workPath = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190506/'


def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

def getListOfFiles(fname):
    fList = readFileToArr(fname)
    if fList[0].rfind('/') == -1:
        fList = [os.path.join(workPath, fileName) for fileName in fList]
    return fList

inList=os.path.join(workPath,'allFits.list')
suffixes = ['','ot','otz','otzf','otzfi','otzfif','otzx','otzxf','otzxfi','otzxfif']

#    copyfile(inList, inList+'bak')
#    silentRemove(inList[:inList.rfind('/')+1]+'*.list')
#    copyfile(inList+'bak', inList)
exptypes = ['BIAS','FLAT','ARC','SCIENCE']
objects = [['*'],['*','Domeflat','Skyflat'],['*'],['*','individual']]
separateFileList(inList, suffixes, exptypes, objects, True)

objectFiles = os.path.join(workPath,'SCIENCE.list')

# subtract overscan and trim all images
for inputList in ['ARC', 'BIAS', 'FLAT', 'SCIENCE']:
    subtractOverscan(getListOfFiles(os.path.join(workPath,inputList+'.list')),
                     overscanSection,
                     trimSection=trimSection,
                     fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_ot.list')),
                     overwrite=True)

# create master bias
masterBias = os.path.join(workPath,'combinedBias_ot.fits')
combinedBias = combine(getListOfFiles(os.path.join(workPath,'BIAS_ot.list')),
                       combinerMethod='median',
                       clippingMethod='sigma',
                       clippingParameters={'niter':0,
                                           'low_thresh':-3.,
                                           'high_thresh':3.,
                                           'func':np.ma.median},
                       scaling=False,
                       fitsOutName=masterBias)
print('average sigma 0: mean(combinedBias) = ',np.mean(combinedBias))

# subtract masterBias from all images
for inputList in ['ARC', 'FLAT', 'SCIENCE']:
    subtractBias(getListOfFiles(os.path.join(workPath,inputList+'_ot.list')),
                 masterBias,
                 fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                 overwrite=True)

# create master DomeFlat
combinedFlat = os.path.join(workPath,'combinedFlat.fits')
print('creating combinedFlat <'+combinedFlat+'>')
flat = combine(getListOfFiles(os.path.join(workPath,'FLATDomeflat_otz.list')),
               combinerMethod='median',
               clippingMethod='sigma',
               clippingParameters={'niter':2,
                                   'low_thresh':-3.,
                                   'high_thresh':3.,
                                   'func':np.ma.median},
               scaling=False,
               minVal=0.0001,
               fitsOutName=combinedFlat)

masterFlat = os.path.join(workPath, 'masterDomeFlat.fits')
smoothedFlat = os.path.join(workPath, 'smoothedDomeFlat.fits')
makeMasterFlat(combinedFlat,
               9,
               80.,
               outFileNameMasterFlat=masterFlat,
               outFileNameMasterFlatSmoothed=smoothedFlat)

# apply master DomeFlat to ARCs, SkyFlats, and SCIENCE frames
for inputList in ['ARC','SCIENCE','FLATSkyflat']:
    flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                masterFlat,
                norm_value = 1.,
                fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzf.list')))

# interpolate images to get straight dispersion and spectral features
for inputList in ['ARC', 'SCIENCE']:
    interpolateTraceIm(getListOfFiles(os.path.join(workPath,inputList+'_otzf.list')),
                       os.path.join(workPath,'database/apvertical_trace'),
                       os.path.join(workPath,'database/aphorizontal_tracer90flipl'))

# create master SkyFlat
combinedSkyFlat = os.path.join(workPath,'combinedSkyFlat.fits')
print('creating combinedSkyFlat <'+combinedSkyFlat+'>')
combine(getListOfFiles(os.path.join(workPath,'FLATSkyflat_otzf.list')),
        combinerMethod='median',
        clippingMethod='sigma',
        clippingParameters={'niter':2,
                            'low_thresh':-3.,
                            'high_thresh':3.,
                            'func':np.ma.median},
        scaling=True,
        minVal=0.0001,
        fitsOutName=combinedSkyFlat)

interpolateTraceIm([combinedSkyFlat],
                    os.path.join(workPath,'database/apvertical_trace'),
                    os.path.join(workPath,'database/aphorizontal_tracer90flipl'))

makeSkyFlat(os.path.join(workPath,'combinedSkyFlati.fits'),
            os.path.join(workPath,'combinedSkyFlati_flattened.fits'),
            7)

flatCorrect(getListOfFiles(os.path.join(workPath,'SCIENCE_otzfi.list')),
            os.path.join(workPath,'combinedSkyFlati_flattened.fits'),
            fitsFilesOut=getListOfFiles(os.path.join(workPath,'SCIENCE_otzfif.list')))
