#from astropy.nddata import CCDData
from astropy.coordinates import EarthLocation
from drUtils import addSuffixToFileName, combine, separateFileList, silentRemove,extractSum
from drUtils import subtractOverscan, subtractBias, cleanCosmic, flatCorrect,interpolateTraceIm
from drUtils import makeSkyFlat, makeMasterFlat, imDivide, extractAndReidentifyARCs, dispCor
from drUtils import readFluxStandardsList,calcResponse,applySensFuncs
import numpy as np
import os
from shutil import copyfile

overscanSection = '[1983:,:]'
trimSection = '[17:1982,38:97]'
#workPath = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190904/'
workPath = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190905/'
refPath = '/Users/azuri/stella/referenceFiles/spupnic'
#workPath = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190506/'

#refVerticalTraceDB = '/Users/azuri/stella/referenceFiles/database/spupnic/apvertical_trace'
refVerticalTraceDB = os.path.join(refPath,'database/aprefVerticalTrace_spupnic_gr7_16_3')
#refHorizontalTraceDB = '/Users/azuri/stella/referenceFiles/database/spupnic/aphorizontal_tracer90flipl'
refHorizontalTraceDB = os.path.join(refPath,'database/aprefHorizontalTrace_spupnic_gr7_16_3_transposed')
refProfApDef = os.path.join(refPath,'database/aprefProfApDef_spupnic_gr7_16_3')
lineList = os.path.join(refPath,'saao_refspec_gr7_angle16_3_lines_identified_good.dat')

print('EarthLocation.get_site_names() = ',EarthLocation.get_site_names())
observatoryLocation = EarthLocation.of_site('SAAO')
print('SAAO.lon = ',observatoryLocation.lon)
print('SAAO.lat = ',observatoryLocation.lat)

fluxStandardNames, fluxStandardDirs, fluxStandardFileNames = readFluxStandardsList()

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
suffixes = ['','ot','otz','otzf','otzfi','otzfif','otzx','otzxf','otzxfi','otzxfif','otzfiEc','otzfifEc','otzfifEcd','otzfifEcdF']

#    copyfile(inList, inList+'bak')
#    silentRemove(inList[:inList.rfind('/')+1]+'*.list')
#    copyfile(inList+'bak', inList)
exptypes = ['BIAS','FLAT','ARC','SCIENCE','FLUXSTDS']
objects = [['*'],['*','DOMEFLAT','SKYFLAT'],['*'],['*','individual'],['*']]
if False:
    separateFileList(inList, suffixes, exptypes, objects, True, fluxStandardNames=fluxStandardNames)
    #STOP
    objectFiles = os.path.join(workPath,'SCIENCE.list')
if False:
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
    flat = combine(getListOfFiles(os.path.join(workPath,'FLATDOMEFLAT_otz.list')),
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
    for inputList in ['ARC','SCIENCE','FLATSKYFLAT']:
        flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                    masterFlat,
                    norm_value = 1.,
                    fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzf.list')))
if False:
    # interpolate images to get straight dispersion and spectral features
    for inputList in ['ARC', 'SCIENCE']:
        interpolateTraceIm(getListOfFiles(os.path.join(workPath,inputList+'_otzf.list')),
                           refVerticalTraceDB,
                           refHorizontalTraceDB)

    # create master SkyFlat
    combinedSkyFlat = os.path.join(workPath,'combinedSkyFlat.fits')
    print('creating combinedSkyFlat <'+combinedSkyFlat+'>')
    combine(getListOfFiles(os.path.join(workPath,'FLATSKYFLAT_otzf.list')),
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
                        refVerticalTraceDB,
                        refHorizontalTraceDB)

    makeSkyFlat(os.path.join(workPath,'combinedSkyFlati.fits'),
                os.path.join(workPath,'combinedSkyFlati_flattened.fits'),
                7)

    flatCorrect(getListOfFiles(os.path.join(workPath,'SCIENCE_otzfi.list')),
                os.path.join(workPath,'combinedSkyFlati_flattened.fits'),
                fitsFilesOut=getListOfFiles(os.path.join(workPath,'SCIENCE_otzfif.list')))

if True:
    # extract and reidentify ARCs
    wavelengthsOrig, wavelengthsResampled = extractAndReidentifyARCs(getListOfFiles(os.path.join(workPath,'ARC_otzf.list')),
                                                                     refProfApDef,
                                                                     lineList)

    # extract science data
    scienceSpectra = [extractSum(fn,'row',fn[:-5]+'Ec.fits') for fn in getListOfFiles(os.path.join(workPath,'SCIENCE_otzfif.list'))]

    dispCor(getListOfFiles(os.path.join(workPath,'SCIENCE_otzfifEc.list')),
            getListOfFiles(os.path.join(workPath,'ARC_otzfiEc.list')),
            wavelengthsOrig,
            getListOfFiles(os.path.join(workPath,'SCIENCE_otzfifEcd.list')),
            observatoryLocation,
            'TARG-RA',
            'TARG-DEC',
            'DATE-OBS')

    sensFuncs = calcResponse(os.path.join(workPath,'FLUXSTDS_otzfif.list'), wavelengthsOrig[0])
    print('sensFuncs = ',sensFuncs)

    applySensFuncs(getListOfFiles(os.path.join(workPath,'SCIENCE_otzfifEcd.list')),getListOfFiles(os.path.join(workPath,'SCIENCE_otzfifEcdF.list')), sensFuncs)

#    response = calcResponse(os.path.join(workPath,'FLUXSTDS_otzfifEcd.list'))
#    fluxCalibrate(os.path.join(workPath,'SCIENCE_LTT7379_a1171120_otzfifEcd.fits'),
#                  os.path.join(workPath,'SCIENCE_Feige110_a1171214_otzfifEcd.fits'))
