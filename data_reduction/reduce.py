#from astropy.nddata import CCDData
from astropy.coordinates import EarthLocation
from drUtils import addSuffixToFileName, combine, separateFileList, silentRemove,extractSum
from drUtils import subtractOverscan, subtractBias, cleanCosmic, flatCorrect,interpolateTraceIm
from drUtils import makeSkyFlat, makeMasterFlat, imDivide, extractAndReidentifyARCs, dispCor
from drUtils import readFluxStandardsList,calcResponse,applySensFuncs,extractObjectAndSubtractSky
from drUtils import scombine,continuum,subtractMedianSky,removeFilesFromListWithAngleNotEqualTo
import numpy as np
import os
from shutil import copyfile

import csvFree, csvData

overscanSection = '[1983:,:]'
trimSection = '[17:1982,38:97]'
#workPath = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190904/'
#workPath = '/Users/azuri/spectra/saao/saao_sep2019/20190907/'
workPath = '/Users/azuri/spectra/saao/saao_aug2018/20180811/'
refPath = '/Users/azuri/stella/referenceFiles/spupnic'
#workPath = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190506/'

#refVerticalTraceDB = '/Users/azuri/stella/referenceFiles/database/spupnic/apvertical_trace'
#refVerticalTraceDB = os.path.join(refPath,'database/aprefVerticalTrace_spupnic_gr7_16_3')
refVerticalTraceDB = os.path.join(refPath,'database/aprefVerticalTrace_spupnic_gr7_15_90')
#refHorizontalTraceDB = '/Users/azuri/stella/referenceFiles/database/spupnic/aphorizontal_tracer90flipl'
refHorizontalTraceDB = os.path.join(refPath,'database/aprefHorizontalTrace_spupnic_gr7_16_3_transposed')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_spupnic_gr7_15_85')#16_3')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_spupnic_gr7_16_3')
refProfApDef = os.path.join(refPath,'database/aprefProfApDef_spupnic_gr7_%d_%d')
#lineList = os.path.join(refPath,'saao_refspec_gr7_angle16_3_lines_identified_good.dat')
#lineList = os.path.join(refPath,'saao_refspec_gr7_angle16_3_may2020_lines_identified_good.dat')
#lineList = os.path.join(refPath,'saao_refspec_gr7_angle15_85_lines_identified_good.dat')
lineList = os.path.join(refPath,'saao_refspec_gr7_angle%d_%d_lines_identified_good_aug2018.dat')

print('EarthLocation.get_site_names() = ',EarthLocation.get_site_names())
observatoryLocation = EarthLocation.of_site('SAAO')
print('SAAO.lon = ',observatoryLocation.lon)
print('SAAO.lat = ',observatoryLocation.lat)

fluxStandardNames, fluxStandardDirs, fluxStandardFileNames = readFluxStandardsList()

if not os.path.exists(os.path.join(workPath,'database')):
    os.makedirs(os.path.join(workPath,'database'))

def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

def getListOfFiles(fname):
    fList = readFileToArr(fname)
    print('fname = ',fname,': fList = ',fList)
    if len(fList) > 0:
        if fList[0].rfind('/') == -1:
            fList = [os.path.join(workPath, fileName) for fileName in fList]
    return fList

inList=os.path.join(workPath,'allFits.list')
suffixes = ['','ot','otz','otzf','otzfi','otzfif','otzx','otzxf','otzxfi','otzxfif','otzfiEc','otzxfiEc','otzfifEc','otzxfifEc','otzfifEcd','otzxfifEcd','otzfifEcdF','otzxfifEcdF']

#    copyfile(inList, inList+'bak')
#    silentRemove(inList[:inList.rfind('/')+1]+'*.list')
#    copyfile(inList+'bak', inList)
exptypes = ['BIAS','FLAT','ARC','SCIENCE','FLUXSTDS']
objects = [['*'],['*','DOMEFLAT','SKYFLAT'],['*'],['*','individual'],['*']]
if False:
#    removeFilesFromListWithAngleNotEqualTo(inList,inList,'15.85')
#    removeFilesFromListWithAngleNotEqualTo(inList,inList,'15.85')
#    STOP
    separateFileList(inList, suffixes, exptypes, objects, True, fluxStandardNames=fluxStandardNames)
#    STOP
    objectFiles = os.path.join(workPath,'SCIENCE.list')
    # subtract overscan and trim all images
    for inputList in ['ARC', 'BIAS', 'FLAT', 'SCIENCE','FLUXSTDS']:
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
    for inputList in ['ARC', 'FLAT', 'SCIENCE','FLUXSTDS']:
        subtractBias(getListOfFiles(os.path.join(workPath,inputList+'_ot.list')),
                     masterBias,
                     fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                     overwrite=True)
if False:
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

    # remove cosmic rays
    for inputList in ['ARC', 'SCIENCE','FLUXSTDS']:
        cosmicParameters = {'niter':3, 'sigclip':4.0, 'objlim':5.0, 'gain':1.145, 'readnoise':2.245, 'sepmed':False, 'cleantype':'meanmask', 'fsmode':'convolve', 'psfmodel':'gaussy', 'psffwhm':0.75}
        #cosmicParameters = {'thresh':5., 'rbox':3}
        cleanCosmic(getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                    #cosmicMethod='lacosmic',
                    #cosmicParameters=cosmicParameters,
                    fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzx.list')),
                    overwrite=True)

    # apply master DomeFlat to ARCs, SkyFlats, and SCIENCE frames
    for inputList in ['ARC','SCIENCE','FLUXSTDS']:
        flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                    masterFlat,
                    norm_value = 1.,
                    fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzf.list')))
    for inputList in ['ARC','SCIENCE','FLUXSTDS']:
        flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otzx.list')),
                    masterFlat,
                    norm_value = 1.,
                    fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzxf.list')))
    for inputList in ['FLATSKYFLAT']:
        flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                    masterFlat,
                    norm_value = 1.,
                    fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzf.list')))
if False:
    # interpolate images to get straight dispersion and spectral features
    for inputList in ['ARC', 'SCIENCE','FLUXSTDS']:
        interpolateTraceIm(getListOfFiles(os.path.join(workPath,inputList+'_otzxf.list')),
                           refVerticalTraceDB,
                           refHorizontalTraceDB)
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

    for inputList in ['ARC', 'SCIENCE','FLUXSTDS']:
        flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otzfi.list')),
                    os.path.join(workPath,'combinedSkyFlati_flattened.fits'),
                    fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzfif.list')))
        flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otzxfi.list')),
                    os.path.join(workPath,'combinedSkyFlati_flattened.fits'),
                    fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzxfif.list')))

if True:
    subtractMedianSky(getListOfFiles(os.path.join(workPath,'SCIENCE_otzfif.list')))
    subtractMedianSky(getListOfFiles(os.path.join(workPath,'SCIENCE_otzxfif.list')))

    # extract and reidentify ARCs
    wavelengthsOrig, wavelengthsResampled = extractAndReidentifyARCs(getListOfFiles(os.path.join(workPath,'ARC_otzxf.list')),
                                                                     refProfApDef,
                                                                     lineList,
                                                                     display=True)

    # extract arcs and science data
    #extractObjectAndSubtractSky(twoDImageFileIn, specOut, yRange, skyAbove, skyBelow, dispAxis)
    inputList = getListOfFiles(os.path.join(workPath,'ARC_otzxfif.list'))
    for inputFile in inputList:
        extractSum(inputFile, 'row', fNameOut = inputFile[:inputFile.rfind('.')]+'Ec.fits')

    areasFileName = os.path.join(workPath,'areas.csv')
    print('reduce: areasFileName = '+areasFileName)
    areas = csvFree.readCSVFile(areasFileName)
    for i in range(areas.size()):
        print('reduce: areas.getData("fName",',i,') = '+areas.getData('fName',i),', object = ',areas.getData('object',i),', skyBelow = ',areas.getData('skyBelow',i),', skyAbove = ',areas.getData('skyAbove',i,),', method = ',areas.getData('method',i),', notes = ',areas.getData('notes',i))
    scienceSpectra = []#extractSum(fn,'row',fn[:-5]+'Ec.fits') for fn in getListOfFiles(os.path.join(workPath,'SCIENCE_otzfif.list'))]
    for i in range(areas.size()):
        if areas.getData('fName',i)[0] != '#':
            skyAbove = None if areas.getData('skyAbove',i) == '' else [int(areas.getData('skyAbove',i)[1:areas.getData('skyAbove',i).find(':')]),int(areas.getData('skyAbove',i)[areas.getData('skyAbove',i).find(':')+1:-1])]
            skyBelow = None if areas.getData('skyBelow',i) == '' else [int(areas.getData('skyBelow',i)[1:areas.getData('skyBelow',i).find(':')]),int(areas.getData('skyBelow',i)[areas.getData('skyBelow',i).find(':')+1:-1])]
            print('reduce: i = ',i,': fName = ',areas.getData('fName',i))
            extractObjectAndSubtractSky(os.path.join(workPath,areas.getData('fName',i)) if '/' not in areas.getData('fName',i) else areas.getData('fName',i),
                                        os.path.join(workPath,areas.getData('fName',i)[:-5]+'Ec.fits'),
                                        [int(areas.getData('object',i)[1:areas.getData('object',i).find(':')]),int(areas.getData('object',i)[areas.getData('object',i).find(':')+1:-1])],
                                        skyAbove = skyAbove,
                                        skyBelow = skyBelow,
                                        extractionMethod = areas.getData('method',i),
                                        dispAxis = 'row',
                                        display = False)

    for inputList in ['ARC', 'SCIENCE']:
        doHelioCor = True if inputList == 'SCIENCE' else False
        print('doHelioCor = ',doHelioCor)
        dispCor(getListOfFiles(os.path.join(workPath,inputList+'_otzxfifEc.list')),
                getListOfFiles(os.path.join(workPath,'ARC_otzxfiEc.list')),
                wavelengthsOrig,
                getListOfFiles(os.path.join(workPath,inputList+'_otzxfifEcd.list')),
                observatoryLocation,
                'TELRA',
                'TELDEC',
                'DATE-OBS',
                doHelioCor = doHelioCor)

    areas = csvFree.readCSVFile(os.path.join(workPath,'areas.csv'))
    sensFuncs = calcResponse(os.path.join(workPath,'FLUXSTDS_otzxfif.list'),
                             getListOfFiles(os.path.join(workPath,'ARC_otzxfiEc.list')),
                             wavelengthsOrig,
                             areas)
    print('sensFuncs = ',sensFuncs)

    applySensFuncs(getListOfFiles(os.path.join(workPath,'SCIENCE_otzxfifEcd.list')),
                   getListOfFiles(os.path.join(workPath,'SCIENCE_otzxfifEcdF.list')),
                   sensFuncs)

if True:
    # combine multiple spectra of the same object
    fileList = getListOfFiles(os.path.join(workPath,'SCIENCE_otzxfifEcdF.list'))
    toCombineLists = []
    combinedList = []
    toContinuumSubtract = []
    continuumSubtracted = []
    for fileName in fileList:
        objectName = fileName[fileName.rfind('SCIENCE_')+8:]
        objectName = objectName[:objectName.find('_a')]
        print('objectName = ',objectName)
        finalName = fileName[:fileName.rfind('/')+1]+objectName+'_SA'+fileName[fileName.rfind('/')-2:fileName.rfind('/')]+fileName[fileName.rfind('/')-4:fileName.rfind('/')-2]+fileName[fileName.rfind('/')-6:fileName.rfind('/')-4]+'.fits'
        print('finalName = <'+finalName+'>')

        nFiles = 0
        for fName in fileList:
            objName = fName[fName.rfind('SCIENCE_')+8:]
            objName = objName[:objName.find('_a')]
            print('objName = ',objName)
            if objectName == objName:
                nFiles += 1
        if nFiles > 1:
            listNameOut = fileName[:fileName.rfind('_a')]+'_toCombine.list'
            if listNameOut not in toCombineLists:
                toCombineLists.append(listNameOut)
                combinedList.append(listNameOut[:listNameOut.rfind('toCombine')]+'combined.fits')
                toContinuumSubtract.append(combinedList[len(combinedList)-1])
                continuumSubtracted.append(finalName)
            print('nFiles = ',nFiles,': creating <'+listNameOut+'>')
            with open(listNameOut,'w') as f:
                for fName in fileList:
                    objName = fName[fName.rfind('SCIENCE_')+8:]
                    objName = objName[:objName.find('_a')]
                    print('objName = ',objName)
                    if objectName == objName:
                        f.write(fName+'\n')
        else:
            toContinuumSubtract.append(fileName)
            continuumSubtracted.append(finalName)
    print('toCombineLists = ',toCombineLists)
    for iList in range(len(toCombineLists)):
        inputList = toCombineLists[iList]
        outFile = combinedList[iList]
        scombine(inputList,
                 outFile,
                 method='median',
                 lowReject=2.,
                 highReject=2.,
                 adjustSigLevels=False,
                 useMean=False,
                 display=True)

    # subtract continuum
    for i in range(len(toContinuumSubtract)):
        fittingFunction = np.polynomial.legendre.legfit
        evalFunction = np.polynomial.legendre.legval
        order = 9
        nIterReject = 2
        nIterFit = 2
        lowReject = 2.
        highReject = 2.
        useMean = False
        continuum(toContinuumSubtract[i],
                  continuumSubtracted[i],
                  fittingFunction,
                  evalFunction,
                  order,
                  nIterReject,
                  nIterFit,
                  lowReject,
                  highReject,
                  type='difference',
                  adjustSigLevels=False,
                  useMean=useMean,
                  display=True)



#    response = calcResponse(os.path.join(workPath,'FLUXSTDS_otzfifEcd.list'))
#    fluxCalibrate(os.path.join(workPath,'SCIENCE_LTT7379_a1171120_otzfifEcd.fits'),
#                  os.path.join(workPath,'SCIENCE_Feige110_a1171214_otzfifEcd.fits'))
