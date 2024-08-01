from astropy.coordinates import EarthLocation
import astropy.io.fits as pyfits
from drUtils import addSuffixToFileName, combine, separateFileList, silentRemove,extractSum
from drUtils import subtractOverscan, subtractBias, cleanCosmic, flatCorrect,interpolateTraceIm
from drUtils import makeSkyFlat, makeMasterFlat, imDivide, extractAndReidentifyARCs, dispCor
from drUtils import readFluxStandardsList,calcResponse,applySensFuncs,extractObjectAndSubtractSky
from drUtils import scombine,continuum,subtractMedianSky,removeFilesFromListWithAngleNotEqualTo
from drUtils import getWavelengthArr,getListOfFiles,getHeaderValue,fixDBSHeaders,writeFits1D,invertY
from drUtils import cleanSpec,setHeaderValue
import numpy as np
import os
from shutil import copyfile

import csvFree, csvData

#trimSection = '[51:2148,18:507]'#'[26:1774,30:115]'#'[17:1982,38:97]'
#trimSection = '[56:2100,135:454]'#DBS May2008 Blue
#trimSection = '[56:2100,133:452]'#DBS May2008 Blue
#trimSection = '[56:2100,75:394]'#'[26:1774,30:115]'#'[17:1982,38:97]'

"""SPUPNIC grating #7, 16.3 deg grating angle"""
observatoryLocation = EarthLocation.of_site('SAAO')
workPaths = ['/Users/azuri/spectra/SAAO_June2024/110624/',]
refPath = '/Users/azuri/stella/referenceFiles/spupnic'

overscanSection = '[1983:,:]'
trimSection = '[16:1981,41:100]'

refVerticalTraceDB = os.path.join(refPath,'database/aprefVerticalTrace_spupnic_gr7_16_3_2024_otzf')#refVerticalTrace_spupnic_gr7_16_3')
refHorizontalTraceDB = os.path.join(refPath,'database/aprefHorizontalTrace_spupnic_gr7_16_3_transposed')
refProfApDef = os.path.join(refPath,'database/aprefProfApDef_spupnic_gr7_16_3')#refProfApDef_spupnic_gr7_16_3_may2020')
lineList = os.path.join(refPath,'refArc_spupnic_gr7_16_30_otzxfifEcd_lines_identified.dat')#saao_refspec_gr7_angle16_3_may2020_lines_identified_good.dat')
referenceSpectrum = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_gr7_16_30_otzfsiEc.fits'

arcListsStartWith = 'arc'#DBS
objectListsStartWith = 'science'#spupnic
fluxstdsListsStartWith = 'fluxstds'#DBS

"""MSSO May 2008"""
#observatoryLocation = EarthLocation.of_site('Mt. Stromlo Observatory')
#workPaths = ['/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-08/',]
#refPath = '/Users/azuri/stella/referenceFiles/dbs'

#MSSSOMay2008 overscanSection = '[10:50,1:562]'#'[4:21,1:133]'#'[1983:,:]'
#MSSSOMay2008 trimSection = '[51:2101,25:367]'#'[26:1774,30:115]'#'[17:1982,38:97]'
#referenceSpectrum = os.path.join(refPath,'refArc_DBS_May2008_red_otzxfiEc.fits')
#referenceSpectrum = os.path.join(refPath,'refARC_DBS_May2008_blue_otzfiEc.fits')
#lineList = os.path.join(refPath,'dbs_refspec_may2008_red.txt')
#lineList = os.path.join(refPath,'dbs_refspec_may2008_blue_CuHe.txt')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_DBS_red_May2008_otzfif')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_DBS_blue_May2008_otzfif_flipped')
#refHorizontalTraceDB = os.path.join(refPath,'database/aprefApTrace_DBS_RED_horizontal_otzxf')#aprefHorizontalTrace_spupnic_2007_transposed')
#refHorizontalTraceDB = os.path.join(refPath,'database/aprefApTrace_horizontal_DBS_blue_May2008_otzf_flipped')#aprefHorizontalTrace_spupnic_2007_transposed')
#refVerticalTraceDB = os.path.join(refPath,'database/aprefApTrace_vertical_DBS_red_May2008_otzf')
#refVerticalTraceDB = os.path.join(refPath,'database/aprefApTrace_vertical_DBS_blue_May2008_otzf_flipped')
#arcListsStartWith = 'ARC'#spupnic
#objectListsStartWith = 'object'#DBS
#fluxstdsListsStartWith = 'FLUXSTDS'#spupnic

#workPath = '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190904/'
#workPath = '/Users/azuri/spectra/saao/saao_sep2019/20190907/'
#workPath = '/Users/azuri/spectra/saao/saao_may2007/RAW/070512/'
#workPaths = ['/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-14/',
#             ]
#workPaths = ['/Users/azuri/spectra/MSSSO_2m3_DBS_aug07/RED/night7/',]
if False:
    workPaths = ['/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-06/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-07/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-08/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-09/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-10/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-11/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-12/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-13/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-14/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/2008-05-15/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-06/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-07/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-08/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-08/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-09/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-10/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-11/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-12/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-13/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-14/',
             '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/2008-05-15/',
            ]
#refPath = '/Users/azuri/stella/referenceFiles/dbs'#spupnic'
#workPath = '/Volumes/work/azuri/spectra/saao/saao_may2019/20190506/'


#refVerticalTraceDB = '/Users/azuri/stella/referenceFiles/database/spupnic/apvertical_trace'
#refVerticalTraceDB = os.path.join(refPath,'database/aprefVerticalTrace_spupnic_gr7_15_90')
#refVerticalTraceDB = os.path.join(refPath,'database/aprefApTrace_DBS_RED_vertical_otzxf')#aprefVerticalTrace_spupnic_2007')
#refVerticalTraceDB = os.path.join(refPath,'database/aprefApTrace_DBS_RED_6230-6725_vertical_otzxf')#aprefVerticalTrace_spupnic_2007')
#refVerticalTraceDB = os.path.join(refPath,'database/aprefApTrace_DBS_BLUE_vertical_otzf')#aprefVerticalTrace_spupnic_2007')
#refVerticalTraceDB = os.path.join(refPath,'database/aprefApTrace_DBS_BLUE_3600-5600_vertical_otzf')#aprefVerticalTrace_spupnic_2007')

#refHorizontalTraceDB = '/Users/azuri/stella/referenceFiles/database/spupnic/aphorizontal_tracer90flipl'
#refHorizontalTraceDB = os.path.join(refPath,'database/aprefApTrace_DBS_RED_6230-6725_horizontal_otzxf')#aprefHorizontalTrace_spupnic_2007_transposed')
#refHorizontalTraceDB = os.path.join(refPath,'database/aprefApTrace_DBS_BLUE_horizontal_otzxf')#aprefHorizontalTrace_spupnic_2007_transposed')
#refHorizontalTraceDB = os.path.join(refPath,'database/aprefApTrace_DBS_BLUE_3600-5600_horizontal_otzf')#aprefHorizontalTrace_spupnic_2007_transposed')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_spupnic_gr7_15_85')#16_3')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_spupnic_gr7_%d_%d_aug2013')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_DBS_RED_otzfif')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_DBS_RED_6230-6725_otzxfif')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_DBS_BLUE_otzfif')
#refProfApDef = os.path.join(refPath,'database/aprefProfApDef_DBS_aug07_Blue_3600-5600_otzfif')

#lineList = os.path.join(refPath,'saao_refspec_gr7_angle16_3_lines_identified_good.dat')
#lineList = os.path.join(refPath,'saao_refspec_gr7_angle15_85_lines_identified_good.dat')
#lineList = os.path.join(refPath,'saao_refspec_gr7_angle%d_%d_lines_identified_good_aug2018.dat')
#lineList = os.path.join(refPath,'saao_refspec_gr7_angle%d_%d_lines_identified_good_mar2014.dat')
#lineList = os.path.join(refPath,'saao_refspec_gr7_angle%d_%d_lines_identified_good_aug2013.dat')
#lineList = os.path.join(refPath,'saao_refspec_may2007.new')
#lineList = os.path.join(refPath,'dbs_refspec_aug2007_CuAr_blue.txt')
#lineList = os.path.join(refPath,'dbs_refspec_aug2007_red.txt')
#lineList = os.path.join(refPath,'dbs_refspec_aug2007_CuAr_red_6230-6725_lines.txt')
#lineList = os.path.join(refPath,'dbs_refspec_aug2007_CuAr_blue_3600-5600_lines.txt')

#referenceSpectrum = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_gr7_15_70_otzxfifEc_aug2018.fits'
#referenceSpectrum = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_otzxfifEc_mar2014.fits'
#referenceSpectrum = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_may2007.fits'
#referenceSpectrum = os.path.join(refPath,'refArc_DBS_aug07_Blue.fits')
#referenceSpectrum = os.path.join(refPath,'refArc_DBS_aug07_Red.fits')
#referenceSpectrum = os.path.join(refPath,'refArc_DBS_RED_62360-6725_otzfifEc.fits')
#referenceSpectrum = os.path.join(refPath,'refArc_DBS_aug07_Blue_3600-5600.fits')

print('EarthLocation.get_site_names() = ',EarthLocation.get_site_names())
print('obsSite.lon = ',observatoryLocation.lon)
print('obsSite.lat = ',observatoryLocation.lat)

fluxStandardNames, fluxStandardDirs, fluxStandardFileNames = readFluxStandardsList()

for workPath in workPaths:
    if not os.path.exists(os.path.join(workPath,'database')):
        os.makedirs(os.path.join(workPath,'database'))

def fixHeaders():
    with open(os.path.join(workPath,'allFits.list'),'r') as f:
        lines = f.readlines()
    for line in lines:
        fName = line.strip()
        print('fixHeaders: fName = ',fName)
        hdulist = pyfits.open(os.path.join(workPath,fName))
        typeName = 'EXPTYPE'
        try:
            expType = hdulist[0].header['EXPTYPE']
        except:
            try:
                expType = hdulist[0].header['IMAGETYP']
                typeName = 'IMAGETYP'
            except:
                print('neither EXPTYPE nor IMAGETYP found in header of file '+fName)
                STOP
        if expType == 'COMPARISON':
            hdulist[0].header[typeName] = 'ARC'
        elif expType == 'flat':
            hdulist[0].header[typeName] = 'FLAT'
        elif expType == 'zero':
            hdulist[0].header[typeName] = 'BIAS'
        elif expType == 'object':
            hdulist[0].header[typeName] = 'SCIENCE'

        hdulist.writeto(os.path.join(workPath,fName), overwrite=True)

def fixHeadersDBS():
    with open(os.path.join(workPath,'allFits.list'),'r') as f:
        lines = f.readlines()
    for line in lines:
        fName = line.strip()
        print('fixHeaders: fName = ',fName)
        hdulist = pyfits.open(os.path.join(workPath,fName))
        typeName = 'EXPTYPE'
        try:
            expType = hdulist[0].header['EXPTYPE']
        except:
            try:
                expType = hdulist[0].header['IMAGETYP']
                typeName = 'IMAGETYP'
            except:
                print('neither EXPTYPE nor IMAGETYP found in header of file '+fName)
                STOP
        if expType == 'COMPARISON':
            hdulist[0].header[typeName] = 'ARC'
        elif expType == 'flat':
            hdulist[0].header[typeName] = 'FLAT'
        elif expType == 'zero':
            hdulist[0].header[typeName] = 'BIAS'
        elif expType == 'object':
            hdulist[0].header[typeName] = 'SCIENCE'

        object = hdulist[0].header['OBJECT']
        if object == 'DOME':
            hdulist[0].header['OBJECT'] = 'DOMEFLAT'
        if hdulist[0].header['OBJECT'] == 'SKY':
            hdulist[0].header['OBJECT'] = 'SKYFLAT'
        if hdulist[0].header['OBJECT'] == 'DOMEFLAT':
            hdulist[0].header[typeName] = 'FLAT'
        if object == 'CuHe':
            hdulist[0].header[typeName] = 'ARC'
        if object == 'BIAS':
            hdulist[0].header[typeName] = 'zero'

        hdulist.writeto(os.path.join(workPath,fName), overwrite=True)

def fixHeadersSAAO():
    with open(os.path.join(workPath,'allFits.list'),'r') as f:
        fNames = f.readlines()
    fNames = [fName.strip() for fName in fNames]
    #obsLogName = os.path.join(workPath,'obslog.txt')
    #obslog = csvFree.readCSVFile(obsLogName)
    #for i in range(obslog.size()):
    #    obslog.setData('FRAME',i,str(1000+i+1))
    #    if obslog.getData('EXPTYPE',i) == 'ARC':
    #        obslog.setData('OBJECT',i,'CuAr')
    #frameIDs = np.array(obslog.getData('FRAME'))
    #print('frameIDs = ',frameIDs)
    #csvFree.writeCSVFile(obslog,obsLogName)
    for fName in fNames:
        #frameID = str(getHeaderValue(fName,'FRAME'))
        #print('fName = ',fName,': frameID = ',type(frameID),' = <'+frameID+'>')
        #idx = np.where(frameIDs == frameID)[0]
        #print('idx = ',idx)
        if getHeaderValue(fName,'EXPTYPE') == 'ARC':
            setHeaderValue(fName,'OBJECT','CuAr')
        if getHeaderValue(fName,'EXPTYPE') == 'HARTMANN':
            setHeaderValue(fName,'OBJECT','HARTMANN')
#        setHeaderValue(fName,'EXPTYPE',obslog.getData('EXPTYPE',idx))
#        command = 'cp %s %s' % (fName,fName[:fName.rfind('/')+1]+getHeaderValue(fName,'EXPTYPE')+'_'+getHeaderValue(fName,'OBJECT')+'_'+fName[fName.rfind('/')+1:])
#        os.system(command)


def fixWCSDim():
    lists = [os.path.join(workPath,'science_otzxfifEcd.list'),os.path.join(workPath,'science_otzxfifEcdF.list'),os.path.join(workPath,'arc_otzxfiEcd.list')]
    for myList in lists:
        with open(myList,'r') as f:
            lines = f.readlines()
        for line in lines:
            fName = line.strip()
            print('fixWCSDim: fName = ',fName)
            hdulist = pyfits.open(fName)
            hdulist[0].header['WCSDIM'] = 1

            hdulist.writeto(fName, overwrite=True)


def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

if __name__ == '__main__':
    #fixHeaders()
    for workPath in workPaths:
        inList=os.path.join(workPath,'allFits.list')
        suffixes = ['','ot','otz','otzf','otzfi','otzfif','otzx','otzxf','otzxfi','otzxfif','otzfiEc','otzxfiEc','otzfifEc','otzxfifEc','otzfifEcd','otzxfifEcd','otzfifEcdF','otzxfifEcdF']

        #    copyfile(inList, inList+'bak')
        #    silentRemove(inList[:inList.rfind('/')+1]+'*.list')
        #    copyfile(inList+'bak', inList)
        """show image types"""
        if False:
            fixHeadersSAAO()# and rename file

        if False:
            with open(inList,'r') as f:
                fitsFiles = f.readlines()
            fitsFiles = [f.strip() for f in fitsFiles]
            for f in fitsFiles:
                imType = getHeaderValue(f,'EXPTYPE')
                if imType == "HARTMANN":
                    setHeaderValue(f,'OBJECT','HARTMANN')
                object = getHeaderValue(f,'OBJECT')
                print(imType,': ',object)
#        STOP
        if True:
            exptypes = ['BIAS','FLAT','ARC','SCIENCE','FLUXSTDS']
            #exptypes = ['ZERO','FLAT','ARC','SCIENCE','FLUXSTDS']
            #exptypes = ['zero','flat','COMPARISON','SCIENCE','FLUXSTDS']
            objects = [['*'],['*','DOMEFLAT','SKYFLAT'],['*'],['*','individual'],['*']]
            masterBias = os.path.join(workPath,'combinedBias_ot.fits')
            combinedFlat = os.path.join(workPath,'combinedFlat.fits')
            masterFlat = os.path.join(workPath, 'masterDomeFlat.fits')
            smoothedFlat = os.path.join(workPath, 'smoothedDomeFlat.fits')
        if False:
            #invertY(readFileToArr(inList))#MSSSO only
            #STOP
        #    removeFilesFromListWithAngleNotEqualTo(inList,inList,'15.85')
        #    removeFilesFromListWithAngleNotEqualTo(inList,inList,'15.85')
        #    STOP
            separateFileList(inList, suffixes, exptypes, objects, True, fluxStandardNames=fluxStandardNames)
        #STOP
        if False:
            objectFiles = os.path.join(workPath,'science.list')
    #        objectFiles = os.path.join(workPath,'SCIENCE.list')
            # subtract overscan and trim all images
            for inputList in ['arc', 'bias', 'flat', 'science','fluxstds']:
                subtractOverscan(getListOfFiles(os.path.join(workPath,inputList+'.list')),
                                overscanSection,
                                trimSection=trimSection,
                                fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_ot.list')),
                                overwrite=True)

            # create master bias
            combinedBias = combine(getListOfFiles(os.path.join(workPath,'bias_ot.list')),
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
            for inputList in ['arc', 'flat', 'science','fluxstds']:
                subtractBias(getListOfFiles(os.path.join(workPath,inputList+'_ot.list')),
                            masterBias,
                            fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                            overwrite=True)
            # create master DomeFlat
            print('creating combinedFlat <'+combinedFlat+'>')
            flat = combine(getListOfFiles(os.path.join(workPath,'flatDOMEFLAT_otz.list')),
                        combinerMethod='median',
                        clippingMethod='sigma',
                        clippingParameters={'niter':2,
                                            'low_thresh':-3.,
                                            'high_thresh':3.,
                                            'func':np.ma.median},
                        scaling=False,
                        minVal=0.0001,
                        fitsOutName=combinedFlat)
            makeMasterFlat(combinedFlat,
                        9,
                        80.,
                        outFileNameMasterFlat=masterFlat,
                        outFileNameMasterFlatSmoothed=smoothedFlat)
            # remove cosmic rays
            for inputList in ['arc', 'science','fluxstds']:
                cosmicParameters = {'niter':2,
                                    'sigclip':10.,
                                    'objlim':5.0,
                                    'gain':1.0,
                                    'readnoise':15.,
                                    'sepmed':False,
                                    'cleantype':'meanmask',
                                    'fsmode':'convolve',
                                    'psfmodel':'gaussy',
                                    'psffwhm':2.5}
                #cosmicParameters = {'thresh':5., 'rbox':3}
                cleanCosmic(getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                            cosmicMethod='lacosmic',
                            cosmicParameters=cosmicParameters,
                            fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzx.list')),
                            overwrite=True)
            #STOP
            # apply master DomeFlat to ARCs, SkyFlats, and SCIENCE frames
        if False:
            for inputList in ['arc','science','fluxstds','flatSKYFLAT']:
                flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otz.list')),
                            masterFlat,
                            norm_value = 1.,
                            fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzf.list')))
            for inputList in ['arc','science','fluxstds']:
                flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otzx.list')),
                            masterFlat,
                            norm_value = 1.,
                            fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzxf.list')))
        if False:
            # interpolate images to get straight dispersion and spectral features
            for inputList in ['arc', 'science','fluxstds']:
                interpolateTraceIm(getListOfFiles(os.path.join(workPath,inputList+'_otzxf.list')),
                                refVerticalTraceDB,
                                refHorizontalTraceDB)
                interpolateTraceIm(getListOfFiles(os.path.join(workPath,inputList+'_otzf.list')),
                                refVerticalTraceDB,
                                refHorizontalTraceDB)

        if False:
            # create master SkyFlat
            with open(os.path.join(workPath,'flatSKYFLAT_otzf.list'),'r') as fl:
                skyFlats = fl.readlines()
            haveSkyFlats = True if len(skyFlats) > 0 else False
            if haveSkyFlats:
                combinedSkyFlat = os.path.join(workPath,'combinedSkyFlat.fits')
                print('creating combinedSkyFlat <'+combinedSkyFlat+'>')
                combine(getListOfFiles(os.path.join(workPath,'flatSKYFLAT_otzf.list')),
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
                for inputList in ['arc', 'science','fluxstds']:
                    flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otzfi.list')),
                                os.path.join(workPath,'combinedSkyFlati_flattened.fits'),
                                fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzfif.list')))
                    flatCorrect(getListOfFiles(os.path.join(workPath,inputList+'_otzxfi.list')),
                                os.path.join(workPath,'combinedSkyFlati_flattened.fits'),
                                fitsFilesOut=getListOfFiles(os.path.join(workPath,inputList+'_otzxfif.list')))
            else:
                for inputList in ['arc', 'science','fluxstds']:
                    for inputFile,outputFile in zip(getListOfFiles(os.path.join(workPath,inputList+'_otzfi.list')),getListOfFiles(os.path.join(workPath,inputList+'_otzfif.list'))):
                        command = 'cp '+inputFile+' '+outputFile
                        print('command = <'+command+'>')
                        os.system(command)
                    for inputFile,outputFile in zip(getListOfFiles(os.path.join(workPath,inputList+'_otzxfi.list')),getListOfFiles(os.path.join(workPath,inputList+'_otzxfif.list'))):
                        command = 'cp '+inputFile+' '+outputFile
                        print('command = <'+command+'>')
                        os.system(command)
#                STOP
            for inputList in ['science','fluxstds']:
                if True:
                    subtractMedianSky(getListOfFiles(os.path.join(workPath,inputList+'_otzfif.list')))
                subtractMedianSky(getListOfFiles(os.path.join(workPath,inputList+'_otzxfif.list')))

        if False:
            # extract and reidentify ARCs
            # arc_otzxfif.list or arc_otzxf.list???
            arc_otzxf_list = getListOfFiles(os.path.join(workPath,'arc_otzxf.list'))
            wavelengthsOrig, wavelengthsResampled = extractAndReidentifyARCs(arc_otzxf_list,
                                                                            refProfApDef,
                                                                            lineList,
                                                                            referenceSpectrum,
                                                                            display=False,
                                                                            chiSquareLimit=0.45,
                                                                            degree=4,
                                                                            minLines=8,
                                                                            maxRMS=0.6,
                                                                            #apOffsetX=-155.5,
                                                                            )
            for i in range(len(arc_otzxf_list)):
                with open(arc_otzxf_list[i][:arc_otzxf_list[i].rfind('.')]+'i_wavelengthsOrig.dat','w') as f:
                    wLenOrig = wavelengthsOrig[i]
                    for wLen in wLenOrig:
                        f.write('%.5f\n' % (wLen))
            # extract arcs and science data
            #extractObjectAndSubtractSky(twoDImageFileIn, specOut, yRange, skyAbove, skyBelow, dispAxis)
            inputList = getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfif.list'))
            for inputFile in inputList:
                if not os.path.isfile(inputFile[:inputFile.rfind('.')]+'Ec.fits'):
                    extractSum(inputFile, 'row', fNameOut = inputFile[:inputFile.rfind('.')]+'Ec.fits')

            print('wavelengthsOrig = ',wavelengthsOrig)
            print('wavelengthsOrig[0] = ',wavelengthsOrig[0])
            inputList = getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfifEc.list'))
            for i in range(len(inputList)):
                with open(inputList[i][:inputList[i].rfind('.')]+'_wLenOrig.dat','w') as f:
                    for wLen in wavelengthsOrig[i]:
                        f.write('%.8f\n' % (wLen))
        #STOP
        if True:
            extractInputList = os.path.join(workPath,'to_extract.csv')
            areasFileName = os.path.join(workPath,'areas.csv')
            print('reduce: areasFileName = '+areasFileName)
            scienceSpectra = []#extractSum(fn,'row',fn[:-5]+'Ec.fits') for fn in getListOfFiles(os.path.join(workPath,'SCIENCE_otzfif.list'))]
            areasFileExists = False
            if os.path.exists(areasFileName):
                areas = csvFree.readCSVFile(areasFileName)
                print('areas.header = ',areas.header)
                areas.header = [key.replace('\r','') for key in areas.header]
                print('areas.header = ',areas.header)
                for i in range(areas.size()):
                    print('reduce: areas.getData("fName",',i,') = '+areas.getData('fName',i),', object = ',areas.getData('object',i),', skyBelow = ',areas.getData('skyBelow',i),', skyAbove = ',areas.getData('skyAbove',i,),', method = ',areas.getData('method',i),', notes = ',areas.getData('notes',i))
                extractList = areas
                areasFileExists = True
            else:
                extractList = csvFree.readCSVFile(extractInputList)
                with open(areasFileExists,'w') as fAreas:
                    fAreas.write('fName,object,skyBelow,skyAbove,method,notes\n')
        if True:
            for i in range(extractList.size()):
                if extractList.getData('fName',i)[0] != '#':
                    skyAbove = None
                    skyBelow = None
                    obsArea = None
                    extractionMethod = None
                    if areasFileExists:
                        obsArea = [int(areas.getData('object',i)[1:areas.getData('object',i).find(':')]),int(areas.getData('object',i)[areas.getData('object',i).find(':')+1:-1])]
                        skyAbove = None if areas.getData('skyAbove',i) == '' else [int(areas.getData('skyAbove',i)[1:areas.getData('skyAbove',i).find(':')]),int(areas.getData('skyAbove',i)[areas.getData('skyAbove',i).find(':')+1:-1])]
                        skyBelow = None if areas.getData('skyBelow',i) == '' else [int(areas.getData('skyBelow',i)[1:areas.getData('skyBelow',i).find(':')]),int(areas.getData('skyBelow',i)[areas.getData('skyBelow',i).find(':')+1:-1])]
                        extractionMethod = areas.getData('method',i)
                    print('reduce: i = ',i,': fName = ',extractList.getData('fName',i))
                    extractObjectAndSubtractSky(os.path.join(workPath,extractList.getData('fName',i)) if '/' not in extractList.getData('fName',i) else extractList.getData('fName',i),
                                                os.path.join(workPath,extractList.getData('fName',i)[:-5]+'-skyEc.fits'),
                                                obsArea,
                                                skyAbove = skyAbove,
                                                skyBelow = skyBelow,
                                                extractionMethod = extractionMethod,
                                                dispAxis = 'row',
                                                display = False,
                                                areasFileOut = areasFileName if not areasFileExists else None)
                    if False:
                        extractObjectAndSubtractSky(os.path.join(workPath,extractList.getData('fName',i)) if '/' not in extractList.getData('fName',i) else extractList.getData('fName',i),
                                                    os.path.join(workPath,extractList.getData('fName',i)[:-5]+'Ec.fits'),
                                                    obsArea,
                                                    skyAbove = None,
                                                    skyBelow = None,
                                                    extractionMethod = extractionMethod,
                                                    dispAxis = 'row',
                                                    display = False,
                                                    areasFileOut = areasFileName if not areasFileExists else None)
        if True:
            inputList = objectListsStartWith
            doHelioCor = True# if inputList == objectListsStartWith else False
            print('doHelioCor = ',doHelioCor)
            ecFiles = []
            secFiles = []
            skyFiles = []
            skyFilesEc = []
            ecdFiles = []
            secdFiles = []
            ecdfFiles = []
            secdfFiles = []
            skydFiles = []
            areasList = csvFree.readCSVFile(os.path.join(workPath,'areas.csv'))

            for fileName in areasList.getData('fName'):
                if fileName[0] != '#':
                    tmpFileName = fileName[:fileName.rfind('.')]+'-skyEc.fits'
                    if os.path.exists(tmpFileName):
                        secFiles.append(tmpFileName)
                    tmpFileName = fileName[:fileName.rfind('.')]+'Ec.fits'
                    if os.path.exists(tmpFileName):
                        ecFiles.append(tmpFileName)
                    tmpFileName = fileName[:fileName.rfind('.')].replace('-MedianSky','')+'MedianSky.fits'
                    if os.path.exists(tmpFileName):
                        skyFiles.append(tmpFileName)
                    tmpFileName = fileName[:fileName.rfind('.')].replace('-MedianSky','')+'MedianSkyEc.fits'
                    if os.path.exists(tmpFileName):
                        skyFilesEc.append(tmpFileName)
                    tmpFileName = fileName[:fileName.rfind('.')]+'Ecd.fits'
                    if os.path.exists(tmpFileName):
                        ecdFiles.append(tmpFileName)
                    else:
                        ecdFiles.append(fileName[:fileName.rfind('.')]+'Ecd.fits')
                    tmpFileName = fileName[:fileName.rfind('.')]+'-skyEcd.fits'
                    if os.path.exists(tmpFileName):
                        secdFiles.append(tmpFileName)
                    else:
                        secdFiles.append(fileName[:fileName.rfind('.')]+'-skyEcd.fits')
                    tmpFileName = fileName[:fileName.rfind('.')].replace('-MedianSky','')+'MedianSkyEcd.fits'
                    if os.path.exists(tmpFileName):
                        skydFiles.append(tmpFileName)
                    tmpFileName = fileName[:fileName.rfind('.')]+'EcdF.fits'
#                    if os.path.exists(tmpFileName):
                    ecdfFiles.append(tmpFileName)
#                    else:
#                        print('did not find ',tmpFileName)
#                        STOP
                    tmpFileName = fileName[:fileName.rfind('.')]+'-skyEcdF.fits'
                    #if os.path.exists(tmpFileName):
                    secdfFiles.append(tmpFileName)
                    # extract MedianSky images
                    if False:
                        skySpec = extractSum(skyFiles[len(skyFiles)-1],'row')

                        writeFits1D(skySpec,
                                    skyFilesEc[len(skyFilesEc)-1],
                                    wavelength=None,
                                    header=skyFiles[len(skyFiles)-1],
                                    CRVAL1=1,
                                    CRPIX1=1,
                                    CDELT1=1)
        if True:
            inputList = getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfifEc.list'))
            wavelengthsOrig = []
            for i in range(len(inputList)):
                wLenStr = readFileToArr(inputList[i][:inputList[i].rfind('.')]+'_wLenOrig.dat')
                wLens = [float(wLen) for wLen in wLenStr]
                wavelengthsOrig.append(np.asarray(wLens))
            print('wavelengthsOrig = ',wavelengthsOrig)
            print('wavelengthsOrig[0] = ',wavelengthsOrig[0])
        if True:
            print('len(skyFilesEc) = ',len(skyFilesEc))
            if len(skyFilesEc) > 0:
                dispCor(skyFilesEc,
                        getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfiEc.list')),
                        wavelengthsOrig,
                        skydFiles,
                        observatoryLocation,
                        'TELRA',#'RA',#
                        'TELDEC',#'DEC',#
                        'DATE-OBS',
                        doHelioCor = False)
            print('len(ecFiles) = ',len(ecFiles),', len(ecdFiles) = ',len(ecdFiles))
            if len(ecFiles) > 0:
                print('ecFiles = ',ecFiles)
                dispCor(ecFiles,
                        getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfiEc.list')),
                        wavelengthsOrig,
                        ecdFiles,
                        observatoryLocation,
                        'TELRA',#RA',#
                        'TELDEC',#DEC',#
                        'DATE-OBS',
                        doHelioCor = doHelioCor)
            print('len(secFiles) = ',len(secFiles),', len(secdFiles) = ',len(secdFiles))
            if len(secFiles) != len(secdFiles):
                print('secFiles = ',secFiles)
                print('secdFiles = ',secdFiles)
                STOP
            if len(secFiles) > 0:
                dispCor(secFiles,
                        getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfiEc.list')),
                        wavelengthsOrig,
                        secdFiles,
                        observatoryLocation,
                        'TELRA',#RA',#
                        'TELDEC',#DEC',#
                        'DATE-OBS',
                        doHelioCor = doHelioCor)

        if True:

            areas = csvFree.readCSVFile(os.path.join(workPath,'areas.csv'))
            sensFuncs = calcResponse(os.path.join(workPath,fluxstdsListsStartWith+'_otzxfif.list'),
                                    #getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfiEc.list')),
                                    #wavelengthsOrig,
                                    areas,
                                    stdStarNameEndsBefore='_a',
                                    display=False)
            print('sensFuncs = ',sensFuncs)
            if False:
                if len(ecdFiles) > 0:
                    applySensFuncs(ecdFiles,
                                ecdfFiles,
                                sensFuncs)
        if True:
            if len(secdFiles) > 0:
                applySensFuncs(secdFiles,
                            secdfFiles,
                            sensFuncs)

        # clean spectra
        if True:
            start = False
            if len(ecdfFiles) > 0:
                for ecdfFile in ecdfFiles:
                    if 'PHR0723+0036' in ecdfFile:
                        start = True
                    if start:
                        cleanSpec(ecdfFile,ecdfFile.replace('EcdF.fits','-sky.fits'),ecdfFile.replace('.fits','_clean.fits'))
            else:
                print('len(ecdfFiles) == 0')
            if len(secdfFiles) > 0:
                for secdfFile in secdfFiles:
                    cleanSpec(secdfFile,secdfFile.replace('-skyEcdF.fits','-sky.fits'),secdfFile.replace('.fits','_clean.fits'))
            else:
                print('len(secdfFiles) == 0')


        if False:
            fixWCSDim()
            # combine multiple spectra of the same object
            fileList = ecdfFiles
            toCombineLists = []
            combinedList = []
            toContinuumSubtract = []
            continuumSubtracted = []
            for fileName in fileList:
                objectName = fileName[fileName.rfind(objectListsStartWith+'_')+8:]
                objectName = objectName[:objectName.find('_a')]
                print('objectName = ',objectName)
                finalName = fileName[:fileName.rfind('/')+1]+objectName+'_SA'+fileName[fileName.rfind('/')-2:fileName.rfind('/')]+fileName[fileName.rfind('/')-4:fileName.rfind('/')-2]+fileName[fileName.rfind('/')-6:fileName.rfind('/')-4]+'.fits'
                print('finalName = <'+finalName+'>')

                nFiles = 0
                for fName in fileList:
                    objName = fName[fName.rfind(objectListsStartWith+'_')+8:]
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
                            objName = fName[fName.rfind(objectListsStartWith+'_')+8:]
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

        if False:
            # subtract continuum
            for i in range(len(toContinuumSubtract)):
                fittingFunction = np.polynomial.legendre.legfit
                evalFunction = np.polynomial.legendre.legval
                order = 5
                nIterReject = 2
                nIterFit = 3
                lowReject = 2.
                highReject = 1.8
                useMean = True
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
                        display=False)



        #    response = calcResponse(os.path.join(workPath,'FLUXSTDS_otzfifEcd.list'))
        #    fluxCalibrate(os.path.join(workPath,'SCIENCE_LTT7379_a1171120_otzfifEcd.fits'),
        #                  os.path.join(workPath,'SCIENCE_Feige110_a1171214_otzfifEcd.fits'))
