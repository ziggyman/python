import numpy as np
import os
from hashUtils import get_IDPNMain_from_name
from drUtils import readFileToArr,getHeaderValue,setHeaderValue,getHeader,merge
import csvFree,csvData

date = '2008-05-06'
path = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/interns/done'

objects_blue = readFileToArr(os.path.join(path,'objects_blue_2008-05-08.list'))
objects_red = readFileToArr(os.path.join(path,'objects_red_2008-05-08.list'))

def checkRMS(night, arm):
    import os
    from drUtils import getHeaderValue
    path = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/'+arm+'/'+night
    arcList = os.listdir(path)
    print('arcList = ',arcList)
    print('dir(arcList) = ',dir(arcList))
    arcList.sort()
    print('arcList = ',arcList)
    for file in arcList:
        if (file[:8] == 'ARC_CuHe') & file.endswith("Ecd.fits"):
            rms = float(getHeaderValue(os.path.join(path,file),'RMS'))
            print(os.path.join(path,file),' has RMS of ',rms)
            nLines = 0
            for i in range(20):
                if getHeaderValue(os.path.join(path,file),'LINE%d' % (i)) is not None:
                    nLines = i+1
            if (rms > 0.02) or (nLines < 12):
                print(os.path.join(path,file),' has RMS of ',rms,' which is too high or nLines = ',nLines,' which is too low')
                scienceList = os.listdir(path)
                scienceList.sort()
                for science in scienceList:
                    if (science[:7] == 'SCIENCE') and science.endswith("fEcd.fits"):
                        ref1 = getHeaderValue(os.path.join(path,science),'REFSPEC1')
                        ref1 = ref1[:ref1.find(' ')]
                        #print('ref1 = <'+ref1+'>')
                        if ref1.replace('Ec.','Ecd.') == file:
                            print('file ',file,' used for calibration of science ',science)
                        ref2 = getHeaderValue(os.path.join(path,science),'REFSPEC2')
                        ref2 = ref2[:ref2.find(' ')]
                        #print('ref2 = <'+ref2+'>')
                        if ref2.replace('Ec.','Ecd.') == file:
                            print('file ',file,' used for calibration of science ',science)
                        #print('science ',science,': ref1 = ',ref1,', ref2 = ',ref2)

def redoArc(date,arm,dbsNo):
    from drUtils import extractAndReidentifyARCs
    workPath=os.path.join('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/',arm,date)

    refApDef = os.path.join('/Users/azuri/stella/referenceFiles/dbs/database','aprefProfApDef_DBS_red_May2008_otzfif' if arm=='R_data' else 'aprefProfApDef_DBS_blue_May2008_otzfif_flipped')
    arcListIn = [os.path.join(workPath,'ARC_CuHe_dbs'+dbsNo+('b' if arm=='B_data' else 'r')+'_otzxf.fits')]
    lineListIn = os.path.join('/Users/azuri/stella/referenceFiles/dbs','dbs_refspec_may2008_red.txt' if arm=='R_data' else 'dbs_refspec_may2008_blue_CuHe.txt')
    xCorSpecIn = os.path.join('/Users/azuri/stella/referenceFiles/dbs','refArc_DBS_May2008_red_otzxfiEc.fits' if arm=='R_data' else 'refARC_DBS_May2008_blue_otzfiEc.fits')

    wavelengthsOrig, wavelengthsResampled = extractAndReidentifyARCs(arcListIn,
                                                                    refApDef,
                                                                    lineListIn,
                                                                    xCorSpecIn,
                                                                    display=True,
                                                                    chiSquareLimit=0.45,
                                                                    degree=5,
                                                                    minLines=13 if arm == 'B_data' else 10,
                                                                    maxRMS=0.3,
                                                                    shiftApDef=True)
    for i in range(len(arcListIn)):#arc_otzxf_list)):
        with open(arcListIn[i][:arcListIn[i].rfind('.')]+'i_wavelengthsOrig.dat','w') as f:
            wLenOrig = wavelengthsOrig[i]
            for wLen in wLenOrig:
                f.write('%.5f\n' % (wLen))
    inputList = [os.path.join(workPath,'ARC_CuHe_dbs'+dbsNo+('b' if arm=='B_data' else 'r')+'_otzxfifEc.fits')]
    for i in range(len(inputList)):
        with open(inputList[i][:inputList[i].rfind('.')]+'_wLenOrig.dat','w') as f:
            for wLen in wavelengthsOrig[i]:
                f.write('%.8f\n' % (wLen))


def redoDisp(date,arm,fName):
    from astropy.coordinates import EarthLocation
    from drUtils import dispCor,getListOfFiles,calcResponse,applySensFuncs
    workPath=os.path.join('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/',arm,date)
    scienceListIn = [os.path.join(workPath,fName)]
    scienceListOut = [os.path.join(workPath,fName[:fName.rfind('.')]+'d.fits')]

    arcListsStartWith = 'arc'
    fluxstdsListsStartWith = 'fluxstds'#DBS

    arcListIn = getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfiEc.list'))
    inputList = getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfifEc.list'))
    wavelengthsOrig = []
    for i in range(len(inputList)):
        wLenStr = readFileToArr(inputList[i][:inputList[i].rfind('.')]+'_wLenOrig.dat')
        wLens = [float(wLen) for wLen in wLenStr]
        wavelengthsOrig.append(np.asarray(wLens))



    observatoryLocation = EarthLocation.of_site('Mt. Stromlo Observatory')
    dispCor(scienceListIn,
            arcListIn,
            wavelengthsOrig,
            scienceListOut,
            observatoryLocation,
            'RA',#TELRA',
            'DEC',#TELDEC',
            'DATE-OBS',
            doHelioCor=True)

    areas = csvFree.readCSVFile(os.path.join(workPath,'areas.csv'))
    sensFuncs = calcResponse(os.path.join(workPath,fluxstdsListsStartWith+'_otzxfif.list'),
#                            getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfiEc.list')),
#                            wavelengthsOrig,
                            areas,
                            stdStarNameEndsBefore='dbs',
                            display=True)
    print('sensFuncs = ',sensFuncs)
    ecdFiles = [os.path.join(workPath,fName[:fName.rfind('.')]+'d.fits')]
    ecdfFiles = [os.path.join(workPath,fName[:fName.rfind('.')]+'dF.fits')]

    applySensFuncs(ecdFiles,
                ecdfFiles,
                sensFuncs)
    print('scienceListOut = ',scienceListOut)
    print('ecdfFiles = ',ecdfFiles)


def fixScienceDisp(fName):
    from drUtils import dispCor,getListOfFiles,extractAndReidentifyARCs
    workPath=os.path.join('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/',arm,date)
    inputList = getListOfFiles(os.path.join(workPath,arcListsStartWith+'_otzxfifEc.list'))
    wavelengthsOrig = []
    for i in range(len(inputList)):
        wLenStr = readFileToArr(inputList[i][:inputList[i].rfind('.')]+'_wLenOrig.dat')
        wLens = [float(wLen) for wLen in wLenStr]
        wavelengthsOrig.append(np.asarray(wLens))


def createCatalogue():
    fitsFiles = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_FitsFiles_020823.csv')
    with open('/Users/azuri/daten/uni/HKU/HASH/dbs_may2007_cat.sql','w') as f:
        f.write('CREATE TABLE IF NOT EXISTS MainPNData.DBS_May2008 (\n')
        f.write('idDBS_Mey2008 INT AUTO_INCREMENT PRIMARY KEY,\n')
        f.write('idPNMain INT,\n')
        f.write('mapFlag VARCHAR(1) NOT NULL\n')
        f.write(');\n')
        f.write("USE `MainPNData`;\n")

        id = 1
        for i in range(fitsFiles.size()):
            if fitsFiles.getData('setname',i) == 'DBS_May2008':
                f.write("INSERT INTO `DBS_May2008`(`idDBS_May2008`,`idPNMain`,`mapflag`)")
                f.write("VALUES (%d,%d,'%s');\n" % (id,int(fitsFiles.getData('idPNMain',i)),'y'))
                id += 1

if __name__ == '__main__':
    if len(objects_blue) != len(objects_red):
        print('different number of files in lists')
        STOP

    for i in range(len(objects_blue)):
        outFileName = os.path.join(path,getHeaderValue(objects_blue[i],'OBJECT')+'_MS080508.fits')
        if os.path.exists(outFileName):
            outFileName = outFileName[:outFileName.rfind('.')]+'a.fits'
            if os.path.exists(outFileName):
                outFileName = outFileName[:outFileName.rfind('a.')]+'b.fits'
                if os.path.exists(outFileName):
                    outFileName = outFileName[:outFileName.rfind('b.')]+'c.fits'
                    if os.path.exists(outFileName):
                        outFileName = outFileName[:outFileName.rfind('c.')]+'d.fits'
        merge(objects_blue[i],objects_red[i],outFileName)

    STOP

    areas_blue = csvFree.readCSVFile(os.path.join(path,'B_data',date,'areas.csv'))
    areas_red = csvFree.readCSVFile(os.path.join(path,'R_data',date,'areas.csv'))
    for i in range(areas_blue.size()):
        fNameBlue = areas_blue.getData('fName',i)
        fNameBlue = fNameBlue[fNameBlue.rfind('/')+1:]
        fNameBlue = fNameBlue[fNameBlue.find('_')+1:]
        fNameBlue = fNameBlue[:fNameBlue.find('_')]
        print('fNameBlue = <'+fNameBlue+'>')
        fNameRed = areas_red.getData('fName',i)
        fNameRed = fNameRed[fNameRed.rfind('/')+1:]
        fNameRed = fNameRed[fNameRed.find('_')+1:]
        fNameRed = fNameRed[:fNameRed.find('_')]
        print('fNameRed = <'+fNameRed+'>')
        if fNameRed != fNameBlue:
            print('line i=',i,': fNameRed = <'+fNameRed+'>, fNameBlue = <'+fNameBlue+'>')
            STOP
        for key in ['RA','DEC','AIRMASS']:
            if getHeaderValue(areas_blue.getData('fName',i),key) is None:#[:areas_blue.getData('fName',i).rfind('.')]+'Ec.fits',key) is None:
                setHeaderValue(areas_blue.getData('fName',i),key,getHeaderValue(areas_red.getData('fName',i)[:areas_red.getData('fName',i).rfind('.')]+'Ec.fits',key))
                if os.path.isfile(areas_blue.getData('fName',i)[:areas_blue.getData('fName',i).rfind('.')]+'Ec.fits'):
                    setHeaderValue(areas_blue.getData('fName',i)[:areas_blue.getData('fName',i).rfind('.')]+'Ec.fits',key,getHeaderValue(areas_red.getData('fName',i)[:areas_red.getData('fName',i).rfind('.')]+'Ec.fits',key))
                print('keyword ',key,' set for ',areas_blue.getData('fName',i))
    STOP

    data_path = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data'
    hash_path = '/Users/azuri/daten/uni/HKU/HASH'
    hash_names_filename = os.path.join(hash_path,'hash_tbCNames_300523.csv')
    hash_PNMain_filename = os.path.join(hash_path,'hash_PNMain_300523.csv')
    hash_PNMain = csvFree.readCSVFile(hash_PNMain_filename)
    fitsFilesOld = csvFree.readCSVFile(os.path.join(hash_path,'hash_FitsFiles_300523.csv'))
    fitsFilesNew = csvFree.readCSVFile(os.path.join(hash_path,'hash_FitsFiles_310523.csv'))

    fluxEcFileNames = readFileToArr(os.path.join(data_path,'fluxstds_otzxfifEc.list'))
    print('fluxEcFileNames = ',fluxEcFileNames)
    scienceEcFileNames = readFileToArr(os.path.join(data_path,'science_otzxfifEc.list'))
    print('scienceEcFileNames = ',scienceEcFileNames)
    objectNames = []
    for fName in scienceEcFileNames:
        if fName not in fluxEcFileNames:
            name = fName[fName.rfind('/SCIENCE')+9:fName.rfind('_dbs')]
            objectNames.append(name)
    idPNMains = get_IDPNMain_from_name(objectNames,hash_names_filename,hash_PNMain_filename)
    print('idPNMains = ',len(idPNMains),': ',idPNMains)
    idPNMainsStr = ''
    for i in range(len(idPNMains)):
        idPNMainsStr += idPNMains[i]+','
    print('name = ',name,': idPNMainStr = ',idPNMainsStr)

    notFounds = []
    for i in range(len(idPNMains)):
        if fitsFilesOld.find('idPNMain',idPNMains[i])[0] < 0:
            notFounds.append(idPNMains[i])
    print('notFounds = ',len(notFounds),': ',notFounds)

    for i in range(len(scienceEcFileNames)):
        scienceEcFileNames[i] = scienceEcFileNames[i][scienceEcFileNames[i].rfind('/')+1:].replace('Ec','EcdF')

    print('len(scienceEcFileNames) = ',len(scienceEcFileNames))
    print('len(idPNMains) = ',len(idPNMains))
    print('scienceEcFileNames = ',scienceEcFileNames)
    print('idPNMains = ',idPNMains)

    for j in range(len(scienceEcFileNames)):
        scienceFName = scienceEcFileNames[j]
        idPNMain = idPNMains[j]
        print('j = ',j,': scienceFName = ',scienceFName,', idPNMain = ',idPNMain)


    with open(os.path.join(hash_path,'sql_dbs_May2007.sql'),'w') as f:
        for i in np.arange(fitsFilesOld.size(),fitsFilesNew.size(),1):
            fName = fitsFilesNew.getData('fileName',i).replace('-MedianSky','')
            idFitsFiles = fitsFilesNew.getData('idFitsFiles',i)
            print('idFitsFiles = ',idFitsFiles,': fName = ',fName)
            found = False
            scienceFName = ''
            for j in range(len(scienceEcFileNames)):
                if fName == scienceEcFileNames[j]:
                    scienceFName = scienceEcFileNames[j]
                    idPNMain = idPNMains[j]
                    found = True
                    print('i = ',i,', j = ',j,': fName = ',fName,', scienceFName = ',scienceFName,', idPNMain = ',idPNMain,', idFitsFiles = ',idFitsFiles)
                    f.write("UPDATE `PNSpectra_Sources`.`FitsFiles` SET `idPNMain` = '"+idPNMain+"' WHERE (`idFitsFiles` = '"+idFitsFiles+"');\n")
                    f.write("UPDATE `PNSpectra_Sources`.`FitsFiles` SET `convToText` = 'y' WHERE (`idFitsFiles` = '"+idFitsFiles+"');\n")
            if not found:
                print('ERROR: fName <'+fName+'> not found')
                STOP

    lines = readFileToArr('/Users/azuri/entwicklung/python/data_reduction/dbs_run_may2007.out')
    with open('/Users/azuri/entwicklung/python/data_reduction/dbs_run_may2007.out2','w') as f:
        for line in lines:
            if ('readCSVFile' not in line) and (line[:8] != 'fileName'):
                f.write(line+'\n')

    STOP
    for fName in scienceEcFileNames:
        if fName not in fluxEcFileNames:
            ra = getHeaderValue(fName,'RA')
            print('ra = ',ra)
            if ra is None:
                setHeaderValue(fName,'RA',hash_PNMain.getData('RAJ2000',hash_PNMain.find('idPNMain',idPNMain)[0]))
            dec = getHeaderValue(fName,'DEC')
            print('dec = ',dec)
            if dec is None:
                setHeaderValue(fName,'DEC',hash_PNMain.getData('DECJ2000',hash_PNMain.find('idPNMain',idPNMain)[0]))
            airmass = getHeaderValue(fName,'AIRMASS')
            print('airmass = ',airmass)
            if airmass is None:
                setHeaderValue(fName,'AIRMASS','1.3')
