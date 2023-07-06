import numpy as np
import os
import csvFree,csvData
#import pandas as pd

path = '/Users/azuri/daten/uni/HKU/HASH/CDS/'

pnMainFName = os.path.join(path,'hash_PNMain_26-06-2023.csv')
pnMain = csvFree.readCSVFile(pnMainFName)

removeKeys = ['refPNG',
              'refCatalogue',
              'userRecord',
              'refDomain',
              'refPNstat',
              'refSimbadID',
              'PNMaincol',
              'created_at',
              'updated_at']
for key in removeKeys:
    pnMain.removeColumn(key)

for i in np.arange(pnMain.size()-1,-1,-1):
    if pnMain.getData('show',i) == 'n':
        print('removing row ',i,' for idPNmain=',pnMain.getData('idPNMain',i),' due to no show')
        pnMain.removeRow(i)

AngDiamFName = os.path.join(path,'hash_AngDiam_26-06-2023.csv')
AngDiam = csvFree.readCSVFile(AngDiamFName)

idsAngDiam = np.array(AngDiam.getData('idPNMain'))
angDiamKeys = [['MajDiam','MajDiam'],
               ['errMajDiam','errMajDiam'],
               ['MinDiam','MinDiam'],
               ['errMinDiam','errMinDiam'],
               ['diamReference','reference'],
               ['diamRefTable','refTable'],
               ['diamRefRecord','refRecord']]
for key in angDiamKeys:
    pnMain.addColumn(key[0])

CNamesFName = os.path.join(path,'hash_CNames_26-06-2023.csv')
CNames = csvFree.readCSVFile(CNamesFName)

idsCNames = np.array(CNames.getData('idPNMain'))
CNamesKeys = [['Name','Name'],
              ['nameReference','reference'],
              ['nameRefTable','refTable'],
              ['nameRefRecord','refRecord']]
for key in CNamesKeys:
    pnMain.addColumn(key[0])

CSCoordsFName = os.path.join(path,'hash_CSCoords_26-06-2023.csv')
CSCoords = csvFree.readCSVFile(CSCoordsFName)

idsCSCoords = np.array(CSCoords.getData('idPNMain'))
CSCoordsKeys = [['CS_RAJ2000','CS_RAJ2000'],
                ['CS_DECJ2000','CS_DECJ2000'],
                ['CSType',None],
                ['CSReference','reference'],
                ['CSRefTable','refTable'],
                ['CSRefRecord','refRecord'],
                ['CSComments','comments']]
for key in CSCoordsKeys:
    pnMain.addColumn(key[0])

PNMorphFName = os.path.join(path,'hash_PNMorph_26-06-2023.csv')
PNMorph = csvFree.readCSVFile(PNMorphFName)

idsPNMorph = np.array(PNMorph.getData('idPNMain'))
PNMorphKeys = [['morphMainClass','mainClass'],
               ['morphSubClass','subClass'],
               ['morphReference','reference'],
               ['morphRefTable','refTable'],
               ['morphRefRecord','refRecord'],
               ['morphComments','comments']]
for key in PNMorphKeys:
    pnMain.addColumn(key[0])

PAFName = os.path.join(path,'hash_PA_26-06-2023.csv')
PA = csvFree.readCSVFile(PAFName)

idsPA = np.array(PA.getData('idPNMain'))
PAKeys = [['flagEPA','flagEPA'],
          ['EPA','EPA'],
          ['GPA','GPA'],
          ['errPA','errPA'],
          ['PAReference','reference'],
          ['PARefTable','refTable'],
          ['PARefRecord','refRecord'],
          ['PAComments','comments']]
for key in PAKeys:
    pnMain.addColumn(key[0])

HrvelFName = os.path.join(path,'hash_Hrvel_26-06-2023.csv')
Hrvel = csvFree.readCSVFile(HrvelFName)

idsHrvel = np.array(Hrvel.getData('idPNMain'))
HrvelKeys = [['H_rad_vel','H_rad_vel'],
             ['err_vel','vel_err'],
             ['velReference','reference'],
             ['velRefTable','refTable'],
             ['velRefRecord','refRecord']]
for key in HrvelKeys:
    pnMain.addColumn(key[0])

IRFluxFName = os.path.join(path,'hash_IRFlux_26-06-2023.csv')
IRFlux = csvFree.readCSVFile(IRFluxFName)

idsIRFlux = np.array(IRFlux.getData('idPNMain'))
IRFluxKeys = [['flagIRFlux','flagFlux'],
              ['IRFlux','Flux'],
              ['errIRFlux','errFlux'],
              ['IRFluxInstrument','instrument'],
              ['IRFluxBand','band'],
              ['IRFluxReference','reference'],
              ['IRFluxRefTable','refTable'],
              ['IRFluxRefRecord','refRecord']]
for key in IRFluxKeys:
    pnMain.addColumn(key[0])

IRMagFName = os.path.join(path,'hash_IRMag_26-06-2023.csv')
IRMag = csvFree.readCSVFile(IRMagFName)

idsIRMag = np.array(IRMag.getData('idPNMain'))
IRMagKeys = [['flagIRMag','flagMag'],
              ['IRMag','Mag'],
              ['errIRMag','errMag'],
              ['IRMagInstrument','instrument'],
              ['IRMagBand','band'],
              ['IRMagReference','reference'],
              ['IRMagRefTable','refTable'],
              ['IRMagRefRecord','refRecord']]
for key in IRMagKeys:
    pnMain.addColumn(key[0])

RadioContFName = os.path.join(path,'hash_RadioCont_26-06-2023.csv')
RadioCont = csvFree.readCSVFile(RadioContFName)

idsRadioCont = np.array(RadioCont.getData('idPNMain'))
RadioContKeys = [['flagRadioFlux','flagFlux'],
              ['RadioFlux','Flux'],
              ['errRadioFlux','errFlux'],
              ['RadioFluxFreq','freq'],
              ['RadioFluxBand','band'],
              ['RadioFluxReference','reference'],
              ['RadioFluxRefTable','refTable'],
              ['RadioFluxRefRecord','refRecord']]
for key in RadioContKeys:
    pnMain.addColumn(key[0])

logFHalphaFName = os.path.join(path,'hash_logFHalpha_26-06-2023.csv')
logFHalpha = csvFree.readCSVFile(logFHalphaFName)

idslogFHalpha = np.array(logFHalpha.getData('idPNMain'))
logFHalphaKeys = [['flaglogFHalpha','flaglogFlux'],
                  ['logFHalpha','logFlux'],
                  ['errlogFHalpha','errlogFlux'],
                  ['logFHalphaInstrument','instrument'],
                  ['logFHalphaReference','reference'],
                  ['logFHalphaRefTable','refTable'],
                  ['logFHalphaRefRecord','refRecord']]
for key in logFHalphaKeys:
    pnMain.addColumn(key[0])

logFOIIIFName = os.path.join(path,'hash_logFOIII_26-06-2023.csv')
logFOIII = csvFree.readCSVFile(logFOIIIFName)

idslogFOIII = np.array(logFOIII.getData('idPNMain'))
logFOIIIKeys = [['flaglogFOIII','flaglogFlux'],
                  ['logFOIII','logFlux'],
                  ['errlogFOIII','errlogFlux'],
                  ['logFOIIIInstrument','instrument'],
                  ['logFOIIIReference','reference'],
                  ['logFOIIIRefTable','refTable'],
                  ['logFOIIIRefRecord','refRecord']]
for key in logFOIIIKeys:
    pnMain.addColumn(key[0])

elcatFName = os.path.join(path,'hash_elcat_available.csv')
elcat = csvFree.readCSVFile(elcatFName)
idselcat = np.array(elcat.getData('idPNMain'))
elcatKeys = [['ELCAT','elcat']]
for key in elcatKeys:
    pnMain.addColumn(key[0])

literatureSpectraLinksFName = os.path.join(path,'hash_spectraLinks.csv')
literatureSpectraLinks = csvFree.readCSVFile(literatureSpectraLinksFName)
idsliteratureSpectraLinks = np.array(literatureSpectraLinks.getData('idPNMain'))
literatureSpectraLinksKeys = [['literatureSpectrum','reference']]
for key in literatureSpectraLinksKeys:
    pnMain.addColumn(key[0])

UsrCommFName = os.path.join(path,'hash_UsrComm_fixed.csv')
UsrComm = csvFree.readCSVFile(UsrCommFName)#,',',True)
idsUsrComm = np.array(UsrComm.getData('idPNMain'))
UsrCommKeys = [['comment','comment']]
for key in UsrCommKeys:
    pnMain.addColumn(key[0])

FitsFilesFName = os.path.join(path,'hash_FitsFiles.csv')
FitsFiles = csvFree.readCSVFile(FitsFilesFName)
idsFitsFiles = np.array(FitsFiles.getData('idPNMain'))
FitsFilesKeys = [['spectrum','parsedIn']]
for key in FitsFilesKeys:
    pnMain.addColumn(key[0])
removeKeysFitsFiles = ['outName',
                        'convToText',
                        'checked',
                        'inUse',
                        'parsedIn',
                        'parserflag',
                        'xfactor',
                        'tempname',
                        'created_at',
                        'updated_at']


for i in range(pnMain.size()):
    idPNMain = pnMain.getData('idPNMain',i)
    print('working on idPNMain = ',idPNMain)

    """Angular Extensions"""
    idxDiam = np.where(idsAngDiam == idPNMain)[0]
    print('idxDiam = ',idxDiam)
    if len(idxDiam) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in AngDiam')
    else:
        idxDiamInUse = -1
        for j in range(len(idxDiam)):
            if AngDiam.getData('InUse',idxDiam[j]) == '1':
                idxDiamInUse = idxDiam[j]
                print('idxDiamInUse = ',idxDiamInUse,': AngDiam[',idxDiamInUse,'] = ',AngDiam.getData(idxDiamInUse))
                for key in angDiamKeys:
                    if AngDiam.getData(key[1],idxDiamInUse) not in ['NULL','0']:
                        pnMain.setData(key[0],i,AngDiam.getData(key[1],idxDiamInUse))
        if idxDiamInUse == -1:
            print('WARNING: could not find an angular extension in use for idPNMain = ',idPNMain)

    """Common Names"""
    idxCNames = np.where(idsCNames == idPNMain)[0]
    print('idxCNames = ',idxCNames)
    if len(idxCNames) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in CNames')
    else:
        idxCNamesInUse = -1
        for j in range(len(idxCNames)):
            print('idxCNames[',j,'] = ',idxCNames[j])
            if CNames.getData('InUse',idxCNames[j]) == '1':
                idxCNamesInUse = idxCNames[j]
                print('idxCNamesInUse = ',idxCNamesInUse,': CNames[',idxCNamesInUse,'] = ',CNames.getData(idxCNamesInUse))
                for key in CNamesKeys:
                    if CNames.getData(key[1],idxCNamesInUse) not in ['NULL','0']:
                        pnMain.setData(key[0],i,CNames.getData(key[1],idxCNamesInUse))
        if idxCNamesInUse == -1:
            print('WARNING: could not find a Name in use for idPNMain = ',idPNMain)

    """CSCoords"""
    idxCSCoords = np.where(idsCSCoords == idPNMain)[0]
    print('idxCSCoords = ',idxCSCoords)
    if len(idxCSCoords) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in CSCoords')
    else:
        idxCSCoordsInUse = -1
        for j in range(len(idxCSCoords)):
            print('idxCSCoords[',j,'] = ',idxCSCoords[j])
            if CSCoords.getData('InUse',idxCSCoords[j]) == '1':
                idxCSCoordsInUse = idxCSCoords[j]
                print('idxCSCoordsInUse = ',idxCSCoordsInUse,': CSCoords[',idxCSCoordsInUse,'] = ',CSCoords.getData(idxCSCoordsInUse))
                for key in CSCoordsKeys:
                    if key[1] is not None:
                        if CSCoords.getData(key[1],idxCSCoordsInUse) not in ['NULL','0']:
                            pnMain.setData(key[0],i,CSCoords.getData(key[1],idxCSCoordsInUse))
        if idxCSCoordsInUse == -1:
            print('WARNING: could not find a CSCoords in use for idPNMain = ',idPNMain)

    """PNMorph"""
    idxPNMorph = np.where(idsPNMorph == idPNMain)[0]
    print('idxPNMorph = ',idxPNMorph)
    if len(idxPNMorph) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in PNMorph')
    else:
        idxPNMorphInUse = -1
        for j in range(len(idxPNMorph)):
            print('idxPNMorph[',j,'] = ',idxPNMorph[j])
            if PNMorph.getData('InUse',idxPNMorph[j]) == '1':
                idxPNMorphInUse = idxPNMorph[j]
                print('idxPNMorphInUse = ',idxPNMorphInUse,': PNMorph[',idxPNMorphInUse,'] = ',PNMorph.getData(idxPNMorphInUse))
                for key in PNMorphKeys:
                    if key[1] is not None:
                        if PNMorph.getData(key[1],idxPNMorphInUse) not in ['NULL','0']:
                            if key[1] == 'subClass':
                                pnMain.setData(key[0],i,PNMorph.getData(key[1],idxPNMorphInUse).replace(';',''))
                            else:
                                pnMain.setData(key[0],i,PNMorph.getData(key[1],idxPNMorphInUse))
        if idxPNMorphInUse == -1:
            print('WARNING: could not find a PNMorph in use for idPNMain = ',idPNMain)

    """PA"""
    idxPA = np.where(idsPA == idPNMain)[0]
    print('idxPA = ',idxPA)
    if len(idxPA) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in PA')
    else:
        idxPAInUse = -1
        for j in range(len(idxPA)):
            print('idxPA[',j,'] = ',idxPA[j])
            if PA.getData('InUse',idxPA[j]) == '1':
                idxPAInUse = idxPA[j]
                print('idxPAInUse = ',idxPAInUse,': PA[',idxPAInUse,'] = ',PA.getData(idxPAInUse))
                for key in PAKeys:
                    if key[1] is not None:
                        if PA.getData(key[1],idxPAInUse) not in ['NULL','0']:
                            pnMain.setData(key[0],i,PA.getData(key[1],idxPAInUse))
        if idxPAInUse == -1:
            print('WARNING: could not find a PA in use for idPNMain = ',idPNMain)

    """Hrvel"""
    idxHrvel = np.where(idsHrvel == idPNMain)[0]
    print('idxHrvel = ',idxHrvel)
    if len(idxHrvel) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in Hrvel')
    else:
        idxHrvelInUse = []
        for j in range(len(idxHrvel)):
            print('idxHrvel[',j,'] = ',idxHrvel[j])
            if Hrvel.getData('InUse',idxHrvel[j]) == '1':
                idxHrvelInUse.append(idxHrvel[j])
                print('idxHrvelInUse = ',idxHrvelInUse,': Hrvel[',idxHrvelInUse,'] = ',Hrvel.getData(idxHrvel[j]))
        if len(idxHrvelInUse) == 1:
            for key in HrvelKeys:
                if key[1] is not None:
                    if Hrvel.getData(key[1],idxHrvelInUse[0]) not in ['NULL','0']:
                        pnMain.setData(key[0],i,Hrvel.getData(key[1],idxHrvelInUse[0]))
        else:
            STOP
            vrad = 0.
            evrad = 0.
            nvrad = 0
            nevrad = 0
            for j in range(len(idxHrvelInUse)):
                thisvrad = Hrvel.getData('H_rad_vel',idxHrvelInUse[j])
#                if thisvrad
        if idxHrvelInUse == -1:
            print('WARNING: could not find a Hrvel in use for idPNMain = ',idPNMain)

    """IRFlux"""
    idxIRFlux = np.where(idsIRFlux == idPNMain)[0]
    print('idxIRFlux = ',idxIRFlux)
    if len(idxIRFlux) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in IRFlux')
    else:
        idxIRFluxInUse = -1
        for j in range(len(idxIRFlux)):
            print('idxIRFlux[',j,'] = ',idxIRFlux[j])
            if IRFlux.getData('InUse',idxIRFlux[j]) == '1':
                idxIRFluxInUse = idxIRFlux[j]
                print('idxIRFluxInUse = ',idxIRFluxInUse,': IRFlux[',idxIRFluxInUse,'] = ',IRFlux.getData(idxIRFluxInUse))
                for key in IRFluxKeys:
                    if key[1] is not None:
                        if IRFlux.getData(key[1],idxIRFluxInUse) not in ['NULL','0']:
                            pnMain.setData(key[0],i,IRFlux.getData(key[1],idxIRFluxInUse))
        if idxIRFluxInUse == -1:
            print('WARNING: could not find a IRFlux in use for idPNMain = ',idPNMain)

    """IRMag"""
    idxIRMag = np.where(idsIRMag == idPNMain)[0]
    print('idxIRMag = ',idxIRMag)
    if len(idxIRMag) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in IRMag')
    else:
        idxIRMagInUse = -1
        for j in range(len(idxIRMag)):
            print('idxIRMag[',j,'] = ',idxIRMag[j])
            if IRMag.getData('InUse',idxIRMag[j]) == '1':
                idxIRMagInUse = idxIRMag[j]
                print('idxIRMagInUse = ',idxIRMagInUse,': IRMag[',idxIRMagInUse,'] = ',IRMag.getData(idxIRMagInUse))
                for key in IRMagKeys:
                    if key[1] is not None:
                        if IRMag.getData(key[1],idxIRMagInUse) not in ['NULL','0']:
                            pnMain.setData(key[0],i,IRMag.getData(key[1],idxIRMagInUse))
        if idxIRMagInUse == -1:
            print('WARNING: could not find a IRMag in use for idPNMain = ',idPNMain)

    """RadioFlux"""
    idxRadioFlux = np.where(idsRadioCont == idPNMain)[0]
    print('idxRadioFlux = ',idxRadioFlux)
    if len(idxRadioFlux) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in RadioFlux')
    else:
        idxRadioFluxInUse = -1
        for j in range(len(idxRadioFlux)):
            print('idxRadioFlux[',j,'] = ',idxRadioFlux[j])
            if RadioCont.getData('InUse',idxRadioFlux[j]) == '1':
                idxRadioFluxInUse = idxRadioFlux[j]
                print('idxRadioFluxInUse = ',idxRadioFluxInUse,': RadioCont[',idxRadioFluxInUse,'] = ',RadioCont.getData(idxRadioFluxInUse))
                for key in RadioContKeys:
                    if key[1] is not None:
                        if RadioCont.getData(key[1],idxRadioFluxInUse) not in ['NULL','0']:
                            pnMain.setData(key[0],i,RadioCont.getData(key[1],idxRadioFluxInUse))
        if idxRadioFluxInUse == -1:
            print('WARNING: could not find a RadioFlux in use for idPNMain = ',idPNMain)

    """logFHalpha"""
    idxlogFHalpha = np.where(idslogFHalpha == idPNMain)[0]
    print('idxlogFHalpha = ',idxlogFHalpha)
    if len(idxlogFHalpha) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in logFHalpha')
    else:
        idxlogFHalphaInUse = -1
        for j in range(len(idxlogFHalpha)):
            print('idxlogFHalpha[',j,'] = ',idxlogFHalpha[j])
            if logFHalpha.getData('InUse',idxlogFHalpha[j]) == '1':
                idxlogFHalphaInUse = idxlogFHalpha[j]
                print('idxlogFHalphaInUse = ',idxlogFHalphaInUse,': logFHalpha[',idxlogFHalphaInUse,'] = ',logFHalpha.getData(idxlogFHalphaInUse))
                for key in logFHalphaKeys:
                    if key[1] is not None:
                        if logFHalpha.getData(key[1],idxlogFHalphaInUse) not in ['NULL','0']:
                            pnMain.setData(key[0],i,logFHalpha.getData(key[1],idxlogFHalphaInUse))
        if idxlogFHalphaInUse == -1:
            print('WARNING: could not find a logFHalpha in use for idPNMain = ',idPNMain)

    """logFOIII"""
    idxlogFOIII = np.where(idslogFOIII == idPNMain)[0]
    print('idxlogFOIII = ',idxlogFOIII)
    if len(idxlogFOIII) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in logFOIII')
    else:
        idxlogFOIIIInUse = -1
        for j in range(len(idxlogFOIII)):
            print('idxlogFOIII[',j,'] = ',idxlogFOIII[j])
            if logFOIII.getData('InUse',idxlogFOIII[j]) == '1':
                idxlogFOIIIInUse = idxlogFOIII[j]
                print('idxlogFOIIIInUse = ',idxlogFOIIIInUse,': logFOIII[',idxlogFOIIIInUse,'] = ',logFOIII.getData(idxlogFOIIIInUse))
                for key in logFOIIIKeys:
                    if key[1] is not None:
                        if logFOIII.getData(key[1],idxlogFOIIIInUse) not in ['NULL','0']:
                            pnMain.setData(key[0],i,logFOIII.getData(key[1],idxlogFOIIIInUse))
        if idxlogFOIIIInUse == -1:
            print('WARNING: could not find a logFOIII in use for idPNMain = ',idPNMain)

    """ELCAT"""
    idxelcat = np.where(idselcat == idPNMain)[0]
    print('idxelcat = ',idxelcat)
    if len(idxelcat) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in elcat')
    elif len(idxelcat) == 1:
        print('idxelcat = ',idxelcat,': elcat[',idxelcat,'] = ',elcat.getData(idxelcat[0]))
        for key in elcatKeys:
            pnMain.setData(key[0],i,elcat.getData(key[1],idxelcat[0]))
    else:
        print('more than 1 elcat entry found for idPNMain = ',idPNMain)
        STOP

    """literature spectra links"""
    idxliteratureSpectraLinks = np.where(idsliteratureSpectraLinks == idPNMain)[0]
    print('idxliteratureSpectraLinks = ',idxliteratureSpectraLinks)
    if len(idxliteratureSpectraLinks) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in elcat')
    elif len(idxelcat) == 1:
        print('idxelcat = ',idxelcat,': elcat[',idxelcat,'] = ',elcat.getData(idxelcat[0]))
        for key in elcatKeys:
            pnMain.setData(key[0],i,elcat.getData(key[1],idxelcat[0]))
    else:
        print('more than 1 literature spectrum links entry found for idPNMain = ',idPNMain)
        STOP

    """UsrComm"""
    idxUsrComm = np.where(idsUsrComm == idPNMain)[0]
    print('idxUsrComm = ',idxUsrComm)
    if len(idxUsrComm) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in UsrComm')
    elif len(idxUsrComm) == 1:
        print('idxUsrComm = ',idxUsrComm,': UsrComm[',idxUsrComm,'] = ',UsrComm.getData(idxUsrComm[0]))
        for key in UsrCommKeys:
            pnMain.setData(key[0],i,UsrComm.getData(key[1],idxUsrComm[0]))
    else:
        comment = ''
        for j in range(len(idxUsrComm)):
            if len(comment) > 0:
                comment += '; '
            comment += UsrComm.getData(UsrCommKeys[0][1],idxUsrComm[j])
        pnMain.setData(UsrCommKeys[0][0],i,comment)
        print('more than 1 UsrComm entry found for idPNMain = ',idPNMain,': comment = <'+comment+'>')

    """FitsFiles"""
    idxFitsFiles = np.where(idsFitsFiles == idPNMain)[0]
    print('idxFitsFiles = ',idxFitsFiles)
    if len(idxFitsFiles) == 0:
        print('WARNING: did not find idPNMain ',idPNMain,' in FitsFiles')
        pnMain.setData(FitsFilesKeys[0][0],i,'n')
    else:
        hasSpec = False
        for j in range(len(idxFitsFiles)):
            if (FitsFiles.getData('parsedIn',idxFitsFiles[j]) == 'y') and (FitsFiles.getData('inUse',idxFitsFiles[j]) == '1'):
                hasSpec = True
        if hasSpec:
            pnMain.setData(FitsFilesKeys[0][0],i,'y')
        print('more than 1 UsrComm entry found for idPNMain = ',idPNMain,': hasSpec = ',hasSpec)



for key in removeKeysFitsFiles:
    FitsFiles.removeColumn(key)
csvFree.writeCSVFile(FitsFiles,os.path.join(path,'hash_FitsFiles_CDS.csv'))

pnMain.addColumn('HASH_link')
for i in range(pnMain.size()):
    pnMain.setData('HASH_link',i,'http://202.189.117.101:8999/gpne/objectInfoPage.php?id='+pnMain.getData('idPNMain',i))

csvFree.writeCSVFile(pnMain,os.path.join(path,'hash_PNMain_CDS.csv'))

pnMain.setData('CSType',0,'    ')
pnMain.setData('IRMagReference',0,'                   ')
pnMain.setData('flaglogFOIII',0,'     ')
pnMain.setData('logFOIII',0,'        ')
pnMain.setData('errlogFOIII',0,'         ')
pnMain.setData('logFOIIIInstrument',0,'      ')
pnMain.setData('logFOIIIReference',0,'                   ')
pnMain.setData('logFOIIIRefTable',0,'                ')
pnMain.setData('logFOIIIRefRecord',0,'    ')
pnMain.setData('literatureSpectrum',0,'                   ')

for key in pnMain.header:
    col = pnMain.getData(key)
    lens = [len(x) for x in col]
    iMax = np.where(np.array(lens) == max(lens))
    print('key = <'+key+'>: max characters = ',max(lens))
    print('key = <'+key+'>: example entry = "'+pnMain.getData(key,iMax[0][0])+'"')

print('pnMain header = ')
for key in pnMain.header:
    print(key)

mainTable = csvData.CSVData()
mainTable.header = ['idPNMain']
mainTable.data = [[i] for i in pnMain.getData('idPNMain')]
keys = [['idPNMain','I5','---','Unique HASH identifier'],
        ['PNG','A11','---','Name of PN, based on Galactic position (1)'],
        ['RAJ2000','','HH:MM:SS.S','Right Ascension J2000'],
        ['DECJ2000','','DD:MM:SS','Declination J2000'],
        ['Glon','F8.4','deg','Galactic longitude'],
        ['Glat','F8.4','deg','Galactic latitude'],
        ['refCoord','A19','---','Coordinate reference'],
        ['Catalogue','A20','---','Catalogue name'],
        ['Name','A31','---','Usual name of the PN (2)'],
        ['nameReference','A20','---','Name reference'],
        ['nameRefTable','A22','---','Name reference table'],
        ['nameRefRecord','I4','---','Name reference table record number'],
        ['SimbadID','A34','---','Simbad ID'],
        ['PNstat','A4','---','PN status (3)'],
        ['domain','A6','---','Domain: Galaxy, LMC, SMC'],
        ['morphMainClass','A1','---','Morphological Main Class (4)'],
        ['morphSubClass','A5','---','Morphologial Sub Class (5)'],
        ['morphReference','A19','---','Morphology reference'],
        ['morphRefTable','A22','---','Morphology reference table'],
        ['morphRefRecord','I4','---','Morphology reference table record number'],
        ['morphComments','A96','---','Comments on Morphology, e.g. from which image'],
        ['MajDiam','I5','arcsec','Major Diameter'],
        ['errMajDiam','I3','arcsec','Major Diameter uncertainty'],
        ['MinDiam','I5','arcsec','Minor Diameter'],
        ['errMinDiam','I3','arcsec','Minor Diameter uncertainty'],
        ['diamReference','A19','---','Diameter reference'],
        ['diamRefTable','A21','---','Diameter reference table'],
        ['diamRefRecord','I4','---','Diameter reference table record number'],
        ['flagEPA','I1','---','EPA measurement quality flag (6)'],
        ['EPA','I3','deg','Measured Equatorial Position Angle'],
        ['GPA','I3','deg','Calculated Galactic Position Angle'],
        ['errPA','F4.1','deg','EPA uncertainty'],
        ['PAReference','A19','---','Position Angle reference'],
        ['PARefTable','A18','---','Position Angle reference table'],
        ['PARefRecord','I4','---','Position Angle reference table record number'],
        ['PAComments','A25','---','Position Angle comments, e.g. from which image'],
        ['H_rad_vel','F6.1','km/s','Heliocentric Radial Velocity'],
        ['err_vel','F5.1','km/s','Heliocentric Radial Velocity uncertainty'],
        ['velReference','A19','---','Velocity reference'],
        ['velRefTable','A16','---','Velocity reference table'],
        ['velRefRecord','I4','---','Velocity reference table record number'],
        ['CS_RAJ2000','A','HH:MM:SS.S','Central Star Right Ascension J2000'],
        ['CS_DECJ2000','A','DD:MM:SS','Central Star Declination J2000'],
        ['CSType','A4','---','Central Star type, e.g. blue star, wels,...'],
        ['CSReference','A19','---','Central Star reference'],
        ['CSRefTable','A21','---','Central Star reference table'],
        ['CSRefRecord','I4','---','Central Star reference table record number'],
        ['CSComments','A62','---','Central Star comments, e.g. from which image'],
        ['flagIRFlux','A4','---','IR Flux flag'],
        ['IRFlux','F12.9','',''],
        ['errIRFlux','F12.9','',''],
        ['IRFluxInstrument','A15','',''],
        ['IRFluxBand','A5','',''],
        ['IRFluxReference','A19','',''],
        ['IRFluxRefTable','A19','',''],
        ['IRFluxRefRecord','I4','',''],
        ['flagIRMag','A4','',''],
        ['IRMag','F6.3','',''],
        ['errIRMag','F5.3','',''],
        ['IRMagInstrument','A3','',''],
        ['IRMagBand','A2','',''],
        ['IRMagReference','A19','',''],
        ['IRMagRefTable','A21','',''],
        ['IRMagRefRecord','I4','',''],
        ['flagRadioFlux','A3','',''],
        ['RadioFlux','F7.1','',''],
        ['errRadioFlux','F7.2','',''],
        ['RadioFluxFreq','F5.3','',''],
        ['RadioFluxBand','A4','',''],
        ['RadioFluxReference','A19','',''],
        ['RadioFluxRefTable','A24','',''],
        ['RadioFluxRefRecord','I6','',''],
        ['flaglogFHalpha','A5','',''],
        ['logFHalpha','F8.4','',''],
        ['errlogFHalpha','F9.6','',''],
        ['logFHalphaInstrument','A6','',''],
        ['logFHalphaReference','A39','',''],
        ['logFHalphaRefTable','A16','',''],
        ['logFHalphaRefRecord','I4','',''],
        ['flaglogFOIII','A5','',''],
        ['logFOIII','F8.4','',''],
        ['errlogFOIII','F9.6','',''],
        ['logFOIIIInstrument','A6','',''],
        ['logFOIIIReference','A19','',''],
        ['logFOIIIRefTable','A20','',''],
        ['logFOIIIRefRecord','I4','',''],
        ['ELCAT','A1','',''],
        ['literatureSpectrum','A19','',''],
        ['comment','A1120','',''],
        ['spectrum','A1','',''],
        ['HASH_link','A60','',''],
        ]
for key in keys[1:]:
    mainTable.addColumn(key[0],pnMain.getData(key[0]))



if False:
    with open(os.path.join(path,'pnMainDescription.txt'),'w') as f:
    #    f.write('  17- 29  A13   ---     Name     Usual name of the PN (2)\n')
        f.write('   1-  5  I5    ---     idPNMain     HASH identifier\n')
        f.write('   7- 18  A11   ---     PNG          Name of PN, based on Galactic position (1)\n')
        f.write('  20- 21  I2    h       RAh          Right Ascension J2000 (hours)\n')
        f.write('  23- 24  I2    min     RAm      Right Ascension J2000 (minutes)\n')
        f.write('  26- 29  F4.1  s       RAs      Right Ascension J2000 (seconds)\n')
        f.write('      31  A1    ---     DE-      Declination J2000 (sign)\n')
        f.write('  32- 33  I2    deg     DEd      Declination J2000 (degrees)\n')
        f.write('  35- 36  I2    arcmin  DEm      Declination J2000 (minutes)\n')
        f.write('  38- 39  I2    arcsec  DEs      Declination J2000 (seconds)\n')
        f.write('  41- 48  F8.4  deg     GLon     Galactic longitude\n')
        f.write('  50- 57  F8.4  deg     GLat     Galactic latitude\n')

        f.write('  73- 78  F6.1  arcsec  MajDiam  ? Major diameter of the PN from H{alpha} image\n')
        f.write('  80- 85  F6.1  arcsec  MinDiam  ? Minor diameter of the PN from H{alpha} image  A11   ---     PNG          Name of PN, based on Galactic position (1)\n')


        f.write('Note (1): This designation follows the accepted IAU convention introduced\n')
        f.write('in the "Strasbourg-ESO Catalogue of Galactic Planetary Nebulae" (Acker\n')
        f.write('et al., 1992, Cat. V/84)\n')

maxLenKeys = max([len(k) for k in keys[:][0]])
print('maxLenKeys = ',maxLenKeys)

pnMaintxt = ['%5s' % x for x in mainTable.getData('idPNMain')]
#print('pnMaintxt = ',pnMaintxt)

with open(os.path.join(path,'pnMainDescription.txt'),'w') as f:
    nChar = 1
    for key in keys:
        col = mainTable.getData(key[0])
        if key[0] == 'RAJ2000':
            strOut = '%4i-%4i  I2    h        RAh                   Right Ascension J2000 (hours)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3

            strOut = '%4i-%4i  I2    min      RAm                   Right Ascension J2000 (minutes)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3

            strOut = '%4i-%4i  F4.1  s        RAs                   Right Ascension J2000 (seconds)\n' % (nChar,nChar+3)
            f.write(strOut)
            nChar += 5
            for i in range(len(pnMaintxt)):
                pnMaintxt[i] = pnMaintxt[i]+' '+col[i][:col[i].find('.')+2]
        elif key[0] == 'DECJ2000':
            strOut = '%9i  A1    ---      DE-                   Declination J2000 (sign)\n' % (nChar)
            f.write(strOut)
            nChar += 1

            strOut = '%4i-%4i  I2    deg      DEd                   Declination J2000 (degrees)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3

            strOut = '%4i-%4i  I2    arcmin   DEm                   Declination J2000 (minutes)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3

            strOut = '%4i-%4i  I2    arcsec   DEs                   Declination J2000 (seconds)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3
            for i in range(len(pnMaintxt)):
                pnMaintxt[i] = pnMaintxt[i]+' '+col[i][:col[i].find('.')]
        elif key[0] == 'CS_RAJ2000':
            strOut = '%4i-%4i  I2    h        RAh                   Central Star: Right Ascension J2000 (hours)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3

            strOut = '%4i-%4i  I2    min      RAm                   Central Star: Right Ascension J2000 (minutes)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3

            strOut = '%4i-%4i  F4.1  s        RAs                   Central Star: Right Ascension J2000 (seconds)\n' % (nChar,nChar+3)
            f.write(strOut)
            nChar += 5
            for i in range(len(pnMaintxt)):
                if col[i] == '':
                    pnMaintxt[i] = pnMaintxt[i]+'           '
                else:
                    pnMaintxt[i] = pnMaintxt[i]+' '+col[i][:col[i].find('.')+2]
        elif key[0] == 'CS_DECJ2000':
            strOut = '%9i  A1    ---      DE-                   Central Star: Declination J2000 (sign)\n' % (nChar)
            f.write(strOut)
            nChar += 1

            strOut = '%4i-%4i  I2    deg      DEd                   Central Star: Declination J2000 (degrees)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3

            strOut = '%4i-%4i  I2    arcmin   DEm                   Central Star: Declination J2000 (minutes)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3

            strOut = '%4i-%4i  I2    arcsec   DEs                   Central Star: Declination J2000 (seconds)\n' % (nChar,nChar+1)
            f.write(strOut)
            nChar += 3
            for i in range(len(pnMaintxt)):
                if col[i] == '':
                    pnMaintxt[i] = pnMaintxt[i]+'          '
                else:
                    pnMaintxt[i] = pnMaintxt[i]+' '+col[i][:col[i].find('.')]
        else:
            col = pnMain.getData(key[0])
            length = key[1][1:]
            if '.' in length:
                length = int(length[:length.find('.')])
            else:
                length = int(length)
            strOut = '%4i-%4i  %s' % (nChar,nChar+length-1,key[1])
            while len(strOut) < 17:
                strOut += ' '
            strOut += '%s' % (key[2])
            while len(strOut) < 25 + 1:
                strOut += ' '
            strOut += key[0]
            while len(strOut) < 25+maxLenKeys+1:
                strOut += ' '
            strOut += key[3]
            strOut += '\n'
            nChar += length+1
            print('key = ',key,': nChar = ',nChar)
            f.write(strOut)
            for i in range(len(pnMaintxt)):
                if col[i].replace(' ','') in ['NULL','']:
                    oldLength = len(pnMaintxt[i])
                    while len(pnMaintxt[i]) < oldLength+length+1:
                        pnMaintxt[i] = pnMaintxt[i]+' '
                else:
                    print('col[',i,'] = ',col[i])
                    str = pnMaintxt[i].replace('%','%%')+' %'
                    if key[1][0] == 'A':
                        char = 's'
                        val = col[i].replace('%','%%')
                    elif key[1][0] == 'I':
                        char = 'i'
                        val = round(float(col[i]))
                    elif key[1][0] == 'F':
                        char = 'f'
                        val = float(col[i])
                    str += key[1][1:]+char
                    print('str = <'+str+'>, val = ',type(val),': ',val)
                    pnMaintxt[i] = str % (val)

    f.write('Note (1): This designation follows the accepted IAU convention introduced\n')
    f.write('    in the "Strasbourg-ESO Catalogue of Galactic Planetary Nebulae" (Acker\n')
    f.write('    et al., 1992, Cat. V/84)\n')
    f.write('Note (2): This name usually uses initials of prime discoverer surnames\n')
    f.write('    followed by truncated equatorial coordinates; the initials are PHR\n')
    f.write('    (Parker, Hartley, Russeil), PPA (Peyaud, Parker, Acker) and FP (Frew,\n')
    f.write('    Parker). The French Amateurs usually use the first 2 letters of the\n')
    f.write('    discoverer surnames and a running number.')
    f.write('Note (3): T: True PN\n')
    f.write('          L: Likely PN\n')
    f.write('          P: Possible PN\n')
    f.write('          agb: AGB Star\n')
    f.write('          agb?: AGB Star candidate\n')
    f.write('          q: artefact\n')
    f.write('          Be*: Be Star\n')
    f.write('          CV*: Cataclysmic Variable Star\n')
    f.write('          cir: CircumStellar Matter\n')
    f.write('          Cl*: Cluster of Stars\n')
    f.write('          CGb: Cometary Globule\n')
    f.write('          Em*: Emission Line Star\n')
    f.write('          EmO: Emission Object\n')
    f.write('          G: Galaxy\n')
    f.write('          HH: Herbig-Haro Object\n')
    f.write('          HII: HII Region\n')
    f.write('          o: Interesting Object\n')
    f.write('          iISM: Ionized ISM\n')
    f.write('          c: New Candidates\n')
    f.write('          X: Not PN\n')
    f.write('          N: Not PN (check)\n')
    f.write('          ooun: Object of unknown nature\n')
    f.write('          s: Objects to check\n')
    f.write('          a: PAGB/Pre-PN\n')
    f.write('          Be?: Possible Be Star\n')
    f.write('          ?Em*: Possible Emission Line Star\n')
    f.write('          G?: Possible Galaxy\n')
    f.write('          HH?: Possible Herbig-Haro Object\n')
    f.write('          b: Possible PrePN\n')
    f.write('          ev?: Possible Transient Event\n')
    f.write('          e: RCrB/eHe/LTP\n')
    f.write('          RNe: Reflection Nebula\n')
    f.write('          r: RV Tau\n')
    f.write('          star: Star\n')
    f.write('          SNR: SuperNova Remnant\n')
    f.write('          SR?: Supernova Remnant Candidate\n')
    f.write('          Sy*: Symbiotic Star\n')
    f.write('          Sy?: Symbiotic Star Candidate\n')
    f.write('          test: Test Object\n')
    f.write('          ev: Transient Event\n')
    f.write('          k: Transition Object\n')
    f.write('          WD*: White Dwarf / Hot Subdwarf\n')
    f.write('          Y*O: Young Stellar Object\n')
    f.write('          Y*?: Young Stellar Object Candidate\n')
    f.write('Note (4): E: Elliptical/oval\n')
    f.write('          R: Round\n')
    f.write('          B: Bipolar\n')
    f.write('          I: Irregular\n')
    f.write('          A: Asymmetric\n')
    f.write('          S: quasi-Stellar\n')
    f.write('Note (5): a: one sided enhancement/asymmetry\n')
    f.write('          m: multiple shells/external structure\n')
    f.write('          p: point symmetry\n')
    f.write('          r: ring structure/annulus\n')
    f.write('          s: internal structure\n')
    f.write('Note (6): 1: very good image\n')
    f.write('          2: pretty good image\n')
    f.write('          3: poor image\n')
