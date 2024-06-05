import numpy as np
import os
import csvFree,csvData
from drUtils import getHeaderValue,dmsToDeg,hmsToDeg
from myUtils import angularDistancePyAsl

path = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08'
pnMain = csvFree.readCSVFile(os.path.join(path,'hash_PNMain.csv'))
cNames = csvFree.readCSVFile(os.path.join(path,'hash_CNames.csv'))
fitsFiles = csvFree.readCSVFile(os.path.join(path,'hash_FitsFiles.csv'))
fitsFilesOrig = csvFree.readCSVFile(os.path.join(path,'hash_FitsFiles_orig.csv'))

pnMain_idPNMains = np.array(pnMain.getData('idPNMain'))

idPNMainsInFitsFiles = []

nFoundInCNames = 0
nFoundInPNMain = 0
with open(os.path.join(path,'make_changes_to_fitsFiles.sql'),'w') as f:
    for i in range(fitsFiles.size()):
        fName = fitsFiles.getData('fileName',i)
        objectName = fName[:fName.find('_MS')].lower().replace('_','')
        print('i = ',i,': objectName = ',objectName)
        setIDPNMain = fitsFiles.getData('idPNMain',i)
        fNameFull = os.path.join('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/REDUCED',fName)

        foundName = False
        for iNames in range(cNames.size()):
            cName = cNames.getData('Name',iNames).lower().replace(' ','')
            if (cName == objectName) or ('png'+objectName == cName):
                print('i = ',i,': cName = ',cName)
                foundName = True
                nFoundInCNames += 1
                idPNMainFromCNames = cNames.getData('idPNMain',iNames)
                if idPNMainFromCNames != setIDPNMain:
                    print('i = ',i,': fName = ',fName,': idPNMainFromCNames = ',idPNMainFromCNames,' != ',setIDPNMain,' = setIDPNMain')
                    fitsFiles.setData('idPNMain',i,idPNMainFromCNames)
                    if (not setIDPNMain in idPNMainsInFitsFiles) and (not setIDPNMain == '33649'):
                        idPNMainsInFitsFiles.append(setIDPNMain)
                    if (not idPNMainFromCNames in idPNMainsInFitsFiles) and (not idPNMainFromCNames == '33649'):
                        idPNMainsInFitsFiles.append(idPNMainFromCNames)
                fitsFiles.setData('object',i,cNames.getData('Name',iNames))
        if not foundName:
            for iPNMain in range(pnMain.size()):
                png = pnMain.getData('PNG',iPNMain)
                if png == objectName:
                    foundName = True
                    nFoundInPNMain += 1
                    idPNMainFromPNMain = pnMain.getData('idPNMain',iPNMain)
                    if idPNMainFromPNMain != setIDPNMain:
                        print('i = ',i,': fName = ',fName,': idPNMainFromPNMain = ',idPNMainFromPNMain,' != ',setIDPNMain,' = setIDPNMain')
                        STOP
                        fitsFiles.setData('idPNMain',i,idPNMainFromPNMain)
                        if not idPNMainFromPNMain in idPNMainsInFitsFiles:
                            idPNMainsInFitsFiles.append(idPNMainFromPNMain)

        if not foundName:
            print('Did not find name ',objectName,' in either list. fName = ',fName,', setIDPNMain = ',setIDPNMain)
            STOP
    #    if objectName == '033.7-02.0':
    #        STOP
        for j in range(len(fitsFiles.header)):
            if fitsFiles.getData(fitsFiles.header[j],i) != fitsFilesOrig.getData(fitsFiles.header[j],i):
                f.write("UPDATE `PNSpectra_Sources`.`FitsFiles` SET `"+fitsFiles.header[j]+"`=")
                if fitsFiles.header[j] != 'idPNMain':
                    f.write("'")
                f.write(fitsFiles.getData(fitsFiles.header[j],i))
                if fitsFiles.header[j] != 'idPNMain':
                    f.write("'")
                f.write(" WHERE `idFitsFiles`="+fitsFiles.getData('idFitsFiles',i)+";\n")

        RA = getHeaderValue(fNameFull,'RA')
        DEC = getHeaderValue(fNameFull,'DEC')
        idx = np.where(pnMain_idPNMains == fitsFiles.getData('idPNMain',i))[0]
        if idx < 0:
            print('idx < 0')
            STOP
        RA_HASH = pnMain.getData('RAJ2000',idx)
        DEC_HASH = pnMain.getData('DECJ2000',idx)
        dist = angularDistancePyAsl(hmsToDeg(RA),dmsToDeg(DEC),hmsToDeg(RA_HASH),dmsToDeg(DEC_HASH)) * 3600
        if dist > 60:
            print('\nRA / DEC = ',RA,DEC)
            print('RA_HASH / DEC_HASH = ',RA_HASH,DEC_HASH)
            print('dist = ',dist,'\n')
#            STOP

print('nFoundInCNames = ',nFoundInCNames,', nFoundInPNMain = ',nFoundInPNMain)

print('need to re-run ',len(idPNMainsInFitsFiles),': ',idPNMainsInFitsFiles)
with open(os.path.join(path,'idPNMainsInFitsFiles.sh'),'w') as f:
    f.write('hashpn spetch all '+idPNMainsInFitsFiles[0])
    for id in idPNMainsInFitsFiles[1:]:
        f.write(','+id)
