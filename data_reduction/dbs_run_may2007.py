import numpy as np
import os
from hashUtils import get_IDPNMain_from_name
from drUtils import readFileToArr,getHeaderValue,setHeaderValue
import csvFree,csvData

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
