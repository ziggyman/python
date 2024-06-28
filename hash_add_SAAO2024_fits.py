import os
import numpy as np
import csvFree,csvData
from drUtils import getHeaderValue
from myUtils import hmsToDeg,dmsToDeg

fitsFiles = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/observing/hash_FitsFiles_050624.csv')
CNames = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/observing/hash_tbCNames_010624.csv')

fitsDir = '/Users/azuri/spectra/SAAO_June2024/110624/hash/'
fileListName = os.path.join(fitsDir,'files.list')

sqlFile = os.path.join(fitsDir,'addDataToFitsFiles.sql')
sqlCatFile = os.path.join(fitsDir,'addDataToCatalogue.sql')

command = 'ls '+fitsDir+'*.fits > '+fileListName
os.system(command)

with open(fileListName,'r') as f:
    fileList = f.readlines()
fileList = [f.strip() for f in fileList]

fileNamesInFitsFiles = np.array(fitsFiles.getData('fileName'))
objNames = CNames.getData('Name')
objNames = np.array([name.lower().replace(' ','') for name in objNames])

idPNMains = []
with open(sqlFile,'w') as f:
    f.write('USE PNSpectra_Sources;\n')
    for fileName in fileList:
        fName = fileName[fileName.rfind('/')+1:]
        objName = fName[:fName.find('_')]
        if objName == 'Mz3-CS':
            objName = 'Mz3'
        elif objName == 'DS1-CS':
            objName = 'DS1'
        elif objName == 'NGC7094-CS':
            objName = 'NGC7094'
        elif objName == 'YP0821-4253-right':
            objName = 'YP0821-4253'
        elif 'PHRJ1900-0014' in objName:
            objName = 'PHRJ1900-0014'
        elif 'YP0821-4253' in objName:
            objName = 'YP0821-4253b'
        elif objName == 'V866ScoA':
            objName = 'V*V866ScoA'
        elif 'Sh2-71' in objName:
            objName = 'Sh2-71'
        elif 'YP0819-2928' in objName:
            objName = 'YP0819-2928'
        elif 'YP1731-3011' in objName:
            objName = 'YP1731-3011'

        idx = np.where(objNames == objName.lower())[0]
        print('objName = ',objName,': idx = ',idx)
        idPNMain = CNames.getData('idPNMain',idx)
        idPNMains.append(idPNMain)
        print("getHeaderValue(fileName,'TELRA') = <"+getHeaderValue(fileName,'TELRA')+'>')
        print("getHeaderValue(fileName,'TELDEC') = <"+getHeaderValue(fileName,'TELDEC')+">")
        print("dmsToDeg(getHeaderValue(fileName,'TELDEC')) = ",dmsToDeg(getHeaderValue(fileName,'TELDEC')))
        print('"%.5f" % (dmsToDeg(getHeaderValue(fileName,\'TELDEC\'))) = ',"%.5f" % (dmsToDeg(getHeaderValue(fileName,'TELDEC'))))
        f.write("UPDATE `FitsFiles` SET `idPNMain` = '"+idPNMain+"', `observer` = 'Andreas Ritter', `instrument` = 'SpUpNIC', `telescope` = 'SAAO 1.9m', `RAJ2000` = '"+getHeaderValue(fileName,'TELRA')+"', `DECJ2000` = '"+getHeaderValue(fileName,'TELDEC')+"', `DRAJ2000` = "+"%.5f" % (dmsToDeg(getHeaderValue(fileName,'TELRA')))+", `DDECJ2000` = "+"%.5f" % (dmsToDeg(getHeaderValue(fileName,'TELDEC')))+", `convToText` = 'y' WHERE `fileName` = '"+fileName[fileName.rfind('/')+1:]+"';\n")

with open(sqlCatFile,'w') as f:
    f.write("CREATE TABLE IF NOT EXISTS MainPNData.SAAO_June2024 (\n")
    f.write("idSAAO_June2024 INT AUTO_INCREMENT PRIMARY KEY,\n")
    f.write("idPNMain INT NOT NULL,\n")
    f.write("mapFlag VARCHAR(1) NOT NULL);\n")
    #INSERT INTO `MainPNData`.`DataInfo` (`idDataInfo`, `Name`, `CatName`, `CatTitle`, `TabName`, `TabTitle`, `Mapped`, `MappedTo`, `MapKey`, `Date`, `Link`, `Comments`, `checked`, `Catalogue`, `checkNew`, `full`) VALUES ('287', 'SAAO_June2024', 'SAAO_June2024', 'SAAO_June2024', 'Table 1', 'SAAO_June2024', 'full', 'MainGPN.PNMain', 'ids', '2024-06-02', '-', '-', 'y', 'y', 'y', 'y');
    for i in np.arange(len(idPNMains)-1,-1,-1):
        if idPNMains[i] in idPNMains[0:i]:
            del idPNMains[i]
    for i in range(len(idPNMains)):
        f.write("INSERT INTO `MainPNData`.`SAAO_June2024`(`idSAAO_June2024`,`idPNMain`,`mapFlag`) VALUES (%d,%d,'y');\n" % (i+49,int(idPNMains[i])))
print('idPNMains = ',idPNMains)
