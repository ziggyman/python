import os
import numpy as np
import csvFree,csvData
import matplotlib.pyplot as plt
from drUtils import getHeaderValue,plotSpec,getWavelengthArr,getImageData
from myUtils import hmsToDeg,dmsToDeg

fitsFiles = csvFree.readCSVFile('/Users/azuri/spectra/saao/saao_nov-2016/hash_FitsFiles_041024.csv')
fileNames = np.array(fitsFiles.getData('fileName'))
print('fileNames = ',fileNames)

CNames = csvFree.readCSVFile('/Users/azuri/spectra/saao/saao_nov-2016/hash_CNames_041024.csv')

fitsDir = '/Users/azuri/spectra/saao/saao_nov-2016/saao_nov2016/reduced2/hash/'
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
        wLen = getWavelengthArr(fileName,0)
        spec = getImageData(fileName,0)
        if spec.shape[0] == 2:
            spec = spec[0][0]
        plt.plot(wLen,spec)
        plt.title(fName)
        plt.show()
        objName = fName[:fName.find('_')]
        if objName == 'Mz3-CS':
            objName = 'Mz3'
        elif objName == 'DS1-CS':
            objName = 'DS1'
        elif objName == 'FP-J1912':
            objName = 'FPJ1912-0331'
        elif objName == 'HaWe13cspn':
            objName = 'HaWe13'
        elif objName == 'HaWe13pn':
            objName = 'HaWe13'
        elif objName == 'Kn-63':
            objName = 'Kn63'
        elif objName == 'M5':
            objName = '[M95]5'
        elif objName == 'NGC1360norm':
            objName = 'NGC1360'
        elif objName == 'NGC220cspn':
            objName = 'NGC2022'
        elif objName == 'NGC220pn':
            objName = 'NGC2022'
        elif objName == 'NGC246norm':
            objName = 'NGC246'
        elif objName == 'PG0109':
            objName = 'PG0109+111'
        elif objName == 'Pa-155':
            objName = 'Pa155'
        elif objName == 'PrTm-1-419':
            objName = 'PRTM1'
        elif objName == 'RE0503':
            objName = 'REJ0503-289'
        elif objName == 'RE J0503-289':
            objName = 'REJ0503-289'
        elif objName == 'RXJ0122':
            objName = 'RXJ0122.9-7521'
        elif objName == 'K2-7cs':
            objName = 'K2-7'
        elif objName == 'NGC1360cs':
            objName = 'NGC1360'
        elif objName == 'NGC7293cs':
            objName = 'NGC7293'
        elif objName == 'PFP1cs':
            objName = 'PFP1'
        elif objName == 'Pa155cs':
            objName = 'Pa155'

        print('fName = ',fName)
        print('objName = ',objName)
        ind = np.where(fileNames == fName)[0]
        print('ind = ',ind)
        if fitsFiles.getData('outName',ind) == 'NULL':
            idx = np.where(objNames == objName.lower())[0]
            print('objName = ',objName,': idx = ',idx)
            idPNMain = CNames.getData('idPNMain',idx)
            idPNMains.append(idPNMain)
            print("getHeaderValue(fileName,'TELRA') = <"+getHeaderValue(fileName,'TELRA')+'>')
            print("getHeaderValue(fileName,'TELDEC') = <"+getHeaderValue(fileName,'TELDEC')+">")
            print("dmsToDeg(getHeaderValue(fileName,'TELDEC')) = ",dmsToDeg(getHeaderValue(fileName,'TELDEC')))
            print('"%.5f" % (dmsToDeg(getHeaderValue(fileName,\'TELDEC\'))) = ',"%.5f" % (dmsToDeg(getHeaderValue(fileName,'TELDEC'))))
            f.write("UPDATE `FitsFiles` SET `idPNMain` = '"+idPNMain+"', `observer` = 'Quentin Parker + Travis Stenborg', `reference` = 'SAAO_Nov2016', `object` = '" + CNames.getData('Name',idx)+"', `instrument` = 'SpUpNIC', `telescope` = 'SAAO 1.9m', `RAJ2000` = '"+getHeaderValue(fileName,'TELRA')+"', `DECJ2000` = '"+getHeaderValue(fileName,'TELDEC')+"', `DRAJ2000` = "+"%.5f" % (dmsToDeg(getHeaderValue(fileName,'TELRA')))+", `DDECJ2000` = "+"%.5f" % (dmsToDeg(getHeaderValue(fileName,'TELDEC')))+", `convToText` = 'y' WHERE `fileName` = '"+fileName[fileName.rfind('/')+1:]+"';\n")

with open(sqlCatFile,'w') as f:
    f.write("CREATE TABLE IF NOT EXISTS MainPNData.SAAO_Nov2016 (\n")
    f.write("idSAAO_Nov2016 INT AUTO_INCREMENT PRIMARY KEY,\n")
    f.write("idPNMain INT NOT NULL,\n")
    f.write("mapFlag VARCHAR(1) NOT NULL);\n")
    #INSERT INTO `MainPNData`.`DataInfo` (`idDataInfo`, `Name`, `CatName`, `CatTitle`, `TabName`, `TabTitle`, `Mapped`, `MappedTo`, `MapKey`, `Date`, `Link`, `Comments`, `checked`, `Catalogue`, `checkNew`, `full`) VALUES ('287', 'SAAO_June2024', 'SAAO_June2024', 'SAAO_June2024', 'Table 1', 'SAAO_June2024', 'full', 'MainGPN.PNMain', 'ids', '2024-06-02', '-', '-', 'y', 'y', 'y', 'y');
    for i in np.arange(len(idPNMains)-1,-1,-1):
        if idPNMains[i] in idPNMains[0:i]:
            del idPNMains[i]
    for i in range(len(idPNMains)):
        f.write("INSERT INTO `MainPNData`.`SAAO_Nov2016`(`idSAAO_Nov2016`,`idPNMain`,`mapFlag`) VALUES (%d,%d,'y');\n" % (i+1,int(idPNMains[i])))
print('idPNMains = ',idPNMains)
