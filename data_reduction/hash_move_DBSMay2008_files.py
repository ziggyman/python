import csvFree,csvData
hashFile = '/Users/azuri/daten/uni/HKU/HASH/hash_FitsFiles_010324.csv'
oldNamesFile = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/REDUCED/science.list'
newNamesFile = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/REDUCED/new.list'

outFile = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/REDUCED/renameFiles.sql'
idFile = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/REDUCED/ids.dat'

hashTable = csvFree.readCSVFile(hashFile)

with open(oldNamesFile,'r') as f:
    oldNames = f.readlines()
with open(newNamesFile,'r') as f:
    newNames = f.readlines()

oldNames = [name.strip() for name in oldNames]
newNames = [name.strip() for name in newNames]

with open(outFile,'w') as f:
    with open(idFile,'w') as id:
        for i in range(hashTable.size()):
            if hashTable.getData('setname',i) == 'DBS_May2008':
                oldName = hashTable.getData('fileName',i)
                oldNameFound = False
                idPNMain = hashTable.getData('idPNMain',i)
                id.write(idPNMain+',')
                for j in range(len(oldNames)):
                    if oldNames[j] == oldName:
                        oldNameFound = True
                        f.write("UPDATE `PNSpectra_Sources`.`FitsFiles` SET `fileName` = '%s' WHERE (`idFitsFiles` = '%s');\n" % (newNames[j],hashTable.getData('idFitsFiles',i)))
                        f.write("UPDATE `PNSpectra_Sources`.`FitsFiles` SET `outName` = NULL WHERE (`idFitsFiles` = '%s');\n" % (hashTable.getData('idFitsFiles',i)))
                        f.write("UPDATE `PNSpectra_Sources`.`FitsFiles` SET `parsedIn` = 'n' WHERE (`idFitsFiles` = '%s');\n" % (hashTable.getData('idFitsFiles',i)))
                        f.write("UPDATE `PNSpectra_Sources`.`FitsFiles` SET `telescope` = 'SSO' WHERE (`idFitsFiles` = '%s');\n" % (hashTable.getData('idFitsFiles',i)))
                if not oldNameFound:
                    print('ERROR: could not find oldName <'+oldName+'>')
                    STOP
