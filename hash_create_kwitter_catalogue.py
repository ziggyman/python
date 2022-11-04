import numpy as np

import csvData
import csvFree
from myUtils import hmsToDeg,dmsToDeg,raDecToLonLat

hashTable = '/Users/azuri/daten/uni/HKU/HASH/hash_fitsfiles_211022.csv'
sqlFileOut = '/Users/azuri/daten/uni/HKU/HASH/add_kwitter_catalogue.sql'

csvHash = csvFree.readCSVFile(hashTable)

i = 1
with open(sqlFileOut,'w') as f:
    f.write('CREATE TABLE IF NOT EXISTS MainPNData.Kwitter (\n')
    f.write('idKwitter INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
    f.write('idPNMain INT,\n')
    f.write('mapFlag VARCHAR(1) NOT NULL\n')
    f.write(');\n')
    f.write("USE `MainPNData`;\n")
    for iHash in range(csvHash.size()):
        if csvHash.getData('setname',iHash) == 'Kwitter':
            f.write("INSERT INTO `Kwitter`(`idKwitter`,`idPNMain`,`mapflag`) ")
            f.write("VALUES (%d,%d,'%s');\n" % (i,
                                                int(csvHash.getData('idPNMain',iHash)),
                                                'y'))
            i += 1

referencesFName = '/Users/azuri/daten/uni/HKU/HASH/KarenKwitter/publications.tsv'
refs = csvFree.readCSVFile(referencesFName,'\t',False)

refTableFName = '/Users/azuri/daten/uni/HKU/HASH/KarenKwitter/References_and_saturation.csv'
refTable = csvFree.readCSVFile(refTableFName)
spectraLinksFName = '/Users/azuri/daten/uni/HKU/HASH/KarenKwitter/hash_spectraLinks.csv'
spectraLinks = csvFree.readCSVFile(spectraLinksFName)

with open('/Users/azuri/daten/uni/HKU/HASH/KarenKwitter/addReferences.sql','w') as sql:
    sql.write("USE `PNSpectra_Sources`;\n")

    print('csvHash.header = ',csvHash.header)
    nNewLines = 1
    for i in range(refTable.size()):
        pnName = refTable.getData('PN NAME',i).replace(' ','').lower()
        print('pnName = <'+pnName+'>')
        for j in range(csvHash.size()):
            if csvHash.getData('setname',j) == 'Kwitter':
                name = csvHash.getData('fileName',j)
                name = name[:name.find('_')].lower()
                if name == pnName:
                    idPNMain = csvHash.getData('idPNMain',j)
                    references = refTable.getData('Published in',i).replace(' ','').replace(';',',').split(',')
                    print('pnName = <'+pnName+'> : idPNMain = ',idPNMain,': references = ',references)
                    idx = spectraLinks.find('idPNMain',idPNMain)
                    print('idx = ',idx)
                    publications = []
                    if idx[0] >= 0:
                        publications = spectraLinks.getData('reference',idx)
                        print('found publication for idPNMain=',idPNMain,': ',publications)
                    for ref in references:
                        ads = refs.getData('ADS Bibcode',refs.find('Reference',ref)[0])
                        print('ads = <'+ads+'>')
                        if ads not in publications:
#                            sql.write("INSERT INTO `spectraLinks`(`idspectraLinks`,`idPNMain`,`reference`,`user`,`date`) ")
#                            sql.write("VALUES (%d,%d,'%s','ziggy','2022-10-24 15:24:09');\n" % (int(spectraLinks.getData('idspectraLinks',spectraLinks.size()-1))+nNewLines,
#                                                                                                int(idPNMain),
#                                                                                                ads))
                            nNewLines += 1
                        else:
                            print('found publication ',ads,' already in spectraLinks')
                    sql.write("UPDATE `PNSpectra_Sources`.`FitsFiles` SET `reference` = '"+ads+"' WHERE (`idPNMain` = '"+idPNMain+"' AND `setname` = 'Kwitter');\n")
