import numpy as np
import astropy.units as u
import csvFree,csvData
from myUtils import angularDistance, hmsToDeg, dmsToDeg

f1name = '/Users/azuri/daten/uni/HKU/interns_projects/simba/Weidman-table1.csv'
f2name = '/Users/azuri/daten/uni/HKU/interns_projects/simba/Weidmann2020-cs-parameters.txt'
sqlFileOut = '/Users/azuri/daten/uni/HKU/interns_projects/simba/Weidman-table1.sql'
hashFile = '/Users/azuri/daten/uni/HKU/HASH/PNMain_full_June-30-2020.csv'
csvHash = csvFree.readCSVFile(hashFile)
simbaFile = '/Users/azuri/daten/uni/HKU/interns_projects/simba/All list v1.csv'
#simbaFile = '/Users/azuri/daten/uni/HKU/interns_projects/simba/simba_table.csv'
csvSimba = csvFree.readCSVFile(simbaFile)

with open(f2name, encoding="utf8", errors='ignore') as f:
    lines = f.readlines()

csv2 = csvData.CSVData()
csv2.header = ['PNG','log g','ref log g','met','log T','ref log T','log (L_star/L_sun)','ref log(L_star/L_sun','mag','ref mag']
for iLine in np.arange(1,len(lines),1):
    data = []
    lines[iLine] = lines[iLine].replace('\r','')
    lines[iLine] = lines[iLine].replace('\n','')
    lines[iLine] = lines[iLine].replace('}','+-')
    lines[iLine] = lines[iLine].replace(' v ',' v')
    lines[iLine] = lines[iLine].replace(' r ',' r')
    lines[iLine] = lines[iLine].replace(' b ',' b')
    lines[iLine] = lines[iLine].replace(' i ',' i')
    data = lines[iLine].split(' ')
    if len(data) < len(csv2.header):
        data.insert(3,' ')
    print('line[',iLine,'] = ',lines[iLine])
    print('data = ',data)
    csv2.append(data)
print('csv2.size() = ',csv2.size())
csvFree.writeCSVFile(csv2,f2name[:-3]+'csv')

csv1 = csvFree.readCSVFile(f1name)
print('csv1 = ',csv1)
print(csv1.header)
head = csv1.header
print('head[',len(head)-1,'] = ',head[len(head)-1])
head[len(head)-1] = 'Ref. spectra'#csv1.header[len(csv1.header)-1].strip('\\r')
csv1.header = head#csv1.header[len(csv1.header)-1].strip('\\r')
print('csv1.header[',len(csv1.header)-1,'] = ',csv1.header[len(csv1.header)-1])
print(csv1.header)
for i in range(csv1.size()):
    csv1.setData('Ref. spectra',i,csv1.getData('Ref. spectra',i).replace('\r',''))
    csv1.setData('Name',i,csv1.getData('Name',i).replace('_','-'))
    csv1.setData('DECJ2000',i,csv1.getData('DECJ2000',i).replace('_','-'))
    csv1.setData('DECJ2000',i,csv1.getData('DECJ2000',i).replace(' ',':'))
    csv1.setData('RAJ2000',i,csv1.getData('RAJ2000',i).replace(' ',':'))
    print(csv1.getData(i))

csvFree.writeCSVFile(csv1,f1name[:-4]+'_clean.csv')

nFound = 0
nNotFound = 0
foundPNe = []
notFoundPNe1 = []
notFoundPNe2 = []
for i1 in range(csv1.size()):
    found = False
    for i2 in range(csv2.size()):
        if csv1.getData('PNG',i1) == csv2.getData('PNG',i2):
            nFound += 1
            found = True
            foundPNe.append(csv1.getData('PNG',i1))
    if not found:
        nNotFound += 1
        notFoundPNe1.append(csv1.getData('PNG',i1))
print('found ',nFound,' CSPN from table 1 in table 2, did not find ',nNotFound)
print('did not find ',notFoundPNe1,' from table 1 in table 2')

nFound = 0
nNotFound = 0
for i2 in range(csv2.size()):
    found = False
    for i1 in range(csv1.size()):
        if csv1.getData('PNG',i1) == csv2.getData('PNG',i2):
            nFound += 1
            found = True
    if not found:
        nNotFound += 1
        notFoundPNe2.append(csv2.getData('PNG',i2))
print('found ',nFound,' CSPN from table 2 in table 1, did not find ',nNotFound)
print('did not find ',notFoundPNe2,' from table 2 in table 1')

with open(sqlFileOut,'w') as f:
    f.write('CREATE TABLE IF NOT EXISTS MainPNData.Weidmann_2020_not_in_HASH (\n')
    f.write('idWeidmann_2020_not_in_HASH INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
    f.write('idPNMain INT,\n')
    f.write('PNG VARCHAR(25),\n')
    f.write('f_PNG VARCHAR(25),\n')
    f.write('Name VARCHAR(25) NOT NULL,\n')
    f.write('RAJ2000 VARCHAR(25) NOT NULL,\n')
    f.write('DECJ2000 VARCHAR(25) NOT NULL,\n')
    f.write('SpT VARCHAR(25),\n')
    f.write('ref_SpT VARCHAR(25),\n')
    f.write('SpT2 VARCHAR(25),\n')
    f.write('ref_SpT2 VARCHAR(25),\n')
    f.write('RB VARCHAR(25),\n')
    f.write('ref_Spectra VARCHAR(25)\n,')
    f.write('PNMainDist INT NOT NULL\n,')
    f.write('mapFlag VARCHAR(1) NOT NULL\n')
    f.write(');\n')

    f.write('CREATE TABLE IF NOT EXISTS MainPNData.Weidmann_2020_catalog_table1 (\n')
    f.write('idWeidmann_2020_catalog_table1 INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
    f.write('idPNMain INT,\n')
    f.write('PNG VARCHAR(25),\n')
    f.write('f_PNG VARCHAR(25),\n')
    f.write('Name VARCHAR(25) NOT NULL,\n')
    f.write('RAJ2000 VARCHAR(25) NOT NULL,\n')
    f.write('DECJ2000 VARCHAR(25) NOT NULL,\n')
    f.write('SpT VARCHAR(25),\n')
    f.write('ref_SpT VARCHAR(25),\n')
    f.write('SpT2 VARCHAR(25),\n')
    f.write('ref_SpT2 VARCHAR(25),\n')
    f.write('RB VARCHAR(25),\n')
    f.write('ref_Spectra VARCHAR(25)\n,')
    f.write('PNMainDist INT NOT NULL\n,')
    f.write('mapFlag VARCHAR(1) NOT NULL\n')
    f.write(');\n')

    f.write('CREATE TABLE IF NOT EXISTS MainPNData.Weidmann_2020_catalog_table2 (\n')
    f.write('idWeidmann_2020_catalog_table2 INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
    f.write('idPNMain INT,\n')
    f.write('PNG VARCHAR(25),\n')
    f.write('f_PNG VARCHAR(25),\n')
    f.write('log_g FLOAT,\n')
    f.write('ref_log_g VARCHAR(25),\n')
    f.write('met VARCHAR(25),\n')
    f.write('log_T FLOAT,\n')
    f.write('ref_log_T VARCHAR(25),\n')
    f.write('log(L_star/L_sun) FLOAT,\n')
    f.write('ref_log(L_star/L_sun VARCHAR(25),\n')
    f.write('mag VARCHAR(25),\n')
    f.write('ref_mag VARCHAR(25),\n')
    f.write('mapFlag VARCHAR(1) NOT NULL\n')
    f.write(');\n')

    f.write("USE `MainPNData`;\n")

    print('csv1.header = ',csv1.header)
    j = 0
    for i in range(csv1.size()):
        found = False
        for iSimba in range(csvSimba.size()):
#            print("csvSimba.getData('Name in the paper',",iSimba,") = <"+csvSimba.getData('Name in the paper',iSimba)+'>, looking for <'+csv1.getData('Name',i)+'>')
            if csvSimba.getData('Name in the paper',iSimba) == csv1.getData('Name',i):
                print('found ',csv1.getData('Name',i))
                found = True
                idPNMain = csvSimba.getData('Hash id',iSimba)
                nfound = 0
                if csvSimba.getData('Does this star have central star coodinate in HASH?',iSimba) == '':
                    for iHash in range(csvHash.size()):
                        if csvHash.getData('idPNMain',iHash) == idPNMain:
                            j += 1
                            nfound += 1
                            print('idPNMain = ',idPNMain)
                            dist = angularDistance(hmsToDeg(csvHash.getData('RAJ2000',iHash))*u.deg,
                                                   dmsToDeg(csvHash.getData('DECJ2000',iHash))*u.deg,
                                                   hmsToDeg(csv1.getData('RAJ2000',i))*u.deg,
                                                   dmsToDeg(csv1.getData('DECJ2000',i))*u.deg)
                            f.write("INSERT INTO `Weidmann_2020_not_in_HASH`(`idWeidmann_2020_not_in_HASH`,`idPNMain`,`PNG`,`f_PNG`,`Name`,`RAJ2000`,`DECJ2000`,`SpT`,`ref_SpT`,`SpT2`,`ref_SpT2`,`RB`,`ref_Spectra`,`PNMainDist`,`mapflag`)")
    #                        f.write("INSERT INTO `Weidmann_2020_catalog_table1`(`idWeidmann_2020_catalog_table1`,`idPNMain`,`PNG`,`f_PNG`,`Name`,`RAJ2000`,`DECJ2000`,`SpT`,`ref_SpT`,`SpT2`,`ref_SpT2`,`RB`,`ref_Spectra`,`PNMainDist`,`mapflag`)")
            #['PNG', 'f_', 'Name', 'RAJ2000', 'DECJ2000', 'SpT', 'SpT_Ref', 'SpT2', 'SpT2_Ref', 'RB', 'Ref spectra\r']
                            f.write("VALUES (%d,%d,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s',%d,'%s');\n" % (j,
                                                                                                                          int(idPNMain),
                                                                                                                          csv1.getData('PNG',i),
                                                                                                                          csv1.getData('f_',i),
                                                                                                                          csv1.getData('Name',i).rstrip(),
                                                                                                                          csv1.getData('RAJ2000',i),
                                                                                                                          csv1.getData('DECJ2000',i),
                                                                                                                          csv1.getData('SpT',i),
                                                                                                                          csv1.getData('SpT_Ref',i),
                                                                                                                          csv1.getData('SpT2',i),
                                                                                                                          csv1.getData('SpT2_Ref',i),
                                                                                                                          csv1.getData('RB',i),
                                                                                                                          csv1.getData('Ref. spectra',i),
                                                                                                                          int(dist),
                                                                                                                          'y'))
                    print('nfound = ',nfound)
        if not found:
            print('ERROR: could not find <'+csv1.getData('Name',i)+'>')
            STOP

if False:
    print('csv2.header = ',csv2.header)
    for i in range(csv2.size()):
        found = False
        for iSimba in range(csvSimba.size()):
            if csvSimba.getData('Name in the paper',iSimba) == csv2.getData('Name',i):
                print('found ',csv1.getData('Name',i))
                found = True
                idPNMain = csvSimba.getData('Hash id',iSimba)
                for iHash in range(csvHash.size()):
                    if csvHash.getData('idPNMain',iHash) == idPNMain:
                        dist = angularDistance(hmsToDeg(csvHash.getData('RAJ2000',iHash))*u.deg,
                                               dmsToDeg(csvHash.getData('DECJ2000',iHash))*u.deg,
                                               hmsToDeg(csv1.getData('RAJ2000',i))*u.deg,
                                               dmsToDeg(csv1.getData('DECJ2000',i))*u.deg)

                        f.write("INSERT INTO `Weidmann_2020_catalog_table1`(`idWeidmann_2020_catalog_table1`,`idPNMain`,`PNG`,`f_PNG`,`Name`,`RAJ2000`,`DECJ2000`,`SpT`,`ref_SpT`,`SpT2`,`ref_SpT2`,`RB`,`ref_Spectra`,`PNMainDist`,`mapflag`)")
        #['PNG', 'f_', 'Name', 'RAJ2000', 'DECJ2000', 'SpT', 'SpT_Ref', 'SpT2', 'SpT2_Ref', 'RB', 'Ref spectra\r']
                        f.write("VALUES (%d,%d,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s',%d,'%s');\n" % (i+1,
                                                                                                                      int(idPNMain),
                                                                                                                      csv1.getData('PNG',i),
                                                                                                                      csv1.getData('f_',i),
                                                                                                                      csv1.getData('Name',i).rstrip(),
                                                                                                                      csv1.getData('RAJ2000',i),
                                                                                                                      csv1.getData('DECJ2000',i),
                                                                                                                      csv1.getData('SpT',i),
                                                                                                                      csv1.getData('SpT_Ref',i),
                                                                                                                      csv1.getData('SpT2',i),
                                                                                                                      csv1.getData('SpT2_Ref',i),
                                                                                                                      csv1.getData('RB',i),
                                                                                                                      csv1.getData('Ref. spectra',i),
                                                                                                                      int(dist),
                                                                                                                      'y'))
        if not found:
            print('ERROR: could not find ',csv1.getData('Name',i))
            STOP
