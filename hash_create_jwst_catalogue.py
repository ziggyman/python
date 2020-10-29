import csv

from myUtils import angularDistance,hmsToDeg,dmsToDeg

hash = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_221020.csv'
inputFileName = '/Users/azuri/daten/uni/HKU/HASH/JWST/list.txt'
csvFileName = inputFileName[:-3]+'csv'

foundObjectsFile = '/Users/azuri/daten/uni/HKU/HASH/JWST/listFound.csv'
sqlFile = '/Users/azuri/daten/uni/HKU/HASH/JWST/listFound.sql'
notInHashFile = '/Users/azuri/daten/uni/HKU/HASH/JWST/listNotFound.csv'

def makeCSVFile():
    with open(inputFileName,'r') as f:
        lines = f.readlines()
    #lines = [line.strip() for line in lines]
    lines = [line.replace('\t',',') for line in lines]
    with open(csvFileName,'w') as f:
        for line in lines:
            print(line)
            f.write(line)

def findInHASH():
    N_minDist_idPNMain = []
    with open(csvFileName,'r') as f:
        listReader = csv.DictReader(f)
        for row in listReader:
            N = row['N']
            ra1 = hmsToDeg(row['ICRS (J2000) RA'].replace(' ',':'))
            dec1 = dmsToDeg(row['ICRS (J2000) DEC'].replace(' ',':'))
            print("row['Otype'] = ",row['Otype'],': RA = ',ra1,', DEC = ',dec1)
            minDist = 100000000
            idPNMain = '-1'
            with open(hash,'r') as fHash:
                hashReader = csv.DictReader(fHash)
                for hashRow in hashReader:
                    ra2 = float(hashRow['DRAJ2000'])
                    dec2 = float(hashRow['DDECJ2000'])
#                    print('N = ',N,': ra2 = ',ra2,', dec2 = ',dec2)
                    dist = angularDistance(ra1,dec1,ra2,dec2)
                    if dist < minDist:
                        minDist = dist
                        idPNMain = hashRow['idPNMain']
                        print('N = ',N,': minDist = ',minDist)
                print('N = ',N,': minDist = ',minDist,', idPNMain = ',idPNMain)
                N_minDist_idPNMain.append([N,minDist,idPNMain])
    return N_minDist_idPNMain

def createCatalog():
    i=1
    with open(sqlFile,'w') as f:
        f.write('CREATE TABLE IF NOT EXISTS MainPNData.JWST_targets (\n')
        f.write('idJWST_targets INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
        f.write('N INT,\n')
        f.write('idPNMain INT,\n')
        f.write('RAJ2000 VARCHAR(25) NOT NULL,\n')
        f.write('DECJ2000 VARCHAR(25) NOT NULL,\n')
        f.write('DRAJ2000 FLOAT NOT NULL,\n')
        f.write('DDECJ2000 FLOAT NOT NULL,\n')
        f.write('Glon FLOAT NOT NULL,\n')
        f.write('Glat FLOAT NOT NULL,\n')
        f.write('PNMainDist FLOAT NOT NULL\n,')
        f.write('mapFlag VARCHAR(1) NOT NULL\n')
        f.write(');\n')
        f.write("USE `MainPNData`;\n")
        with open(foundObjectsFile,'r') as fIn:
            reader = csv.DictReader(fIn)
            for row in reader:
                with open(hash,'r') as fHash:
                    hashReader = csv.DictReader(fHash)
                    for hashRow in hashReader:
                        if hashRow['idPNMain'] == row['idPNMain']:
                            f.write("INSERT INTO `JWST_targets`(`idJWST_targets`,`N`,`idPNMain`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDECJ2000`,`Glon`,`Glat`,`PNMainDist`,`mapflag`) ")
                            f.write("VALUES (%d,%d,%d,'%s','%s',%.5f,%.5f,%.5f,%.5f,%.5f,'%s');\n" % (i,
                                                                                                      int(row['N']),
                                                                                                      int(row['idPNMain']),
                                                                                                      hashRow['RAJ2000'],
                                                                                                      hashRow['DECJ2000'],
                                                                                                      float(hashRow['DRAJ2000']),
                                                                                                      float(hashRow['DDECJ2000']),
                                                                                                      float(hashRow['Glon']),
                                                                                                      float(hashRow['Glat']),
                                                                                                      float(row['minDist']),
                                                                                                      'y'))
                            i += 1

def findObjectsNotInHash():
    with open(notInHashFile,'w') as fOut:
        fOut.write('N,type\n')
        with open(inputFileName[:-3]+'csv','r') as fIn:
            reader = csv.DictReader(fIn)
            for row in reader:
                found = False
                with open(foundObjectsFile,'r') as fFound:
                    foundReader = csv.DictReader(fFound)
                    for foundRow in foundReader:
                        if row['N'] == foundRow['N']:
                            found = True
                if not found:
                    fOut.write(row['N']+','+row['Otype']+'\n')

if __name__ == '__main__':
    if False:
        makeCSVFile()
        N_minDist_idPNMain = findInHASH()
        with open(foundObjectsFile,'w') as f:
            f.write('N,minDist,idPNMain\n')
            for row in N_minDist_idPNMain:
                f.write(row[0]+','+str(row[1])+','+row[2]+'\n')
    createCatalog()
    findObjectsNotInHash()
