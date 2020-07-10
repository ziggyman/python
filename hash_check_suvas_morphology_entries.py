import csv

pnMorphFName = '/Users/azuri/daten/uni/HKU/interns_projects/suvas/hash-PNMorph_July-09-2020.csv'
sqlFile = '/Users/azuri/daten/uni/HKU/interns_projects/suvas/removeSuvasEntries.sql'
onlySuvasFile = '/Users/azuri/daten/uni/HKU/interns_projects/suvas/onlySuvas.sql'

with open(onlySuvasFile,'w') as f:
    f.write('CREATE TABLE IF NOT EXISTS MainPNData.PNMorphChandraOnly (\n')
    f.write('idPNMorphChandraOnly INT UNIQUE AUTO_INCREMENT PRIMARY KEY,\n')
    f.write('idPNMain INT UNIQUE,\n')
    f.write('mapFlag  VARCHAR(1) NOT NULL\n')
    f.write(');\n')
    f.write('USE MainPNData;\n')
with open(sqlFile,'w') as fs:
    fs.write('USE MainGPN;\n')

pnMorph = csv.DictReader(open(pnMorphFName))

iSuvas = 1

for row in pnMorph:
    if row['userRecord'] == 'chandra0101':
        temp = csv.DictReader(open(pnMorphFName))
        nFound = 0
        idsFound = []
        for tempRow in temp:
            if (tempRow['idPNMain'] == row['idPNMain']) and (tempRow['idtbPNMorph'] != row['idtbPNMorph']):
                #print('row = ',row)
                #print('tempRow = ',tempRow)
                nFound += 1
                idsFound.append(tempRow['idtbPNMorph'])

        if nFound > 0:
            with open(sqlFile,'a') as fs:
                fs.write('DELETE FROM `tbPNMorph` WHERE `idtbPNMorph` = %d;\n' % int(row['idtbPNMorph']))
        if nFound == 0:
            with open(onlySuvasFile,'a') as f:
                f.write("INSERT INTO `PNMorphChandraOnly`(`idPNMorphChandraOnly`,`idPNMain`,`mapflag`)\n")
                f.write("VALUES (%d,%d,'%s');\n" % (iSuvas,
                                                    int(row['idPNMain']),
                                                    'y'))
            iSuvas += 1
        elif nFound == 1:
            with open(sqlFile,'a') as f:
                f.write('UPDATE `tbPNMorph` SET `InUse` = 1 WHERE `idtbPNMorph` = %d;\n' % (int(idsFound[0])))
        else:
            print('found idPNMain=',row['idPNMain'],' ',nFound,' times:')
            isActive = 0
            for id in idsFound:
                temp = csv.DictReader(open(pnMorphFName))
                for tempRow in temp:
                    if (tempRow['idtbPNMorph'] == id) and (tempRow['InUse'] == '1') and (tempRow['userRecord'] != 'chandra0101'):
                        isActive += 1
                        print('found ',id,' active at ',tempRow['idtbPNMorph'])
            if isActive == 0:
                print('found idPNMain ',row['idPNMain'],' also at idtbPNMorph = ',id)
            elif isActive > 1:
                print('PROBLEM: More than 1 entry for idPNMain = ',row['idPNMain'],' is active')
