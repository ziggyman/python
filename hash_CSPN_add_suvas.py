import csv
from myUtils import hmsToDeg, dmsToDeg, raDecToLonLat

#fNameSuvas = '/Users/azuri/daten/uni/HKU/interns_projects/simba/hash-CSCoords_chandra0101.csv'
fNameSimba = '/Users/azuri/daten/uni/HKU/interns_projects/simba/hashCSCoords_add_chandra0101.csv'
sqlOut = '/Users/azuri/daten/uni/HKU/interns_projects/simba/sql_add_chandra0101.sql'
sqlOutSuvas = '/Users/azuri/daten/uni/HKU/interns_projects/simba/sql_add_suvas.sql'
#suvas = csv.DictReader(open(fNameSuvas))
simba = csv.DictReader(open(fNameSimba))

iSuvas = 1
idtbCSCoords = 7828

with open(sqlOut,'w') as f:
    f.write('USE MainGPN;\n')

with open(sqlOutSuvas,'w') as fs:
    fs.write('CREATE TABLE IF NOT EXISTS MainPNData.CSPN_suvas (\n')
    fs.write('idCSPN_suvas INT UNIQUE AUTO_INCREMENT PRIMARY KEY,\n')
    fs.write('idPNMain INT UNIQUE,\n')
    fs.write('mapFlag  VARCHAR(1) NOT NULL\n')
    fs.write(');\n')
    fs.write('USE MainPNData;\n')
for row in simba:
    if (row['userRecord'] == 'chandra0101') and (':' in row['CS_RA_Simba']):
        draCalc = hmsToDeg(row['CS_RA_Simba'])
        ddecCalc = dmsToDeg(row['CS_Dec_Simba'])
        lCalc, bCalc = raDecToLonLat(draCalc,ddecCalc)
        with open(sqlOut,'a') as f:
            f.write("INSERT INTO `tbCSCoords`(`idtbCSCoords`,`CS_RAJ2000`,`CS_DECJ2000`,`CS_DRAJ2000`,`CS_DDECJ2000`,`CS_Glon`,`CS_Glat`,`InUse`,`userRecord`,`idPNMain`)\n")
            f.write("VALUES (%d,'%s','%s',%.5f,%.5f,%.5f,%.5f,%d,'%s',%d);\n" % (idtbCSCoords,
                                                                                 row['CS_RA_Simba'],
                                                                                 row['CS_Dec_Simba'],
                                                                                 draCalc,
                                                                                 ddecCalc,
                                                                                 lCalc,
                                                                                 bCalc,
                                                                                 0,
                                                                                 'Simba',
                                                                                 int(row['idPNMain'])))
        with open(sqlOutSuvas,'a') as fs:
            fs.write("INSERT INTO `CSPN_suvas`(`idCSPN_suvas`,`idPNMain`,`mapflag`)\n")
            fs.write("VALUES (%d,%d,'%s');\n" % (iSuvas,
                                                int(row['idPNMain']),
                                                'y'))
        idtbCSCoords += 1
        iSuvas += 1


