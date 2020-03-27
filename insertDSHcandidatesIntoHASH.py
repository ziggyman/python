import numpy as np

import csvData
import csvFree
from myUtils import hmsToDeg,dmsToDeg,raDecToLonLat

inList = '/Users/azuri/daten/uni/HKU/HASH/DSH/DSH_PNe_20160303.csv'

idPNMainStart = 32404
idtbCNamesStart = 121338
idtbAngDiamStart = 39150
hashIDs = []

def readFRAFile(fname):
    csv = csvFree.readCSVFile(fname,',',False)
    if csv.header[0] == '\ufeffNOM':
        csv.renameColumn('\ufeffNOM','NOM')
        print('renamed column "\ufeffNOM" to <'+csv.header[0]+'>')
    if csv.header[len(csv.header)-1] == 'HASH ID\r':
        csv.renameColumn('HASH ID\r','HASH ID')
        print('renamed column "HASH ID\r" to <'+csv.header[len(csv.header)-1]+'>')
    print('csv = ',csv)
    print('csv.header = ',csv.header)
    return csv

def createMySQLCommands():
    csv = csvFree.readCSVFile(inList)
    with open(inList[:inList.rfind('.')]+'.mysql','w') as f:
        f.write("USE `MainGPN`;\n")
        for i in np.arange(0,csv.size(),1):
            hashIDs.append(idPNMainStart+i)
            ra = hmsToDeg(csv.getData("RAJ2000",i))
            dec = dmsToDeg(csv.getData("DECJ2000",i))
            refined = csv.getData('Refined coordinates',i)
            if refined != '':
                ra, dec = refined.split(' ')
            lon, lat = raDecToLonLat(ra,dec)
            f.write("INSERT INTO `PNMain`(`idPNMain`,`PNG`,`refPNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDecJ2000`,")
            f.write("`Glon`,`Glat`,`refCoord`,`Catalogue`,`refCatalogue`,`userRecord`,`domain`,`refDomain`,`PNstat`,`refPNstat`,`refSimbadID`,`show`) ")
            f.write("VALUES (%d,'%s','%s','%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n" % (idPNMainStart+i,
                                                                                                                                  csv.getData('Coordonnées galactiques',i),
                                                                                                                                  'sys',
                                                                                                                                  csv.getData("AD:(J2000)",i),
                                                                                                                                  csv.getData("DEC (J2000)",i),
                                                                                                                                  hmsToDeg(csv.getData("AD:(J2000)",i)),
                                                                                                                                  dmsToDeg(csv.getData("DEC (J2000)",i)),
                                                                                                                                  lon,
                                                                                                                                  lat,
                                                                                                                                  'FrenchAmateurs',
                                                                                                                                  'FrenchAmateurs',
                                                                                                                                  'ziggy',
                                                                                                                                  'ziggy',
                                                                                                                                  'Galaxy',
                                                                                                                                  'ziggy',
                                                                                                                                  'c',
                                                                                                                                  'ziggy',
                                                                                                                                  'sys',
                                                                                                                                  'y'))

            f.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`simbadID`,`flag`)")
            f.write("VALUES (%d,'%s','%s',%d,'%s','%s',%d,'%s','%s');\n" % (idtbCNamesStart + i,
                                                                            csv.getData('NOM',i),
                                                                            'FrenchAmateurs',
                                                                            1,
                                                                            'sys',
                                                                            'ziggy',
                                                                            idPNMainStart+i,
                                                                            'n',
                                                                            'n'))

            f.write("INSERT INTO `PNMain_tbCNames`(`idPNMain`,`idtbCNames`) VALUES (%d,%d);\n" % (idPNMainStart+i,idtbCNamesStart+i))

            diamStr = csv.getData("Dimension en minute d'arc (')",i).strip(' ').rstrip(' ')
            if 'x' in diamStr:
                f.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`MinDiam`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`tempflag`) ")
                f.write("VALUES (%d,%.1f,%.1f,'%s',%d,'%s','%s',%d,'%s');\n" % (idtbAngDiamStart+i,
                                                                                float(diamStr[:diamStr.find(' ')]),
                                                                                float(diamStr[diamStr.rfind(' ')+1:]),
                                                                                'FrenchAmateurs',
                                                                                1,
                                                                                'sys',
                                                                                'ziggy',
                                                                                idPNMainStart+i,
                                                                                'n'))
            else:
                f.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`tempflag`) ")
                f.write("VALUES (%d,%.1f,'%s',%d,'%s','%s',%d,'%s');\n" % (idtbAngDiamStart+i,
                                                                           float(diamStr),
                                                                           'FrenchAmateurs',
                                                                           1,
                                                                           'sys',
                                                                           'ziggy',
                                                                           idPNMainStart+i,
                                                                           'n'))

            f.write("INSERT INTO `PNMain_tbAngDiam`(`idPNMain`,`idtbAngDiam`) VALUES (%d,%d);\n" % (idPNMainStart+i,idtbAngDiamStart+i))
    with open(hashOutList[:hashOutList.rfind('.')]+'_not_in_HASH.fetch','w') as f:
        f.write('hashpn fetch all '+str(hashIDs[0]))
        for id in np.arange(1,len(hashIDs),1):
            f.write(','+str(hashIDs[id]))
        f.write(' -w force\n')
        f.write('hashpn brew all '+str(hashIDs[0]))
        for id in np.arange(1,len(hashIDs),1):
            f.write(','+str(hashIDs[id]))
        f.write(' -w\n')


def main():
    findPNeNotInHASH()
    createMySQLCommands()

if __name__ == "__main__":
    main()

