import numpy as np

import csvData
import csvFree
from myUtils import hmsToDeg,dmsToDeg,raDecToLonLat

inList = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/Oct2021/NewCandidates_Tables_I+II_22102021.csv'

idPNMainStart = 33424
idtbCNamesStart = 122516
idtbAngDiamStart = 39803
idFRAStart = 11
#idtbUsrCommStart = 7839
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
    with open(inList[:inList.rfind('.')]+'.sql','w') as f:
        f.write('DROP TABLE IF EXISTS MainPNData.FrenchAmateurs_Oct2021;\n')
        f.write('CREATE TABLE IF NOT EXISTS MainPNData.FrenchAmateurs_Oct2021 (\n')
        f.write('idFrenchAmateurs_Oct2021 INT AUTO_INCREMENT PRIMARY KEY,\n')
        f.write('idPNMain INT,\n')
        f.write('mapFlag VARCHAR(1) NOT NULL\n')
        f.write(');\n')

        for i in np.arange(0,csv.size(),1):
            f.write("USE `MainGPN`;\n")
            hashIDs.append(idPNMainStart+i)
            raHMS = csv.getData("AD:(J2000)",i)
            print('raHMS = <'+raHMS+'>')
            ra = hmsToDeg(raHMS)
            decDMS = csv.getData("DEC (J2000)",i)
            dec = dmsToDeg(decDMS)
            print('raHMS = <'+raHMS+'> decDMS = <'+decDMS+'>')
            print('ra = ',ra,' dec = ',dec)
            lon, lat = raDecToLonLat(ra,dec)
            print('lon=',lon,', lat=',lat)
            png = ''
            if lon < 100:
                png += '0'
            if lon < 10:
                png += '0'
            png += str(lon)
            png = png[:png.rfind('.')+2]
#            print('png = <'+png+'>')
            if lat >= 0.:
                png += '+'
#                print('png = <'+png+'>')
            png += str(lat)
#            print('png = <'+png+'>')
            png = png[:png.rfind('.')+2]
#            print('png = <'+png+'>')

            if (lat < 10.) and (lat >= 0.):
                png = png[:png.rfind('+')+1]+'0'+png[png.rfind('+')+1:]
            if (lat  > -10.) and (lat < 0.):
                png = png[:png.rfind('-')+1]+'0'+png[png.rfind('-')+1:]
#                print('png = <'+png+'>')
            print('PNG '+png)
#            STOP
            f.write("INSERT INTO `PNMain`(`idPNMain`,`PNG`,`refPNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDecJ2000`,")
            f.write("`Glon`,`Glat`,`refCoord`,`Catalogue`,`refCatalogue`,`userRecord`,`domain`,`refDomain`,`PNstat`,`refPNstat`,`refSimbadID`,`show`) ")
            f.write("VALUES (%d,'%s','%s','%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n" % (idPNMainStart+i,
                                                                                                                                  png,
                                                                                                                                  'sys',
                                                                                                                                  raHMS,
                                                                                                                                  decDMS,
                                                                                                                                  ra,
                                                                                                                                  dec,
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

            majDiamStr = csv.getData("Dimension en minute d'arc (')",i).strip('~').rstrip('"').rstrip(':')
            if 'x' in majDiamStr:
                majDiamStr,minDiamStr = majDiamStr.split(' x ')
                print('majDiamStr = <',majDiamStr,'>, minDiamStr = <',minDiamStr,'>')
                if float(majDiamStr) < float(minDiamStr):
                    minDiamStr,majDiamStr = [majDiamStr,minDiamStr]
            else:
                minDiamStr = ''
            if majDiamStr in ['stellar?','stellar']:
                majDiamStr = '1'
                minDiamStr = '1'
#            if 'x' in diamStr:
            print('majDiamStr = <',majDiamStr,'>, minDiamStr = <',minDiamStr,'>')
            if minDiamStr != '':
                f.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`MinDiam`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`tempflag`) ")
                f.write("VALUES (%d,%.1f,%.1f,'%s',%d,'%s','%s',%d,'%s');\n" % (idtbAngDiamStart+i,
                                                                                float(majDiamStr)*60.,
                                                                                float(minDiamStr)*60.,
                                                                                'FrenchAmateurs',
                                                                                1,
                                                                                'sys',
                                                                                'ziggy',
                                                                                idPNMainStart+i,
                                                                                'n'))
            else:
                f.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`tempflag`) ")
                f.write("VALUES (%d,%.1f,'%s',%d,'%s','%s',%d,'%s');\n" % (idtbAngDiamStart+i,
                                                                           float(majDiamStr)*60.,
                                                                           'FrenchAmateurs',
                                                                           1,
                                                                           'sys',
                                                                           'ziggy',
                                                                           idPNMainStart+i,
                                                                           'n'))

            f.write("INSERT INTO `PNMain_tbAngDiam`(`idPNMain`,`idtbAngDiam`) VALUES (%d,%d);\n" % (idPNMainStart+i,idtbAngDiamStart+i))

#            name=csv.getData('Name',i)
#            if name[:4] == 'Pa J':
#                name = name[3:]
#            notes = csv.getData('comment1',i)
#            if notes == []:
#                notes = csv.getData('comment2',i)
#            else:
#                notes += ', '+csv.getData('comment2',i)
#            print('notes = <'+notes+'>')
#            f.write("INSERT INTO `tbUsrComm`(`idtbUsrComm`,`idPNMain`,`user`,`public`,`comment`,`date`) ")
#            f.write("VALUES (%d,%d,'%s','%s','%s','%s');\n" % (idtbUsrCommStart+i,
#                                                               idPNMainStart+i,
#                                                               'ziggy',
#                                                               'y',
#                                                               notes,
#                                                               '2020-03-31 19:30:00'))

            f.write("USE `MainPNData`;\n")
            f.write("INSERT INTO `FrenchAmateurs_Oct2021`(`idFrenchAmateurs_Oct2021`,`idPNMain`,`mapFlag`) ")
            f.write("VALUES (%d,%d,'%s');\n"
                                                                         % (idFRAStart+i,
                                                                            idPNMainStart+i,
                                                                            'y',
                                                                            ))

    with open(inList[:inList.rfind('.')]+'.fetch','w') as f:
        f.write('hashpn fetch all '+str(hashIDs[0]))
        for id in np.arange(1,len(hashIDs),1):
            f.write(','+str(hashIDs[id]))
        f.write(' -w force\n')
        f.write('hashpn brew all '+str(hashIDs[0]))
        for id in np.arange(1,len(hashIDs),1):
            f.write(','+str(hashIDs[id]))
        f.write(' -w\n')


if __name__ == "__main__":
    createMySQLCommands()

