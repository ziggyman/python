import numpy as np

import csvData
import csvFree
from myUtils import hmsToDeg,dmsToDeg,raDecToLonLat

inList = '/Users/azuri/daten/uni/HKU/HASH/DSH/DSH_PNe_20160303.csv'

idPNMainStart = 32404
idtbCNamesStart = 121338
idtbAngDiamStart = 39150
idDSHPNeStart = 99
idtbUsrCommStart = 7839
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
        for i in np.arange(0,csv.size(),1):
            f.write("USE `MainGPN`;\n")
            hashIDs.append(idPNMainStart+i)
            raHMS = csv.getData("RAJ2000",i)
            ra = hmsToDeg(raHMS)
            decDMS = csv.getData("DECJ2000",i)
            dec = dmsToDeg(decDMS)
            refined = csv.getData('Refined coordinates',i)
            refined = refined.replace('  ',' ')
            if refined != '':
                print('refined = <'+refined+'>')
                raH,raM,raS, decD,decM,decS = refined.split(' ')
                raHMS = raH+':'+raM+':'+raS
                ra = hmsToDeg(raHMS)
                decDMS = decD+':'+decM+':'+decS
                dec = dmsToDeg(decDMS)
                print('refined coordinates: raHMS = <'+raHMS+'> decDMS = <'+decDMS+'>')
                print('refined coordinates: ra = ',ra,' dec = ',dec)
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
                                                                                                                                  'DSHPNe',
                                                                                                                                  'DSHPNe',
                                                                                                                                  'ziggy',
                                                                                                                                  'ziggy',
                                                                                                                                  'Galaxy',
                                                                                                                                  'ziggy',
                                                                                                                                  csv.getData('PNstat',i),
                                                                                                                                  'ziggy',
                                                                                                                                  'sys',
                                                                                                                                  'y'))

            f.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`simbadID`,`flag`)")
            f.write("VALUES (%d,'%s','%s',%d,'%s','%s',%d,'%s','%s');\n" % (idtbCNamesStart + i,
                                                                            csv.getData('Name',i),
                                                                            'DSHPNe',
                                                                            1,
                                                                            'sys',
                                                                            'ziggy',
                                                                            idPNMainStart+i,
                                                                            'n',
                                                                            'n'))

            f.write("INSERT INTO `PNMain_tbCNames`(`idPNMain`,`idtbCNames`) VALUES (%d,%d);\n" % (idPNMainStart+i,idtbCNamesStart+i))

            majDiamStr = csv.getData("MajDiam",i).strip('~').rstrip('"').rstrip(':')
            minDiamStr = csv.getData("MinDiam",i).strip('~').rstrip('"').rstrip(':')
            if majDiamStr in ['stellar?','stellar']:
                majDiamStr = '1'
                minDiamStr = '1'
#            if 'x' in diamStr:
            f.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`MinDiam`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`tempflag`) ")
            f.write("VALUES (%d,%.1f,%.1f,'%s',%d,'%s','%s',%d,'%s');\n" % (idtbAngDiamStart+i,
                                                                            float(majDiamStr),
                                                                            float(minDiamStr),
                                                                            'DSHPNe',
                                                                            1,
                                                                            'sys',
                                                                            'ziggy',
                                                                            idPNMainStart+i,
                                                                            'n'))
#            else:
#                f.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`tempflag`) ")
#                f.write("VALUES (%d,%.1f,'%s',%d,'%s','%s',%d,'%s');\n" % (idtbAngDiamStart+i,
#                                                                           float(diamStr),
#                                                                           'FrenchAmateurs',
#                                                                           1,
#                                                                           'sys',
#                                                                           'ziggy',
#                                                                           idPNMainStart+i,
#                                                                           'n'))

            f.write("INSERT INTO `PNMain_tbAngDiam`(`idPNMain`,`idtbAngDiam`) VALUES (%d,%d);\n" % (idPNMainStart+i,idtbAngDiamStart+i))

            name=csv.getData('Name',i)
            if name[:4] == 'Pa J':
                name = name[3:]
            notes = csv.getData('comment1',i)
            if notes == []:
                notes = csv.getData('comment2',i)
            else:
                notes += ', '+csv.getData('comment2',i)
            print('notes = <'+notes+'>')
            f.write("INSERT INTO `tbUsrComm`(`idtbUsrComm`,`idPNMain`,`user`,`public`,`comment`,`date`) ")
            f.write("VALUES (%d,%d,'%s','%s','%s','%s');\n" % (idtbUsrCommStart+i,
                                                               idPNMainStart+i,
                                                               'ziggy',
                                                               'y',
                                                               notes,
                                                               '2020-03-31 19:30:00'))

            f.write("USE `MainPNData`;\n")
            f.write("INSERT INTO `DSHPNe`(`idDSHPNe`,`Discoverer`,`ID`,`PNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDECJ2000`,`Glon`,`Glat`,`MajDiam`,`MinDiam`,`status`,`discovery`,`narrowImag`,`broadband`,`echelle`,`notes`,`PNMainDist`,`mapFlag`,`idPNMain`) ")
            f.write("VALUES (%d,'%s','%s','%s','%s','%s',%.4f,%.4f,%.2f,%.2f,%d,%d,'%s','%s','%s','%s','%s','%s',%d,'%s',%d);\n"
                                                                         % (idDSHPNeStart+i,
                                                                            'Pa',
                                                                            name,
                                                                            png,
                                                                            raHMS,
                                                                            decDMS,
                                                                            ra,
                                                                            dec,
                                                                            lon,
                                                                            lat,
                                                                            float(majDiamStr),
                                                                            float(minDiamStr),
                                                                            'New Candidates',
                                                                            '',
                                                                            '',
                                                                            '',
                                                                            '',
                                                                            notes,
                                                                            -1,
                                                                            'y',
                                                                            idPNMainStart+i))

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

