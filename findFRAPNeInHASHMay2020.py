import numpy as np

import csvData
import csvFree
from myUtils import hmsToDeg,dmsToDeg,raDecToLonLat

longList = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/May2020/TableI_605PNG_21052020.csv'#'/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/TableI_569PNG_10042020_HASH.csv'
shortList = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/May2020/TableI I_ 124objets_21052020.csv'#'/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/TableI I_ 120objets_25032020.csv'

outList = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/May2020/FRA_all.csv'
hashOutList = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/May2020/FRAall_HASHoutput.csv'

idPNMainStart = 32567
idtbCNamesStart = 121506
idtbAngDiamStart = 39344
hashIDs = []

def readFRAFile(fname):
    csv = csvFree.readCSVFile(fname,',',False)
    print('csv.header[0] = <'+csv.header[0]+'>')
#    if csv.header[0] == '\ufeffNOM':
#        csv.renameColumn('\ufeffNOM','NOM')
#        print('renamed column "\ufeffNOM" to <'+csv.header[0]+'>')
#    if csv.header[len(csv.header)-1] == 'HASH ID\r':
#        csv.renameColumn('HASH ID\r','HASH ID')
#        print('renamed column "HASH ID\r" to <'+csv.header[len(csv.header)-1]+'>')
    print('csv = ',csv)
    print('csv.header = ',csv.header)
    return csv

def combineLists():
    keepColumns= ['NOM','AD:(J2000)','DEC (J2000)']
    print('keepColumns = ',keepColumns)
    csvAll = None
    for iList in [longList, shortList]:
        csv = readFRAFile(iList)
        print('csv.header = ',csv.header)
#        print('csv.data = ',csv.data)
#        STOP
        for keyword in csv.header:
            print('checking keyword <'+keyword+'>')
            if keyword not in keepColumns:
                print('keyword <'+keyword+'> not found in keepColumns => removing column')
                csv.removeColumn(keyword)
        print('csv.header = ',csv.header)
        print('csv.data = ',csv.data)

        for i in np.arange(0,csv.size(),1):
            print('changing ',csv.getData('NOM',i))
            csv.setData('NOM',i,csv.getData('NOM',i).replace(' ',''))
            print('... to ',csv.getData('NOM',i))

        if not csvAll:
            csvAll = csv
        else:
            csvAll.append(csv)
        print('csvAll.size() = ',csvAll.size())

    csvFree.writeCSVFile(csvAll, outList)

def findPNeNotInHASH():
    csvHASH = csvFree.readCSVFile(hashOutList)
    print('csvAll.size() = ',csvHASH.size())
    print('csvAll.header = ',csvHASH.header)

    found = []
    for i in np.arange(0,csvHASH.size(),1):
        if csvHASH.getData('pndb',i) != '':
            found.append(i)
    print('found = ',found)
    found.reverse()
    print('found = ',found)
    for i in found:
        csvHASH.removeRow(i)
    print('csvAll.size() = ',csvHASH.size())

    keepColumns= ["NOM","AD:(J2000)","DEC (J2000)","Dimension en minute d'arc (')","Coordonnées galactiques"]
    print('keepColumns = ',keepColumns)
    csvOut = None
    for iList in [longList, shortList]:
        csv = readFRAFile(iList)
        for keyword in csv.header:
            print('checking keyword <'+keyword+'>')
            if keyword not in keepColumns:
                print('keyword <'+keyword+'> not found in keepColumns => removing column')
                csv.removeColumn(keyword)
        print('csv.header = ',csv.header)
        print('csv.data = ',csv.data)

        hashNames = csvHASH.getData('id')
        remove = []
        for i in np.arange(0,csv.size(),1):
            if csv.getData('NOM',i).replace(' ','') not in hashNames:
                remove.append(i)
        remove.reverse()
        for i in remove:
            csv.removeRow(i)

        if "Coordonnées galactiques" not in csv.header:
            csv.addColumn("Coordonnées galactiques")
        for i in np.arange(0,csv.size(),1):
            lon, lat = raDecToLonLat(hmsToDeg(csv.getData("AD:(J2000)",i)),dmsToDeg(csv.getData("DEC (J2000)",i)))
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
            csv.setData("Coordonnées galactiques",i,png)
        if not csvOut:
            csvOut = csv
        else:
            csvOut.append(csv)

    # convert diameters from arcmin to arcsec
    for i in np.arange(0,csvOut.size(),1):
        diamStr = csvOut.getData("Dimension en minute d'arc (')",i).rstrip(' ')
        print('diamStr = <'+diamStr+'>')
        diamStrOut = ''
        if 'x' in diamStr:
            diamA = float(diamStr[:diamStr.find(' ')]) * 60.
            diamB = float(diamStr[diamStr.rfind(' ')+1:]) * 60.
            diamStrOut = '%.1f x %.1f' % (diamA, diamB)
        else:
            diamStrOut = '%.1f' % (float(diamStr) * 60.)
        csvOut.setData("Dimension en minute d'arc (')",i,diamStrOut)

    # write output
    csvFree.writeCSVFile(csvOut,hashOutList[:hashOutList.rfind('.')]+'_not_in_HASH.csv')

def createMySQLCommands():
    csv = csvFree.readCSVFile(hashOutList[:hashOutList.rfind('.')]+'_not_in_HASH.csv')
    print('csv.header = ',csv.header)
    with open(hashOutList[:hashOutList.rfind('.')]+'_not_in_HASH.sql','w') as f:
        f.write("USE `MainGPN`;\n")
        for i in np.arange(0,csv.size(),1):
            hashIDs.append(idPNMainStart+i)
            lon, lat = raDecToLonLat(hmsToDeg(csv.getData("AD:(J2000)",i)),dmsToDeg(csv.getData("DEC (J2000)",i)))
            f.write("INSERT INTO `PNMain`(`idPNMain`,`PNG`,`refPNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDecJ2000`,")
            f.write("`Glon`,`Glat`,`refCoord`,`Catalogue`,`refCatalogue`,`userRecord`,`domain`,`refDomain`,`PNstat`,`refPNstat`,`refSimbadID`,`show`) ")
            f.write("VALUES (%d,'%s','%s','%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n" % (idPNMainStart+i,
                                                                                                                                  csv.getData('Coordonnes galactiques',i),
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
#    combineLists()
    findPNeNotInHASH()
    createMySQLCommands()

if __name__ == "__main__":
    main()

