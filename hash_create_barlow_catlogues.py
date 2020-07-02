import numpy as np

import csvData
import csvFree
from myUtils import hmsToDeg,dmsToDeg,raDecToLonLat

barlowTables = ['/Users/azuri/daten/uni/HKU/HASH/Barlow1.csv','/Users/azuri/daten/uni/HKU/HASH/Barlow2.csv']
hashTable = '/Users/azuri/daten/uni/HKU/HASH/PNMain_full_June-25-2020.csv'#April-28-2020.csv'#
cNamesTable = '/Users/azuri/daten/uni/HKU/HASH/hash-commonNames.csv'

csvHash = csvFree.readCSVFile(hashTable)
print('csvHash = ',csvHash)

cNames = csvFree.readCSVFile(cNamesTable)

def getRAandDECFromIPHAS(iphas):
    iphas = iphas[iphas.find(' ')+1:]
    iphasTmp = iphas.replace('J','')
    iphasTmp = iphasTmp.replace('X','')
    iphasmp = iphasTmp.replace('-','+')
    raStr = iphasmp[:iphasmp.find('+')]
    decStr = iphasTmp[iphasmp.find('+'):]
    print('iphas = <'+iphas+': raStr = <'+raStr+'>, decStr = <'+decStr+'>')
    ras = raStr[raStr.find('.')-2:]
    raStr = raStr[:raStr.find('.')-2]
    ram = raStr[-2:]
    rah = raStr[:-2]
    print('h:m:s = '+rah+':'+ram+':'+ras)
    if '.' in decStr:
        decs = decStr[decStr.find('.')-2:]
        decStr = decStr[:decStr.find('.')-2]
        decm = decStr[-2:]
        decd = decStr[:-2]
    else:
        decs = decStr[-2:]
        decStr = decStr[:-2]
        decm = decStr[-2:]
        decd = decStr[:-2]
#        if '-' in iphasTmp:
#            decd = '-'+decd
    print('d:m:s = '+decd+':'+decm+':'+decs)
    return [rah+':'+ram+':'+ras,decd+':'+decm+':'+decs]

print(getRAandDECFromIPHAS('IPHAS J123456.7-123456'))
print(getRAandDECFromIPHAS('IPHAS J123456.7+123456'))
print(getRAandDECFromIPHAS('IPHAS J123456.7-123456.7'))
print(getRAandDECFromIPHAS('IPHAS J123456.7+123456.7'))
print(getRAandDECFromIPHAS('IPHASX J123456.7-123456'))
print(getRAandDECFromIPHAS('IPHASX J123456.7+123456'))
print(getRAandDECFromIPHAS('IPHASX J123456.7-123456.7'))
print(getRAandDECFromIPHAS('IPHASX J123456.7+123456.7'))
#STOP
i='a'
iBarlow = 0
hashIDs = []
hashIDsNew = []
idPNMainStart = 32611
idtbCNamesStart = 121547

for b in barlowTables:
    iBarlow += 1

    nFound = 0
    nNotFound = 0
    csvBarlow = csvFree.readCSVFile(b)
    print('iBarlow = ',iBarlow,': header = ',csvBarlow.header)
    csvBarlow.header = [h.strip() for h in csvBarlow.header]
    print('iBarlow = ',iBarlow,': header = ',csvBarlow.header)
    with open('/Users/azuri/daten/uni/HKU/HASH/mysql_commands_add_barlow'+str(iBarlow)+'.sql','w') as f:
        with open('/Users/azuri/daten/uni/HKU/HASH/mysql_commands_add_barlow'+str(iBarlow)+'_table.sql','w') as fb:
            fb.write("USE `MainPNData`;\n")
            for iRow in range(csvBarlow.size()):
                found = False
                foundIDx = 0
                if csvBarlow.getData('Name',iRow)[0] in ['J','X']:
                    csvBarlow.setData('Name',iRow,'IPHAS '+csvBarlow.getData('Name',iRow))
                tmpName = csvBarlow.getData('Name',iRow)
                if 'X' in tmpName:
                    tmpName = tmpName.replace(' ','')
                    csvBarlow.setData('Name',iRow,tmpName[:tmpName.find('X')+1]+' '+tmpName[tmpName.find('X')+1:])
#                    print('csvBarlow.getData("Name",iRow) = <'+csvBarlow.getData("Name",iRow))
    #            print('row ',i,': PNG '+csvBarlow.getData('PNG',iRow),', name = ',csvBarlow.getData('Name',iRow))
                for iHash in range(csvHash.size()):
                    if csvHash.getData('PNG',iHash) == csvBarlow.getData('PNG',iRow):
                        found = True
                        idPNMain = int(csvHash.getData('idPNMain',iHash))
                        foundIDx = iHash
                        print('found PNG '+csvHash.getData('PNG',iHash)+', name = '+csvBarlow.getData('Name',iRow)+' in HASH at ID ',idPNMain)
                if not found:
                    for iName in range(cNames.size()):
                        if csvBarlow.getData('Name',iRow).rstrip() == cNames.getData('Name',iName):
                            idPNMain = int(cNames.getData('idPNMain',iName))
                            print(csvBarlow.getData('Name',iRow).rstrip(),' ',cNames.getData('Name',iName),': idPNMain = ',idPNMain)
                            foundIDx = -1
                            for iHash in range(csvHash.size()):
                                if int(csvHash.getData('idPNMain',iHash)) == idPNMain:
                                    foundIDx = iHash
                                    found = True
                                    print('found Barlow PNG '+csvBarlow.getData('PNG',iRow)+', HASH PNG = '+csvHash.getData('PNG',foundIDx)+', name = <'+cNames.getData('Name',iName)+'> at idPNMain = ',idPNMain)
                if not found:

                    ra, dec = getRAandDECFromIPHAS(csvBarlow.getData('Name',iRow))
                    print('ra = <'+ra+'>, dec = <'+dec+'>')
                    lon, lat = raDecToLonLat(hmsToDeg(ra),dmsToDeg(dec))
                    png = '%010.6f' % lon
                    png = png[:png.find('.')+2]
                    #png = png.zfill(3)
                    if lat > 0:
                        png = png+'+'
                    png = png + '%08.6g' % lat
                    png = png[:png.rfind('.')+2]
                    print('lon = ',lon,', lat = ',lat,', png = <'+png+'>, PNG = <'+csvBarlow.getData('PNG',iRow)+'>')
                    for iHash in range(csvHash.size()):
                        if csvHash.getData('PNG',iHash) == png:
                            found = True
                            idPNMain = int(csvHash.getData('idPNMain',iHash))
                            foundIDx = iHash
                            print('found PNG '+csvHash.getData('PNG',iHash)+', name = '+csvBarlow.getData('Name',iRow)+' in HASH at ID ',idPNMain)
                            hashIDsNew.append(idPNMain)
                if found:
                    nFound += 1
                    fb.write("INSERT INTO `Barlow%d`(`idBarlow%d_2020`,`idPNMain`,`PNG`,`Name`,`mapflag`)" % (iBarlow,iBarlow))
                    fb.write("VALUES (%d,'%d','%s','%s','%s');\n" % (nFound,
                                                                     idPNMain,
                                                                     csvBarlow.getData('PNG',iRow),
                                                                     csvBarlow.getData('Name',iRow).rstrip(),
                                                                     'y'))
                else:
                    nNotFound += 1
                    print('COULD NOT FIND PNG '+csvBarlow.getData('PNG',iRow),', name = ',csvBarlow.getData('Name',iRow))

                    # add to HASH
                    hashIDs.append(idPNMainStart+nNotFound)
                    hashIDsNew.append(idPNMainStart+nNotFound)

                    f.write("INSERT INTO `PNMain`(`idPNMain`,`PNG`,`refPNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDecJ2000`,")
                    f.write("`Glon`,`Glat`,`refCoord`,`Catalogue`,`refCatalogue`,`userRecord`,`domain`,`refDomain`,`PNstat`,`refPNstat`,`refSimbadID`,`show`) ")
                    f.write("VALUES (%d,'%s','%s','%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n" % (idPNMainStart+nNotFound,
                                                                                                                                          png,#csvBarlow.getData('PNG',iRow),
                                                                                                                                          'sys',
                                                                                                                                          ra,
                                                                                                                                          dec,
                                                                                                                                          hmsToDeg(ra),
                                                                                                                                          dmsToDeg(dec),
                                                                                                                                          lon,
                                                                                                                                          lat,
                                                                                                                                          'ziggy',
                                                                                                                                          'ziggy',
                                                                                                                                          'ziggy',
                                                                                                                                          'ziggy',
                                                                                                                                          'Galaxy',
                                                                                                                                          'ziggy',
                                                                                                                                          'c',
                                                                                                                                          'ziggy',
                                                                                                                                          'sys',
                                                                                                                                          'y'))

                    f.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`simbadID`,`flag`)")
                    f.write("VALUES (%d,'%s','%s',%d,'%s','%s',%d,'%s','%s');\n" % (idtbCNamesStart + nNotFound,
                                                                                    csvBarlow.getData('Name',iRow),
                                                                                    'ziggy',
                                                                                    1,
                                                                                    'sys',
                                                                                    'ziggy',
                                                                                    idPNMainStart+nNotFound,
                                                                                    'n',
                                                                                    'n'))

                    f.write("INSERT INTO `PNMain_tbCNames`(`idPNMain`,`idtbCNames`) VALUES (%d,%d);\n" % (idPNMainStart+nNotFound,idtbCNamesStart+nNotFound))

    print('nFound = ',nFound,', nNotFound = ',nNotFound)
    if nNotFound > 0:
        with open('/Users/azuri/daten/uni/HKU/HASH/mysql_commands_add_barlow'+str(iBarlow)+'_not_in_HASH.fetch','w') as f:
            f.write('hashpn fetch all '+str(hashIDs[0]))
            for id in np.arange(1,len(hashIDs),1):
                f.write(','+str(hashIDs[id]))
            f.write(' -w force\n')
            f.write('hashpn brew all '+str(hashIDs[0]))
            for id in np.arange(1,len(hashIDs),1):
                f.write(','+str(hashIDs[id]))
            f.write(' -w\n')

print(hashIDsNew)
