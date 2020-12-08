import numpy as np

import csvFree
import csvData
from myUtils import hmsToDeg, dmsToDeg, angularDistance, degToHMS, degToDMS, raDecToLonLat, getPNG

inFileName = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort.csv'
hashInFileName = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_check_hash'
hashOutFileName = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_hash_output.csv'
hashPNMain_tbAngDiamFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_tbAngDiam_111120.csv'
hashAngDiamFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbAngDiam_111120.csv'
hashPNMain_CNamesFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_tbCNames_111120.csv'
hashCNamesFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbCNames_111120.csv'
hashAngExtFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbAngExt_290920.csv'
hashPNMainFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_111120.csv'

minDistFile = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_hashDists.csv'
newPNeFile = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_new_grouped.csv'
newHashOutFile = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_new_grouped_hashout.csv'
sqlCommandsFile = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_new_grouped_add_to_HASH.sql'
iphasTable = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_new_grouped_table.sql'
catalogName = 'Sabin_IPHASPNe_Nov2020'
hashpnFile = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_new_grouped_table.hash'

minDist = 100.
csv = csvFree.readCSVFile(hashPNMainFileName)
idPNMainStart = int(csv.getData('idPNMain',csv.size()-1))+1
csv = csvFree.readCSVFile(hashPNMain_CNamesFileName)
idtbCNamesStart = int(csv.getData('idtbCNames',csv.size()-1)) + 1
csv = csvFree.readCSVFile(hashPNMain_tbAngDiamFileName)
idtbAngDiamStart = int(csv.getData('idtbAngDiam',csv.size()-1)) + 1

def createHashInFile():
    csvIphas = csvFree.readCSVFile(inFileName)
    prevRA = '10:10:10.0'
    prevDEC = '10:10:10.0'
    with open(hashInFileName,'w') as f:
        for i in range(csvIphas.size()):
            ra = csvIphas.getData('RA_H',i)+':'+csvIphas.getData('RA_M',i)+':'+csvIphas.getData('RA_S',i)
            dec = csvIphas.getData('DEC_D',i)+':'+csvIphas.getData('DEC_M',i)+':'+csvIphas.getData('DEC_S',i)
            print(i,': ra = ',ra,', dec = ',dec)
#            dist = angularDistance(hmsToDeg(ra), dmsToDeg(dec), hmsToDeg(prevRA), dmsToDeg(prevDEC))
#            distA = pyasl.getAngDist(hmsToDeg(ra), dmsToDeg(dec), hmsToDeg(prevRA), dmsToDeg(prevDEC))
#            print(i,': dist = ',dist,', distA = ',distA)
#            if dist > minDist:
            f.write(str(i)+','+ra+','+dec+'\n')#+','+str(dist)
#            prevRA = ra
#            prevDEC = dec

def getMinDistPN(ra,dec,PNMain):
    print('getMinDistPN: checking for ra = <'+ra+'>, dec = <'+dec+'>')
    print('type(ra) = ',type(ra),', type(dec) = ',type(dec))
    dra = hmsToDeg(ra)
    ddec = dmsToDeg(dec)
    print('dra = ',dra,', ddec = ',ddec)
    minDist = 10000.
    id = -1
    for i in range(PNMain.size()):
        dist = angularDistance(dra, ddec, hmsToDeg(PNMain.getData('RAJ2000',i)), dmsToDeg(PNMain.getData('DECJ2000',i)))
        if dist < minDist:
            minDist = dist
            id = PNMain.getData('idPNMain',i)
            print('closest PN so far: dist = ',minDist,', id = ',id)
    return [id, minDist]

def findInsidePN():
    angDiams = csvFree.readCSVFile(hashAngDiamFileName)
    angExt = csvFree.readCSVFile(hashAngExtFileName)
    PNMain = csvFree.readCSVFile(hashPNMainFileName)
    csvIphas = csvFree.readCSVFile(inFileName)

    with open(minDistFile,'w') as f:
        f.write('iphasRow,RAJ2000,DECJ2000,closestHASHPNid,distance,inside\n')
        for i in range(csvIphas.size()):
            ra = csvIphas.getData('RA_H',i)+':'+csvIphas.getData('RA_M',i)+':'+csvIphas.getData('RA_S',i)
            dec = csvIphas.getData('DEC_D',i)+':'+csvIphas.getData('DEC_M',i)+':'+csvIphas.getData('DEC_S',i)
            id, dist = getMinDistPN(ra, dec, PNMain)
            print('i = ',i,': closest PN id = ',id,': dist = ',dist)
            majDiamFound = False
            inside = False
            for j in range(angDiams.size()):
                if (angDiams.getData('idPNMain',j) == id) and (angDiams.getData('InUse',j) == '1'):
                    majDiamFound = True
                    if dist < (float(angDiams.getData('MajDiam',j)) / 2.):
                        inside = True
                        print('i = ',i,' is inside idPNMain ',id)
                    else:
                        print('i = ',i,' is NOT inside idPNMain ',id)
            f.write('%d,%s,%s,%s,%d,' % (i,ra,dec,id,int(dist)))
            if majDiamFound and inside:
                f.write('isInside')
            elif majDiamFound and not inside:
                f.write('notInside')
            else:
                f.write('unsure')
            f.write('\n')

def checkHashOutFileRadii():
    if False:
        csvHashList = csvFree.readCSVFile(hashOutFileName)
        for i in range(csvHashList.size()):
            idPNMain = csvHashList.getData('pndb',i)
            diam = 0.
            ext = 0.
            for j in range(angDiams.size()):
                if angDiams.getData('idPNMain',j) == idPNMain:
                    if angDiams.getData('InUse',j) == '1':
                        diam = float(angDiams.getData('MajDiam',j))
            for j in range(angExt.size()):
                if angExt.getData('idPNMain',j) == idPNMain:
                    if angExt.getData('InUse',j) == '1':
                        ext = float(angExt.getData('MajExt',j))
            print('idPNMain = ',idPNMain,': diam = ',diam,', ext = ',ext)

def groupNewPNeAndAverageCoordinates():
    minDists = csvFree.readCSVFile(minDistFile)

    idsAlreadyGrouped = []
    with open(newPNeFile,'w') as f:
        f.write('DRAJ2000,DDECJ2000\n')
        for i in range(minDists.size()):
            group = []
            dras = []
            ddecs = []
            if (minDists.getData('inside',i) != 'isInside') and (float(minDists.getData('distance',i)) > minDist):
                if minDists.getData('iphasRow',i) not in idsAlreadyGrouped:
                    idsAlreadyGrouped.append(minDists.getData('iphasRow',i))
                    ra = minDists.getData('RAJ2000',i)
                    dec = minDists.getData('DECJ2000',i)
                    dra = hmsToDeg(ra)
                    ddec = dmsToDeg(dec)
                    group.append(i)
                    dras.append(dra)
                    ddecs.append(ddec)
                    for j in np.arange(i+1,minDists.size(),1):
                        if minDists.getData('iphasRow',j) not in idsAlreadyGrouped:
                            ra = minDists.getData('RAJ2000',j)
                            dec = minDists.getData('DECJ2000',j)
                            thisdra = hmsToDeg(ra)
                            thisddec = dmsToDeg(dec)
                            if angularDistance(dra, ddec, thisdra, thisddec) < minDist:
                                group.append(minDists.getData('iphasRow',j))
                                idsAlreadyGrouped.append(minDists.getData('iphasRow',j))
                                dras.append(thisdra)
                                ddecs.append(thisddec)
                    print('found ',len(group),' different coordinates for apparently the same PN')
                    f.write('%.5f,%.5f\n' % (np.mean(dras),np.mean(ddecs)))

def makeHashFile():
    inputData = csvFree.readCSVFile(newPNeFile)
    with open(newPNeFile[:newPNeFile.rfind('.')]+'.hash','w') as fh:
        for i in range(inputData.size()):
            dra = float(inputData.getData('DRAJ2000',i))
            ddec = float(inputData.getData('DDECJ2000',i))
            lon, lat = raDecToLonLat(dra, ddec)
            ra = degToHMS(dra)
            dec = degToDMS(ddec)
            print('lon = ',lon,', lat = ',lat,', dra = ',dra,', ddec = ',ddec,', ra = ',ra,', dec = ',dec)
            png = getPNG(lon,lat)
            fh.write(png+','+ra+','+dec+'\n')

def getIPHASName(ra,dec):
    hRa,mRa,sRa = ra.split(':')
    dDec,mDec,sDec = dec.split(':')
    sign = '+' if int(dDec) >= 0 else ''
    return 'J'+hRa+mRa+sRa[:sRa.find('.')+3]+sign+dDec+mDec+sDec[:sDec.find('.')+2]

def addNewPNeToHASH():
    inputData = csvFree.readCSVFile(newPNeFile)
    PNMain = csvFree.readCSVFile(hashPNMainFileName)
    hashOut = csvFree.readCSVFile(newHashOutFile)
    pneInHash = hashOut.getData('id')
    pngs = PNMain.getData('PNG')
    pngsInHash = []
    with open(sqlCommandsFile,'w') as f:
        with open(iphasTable,'w') as ft:
            ft.write('CREATE TABLE IF NOT EXISTS MainPNData.'+catalogName+' (\n')
            ft.write('id'+catalogName+' INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
            ft.write('idPNMain INT,\n')
            ft.write('mapFlag VARCHAR(1) NOT NULL\n')
            ft.write(');\n')

            ft.write("USE `MainPNData`;\n")

            ids = []
            f.write("USE `MainGPN`;\n")
            for i in range(inputData.size()):
                dra = float(inputData.getData('DRAJ2000',i))
                ddec = float(inputData.getData('DDECJ2000',i))
                lon, lat = raDecToLonLat(dra, ddec)
                ra = degToHMS(dra)
                dec = degToDMS(ddec)
                print('lon = ',lon,', lat = ',lat)
                png = getPNG(lon,lat)
                print('png = <'+png+'>')
                if png in pngs:
                    print('PNG '+png+' already in HASH')
                    pngsInHash.append(png)
                    png = png+'a'
                if png in pngs:
                    pngsInHash.append(png)
                    print('PNG '+png+' already in HASH')
                    png = png[:-1]+'b'
                if png in pngs:
                    pngsInHash.append(png)
                    print('PNG '+png+' already in HASH')
                    png = png[:-1]+'c'
                if png in pngs:
                    pngsInHash.append(png)
                    print('PNG '+png+' already in HASH')
                    png = png[:-1]+'d'
                if png in pngs:
                    pngsInHash.append(png)
                    print('PNG '+png+' already in HASH')
                    png = png[:-1]+'e'

                if (png in pneInHash) and (hashOut.getData('pndb',pneInHash.index(png)) != ''):
                    print('PNG '+png,' found in pneInHash: pneInHash.index(',png,') = ',pneInHash.index(png))
                    idPNMain = int(hashOut.getData('pndb',pneInHash.index(png)))
                    # add IPHAS name to common names
                else:
                    idPNMain = idPNMainStart+i
                    ids.append(idPNMain)
                    f.write("INSERT INTO `PNMain`(`idPNMain`,`PNG`,`refPNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDecJ2000`,")
                    f.write("`Glon`,`Glat`,`refCoord`,`Catalogue`,`refCatalogue`,`userRecord`,`domain`,`refDomain`,`PNstat`,`refPNstat`,`refSimbadID`,`show`) ")
                    f.write("VALUES (%d,'%s','%s','%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n" % (idPNMain,
                                                                                                                                          png,#csvBarlow.getData('PNG',iRow),
                                                                                                                                          'sys',
                                                                                                                                          ra,
                                                                                                                                          dec,
                                                                                                                                          dra,
                                                                                                                                          ddec,
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

                ft.write("INSERT INTO `"+catalogName+"`(`idPNMain`,`mapflag`) ")
                ft.write("VALUES (%d,'%s');\n" % (idPNMain,
                                                  'y'))
                iphasName = getIPHASName(ra,dec)
                f.write("INSERT INTO `PNMain_tbCNames`(`idPNMain`,`idtbCnames`) VALUES (%d,%d);\n" % (idPNMain,idtbCNamesStart + i))
                f.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`simbadID`,`flag`) ")
                f.write("VALUES (%d,'%s',%d,'%s','%s',%d,'%s','%s');\n" % (idtbCNamesStart+1,
                                                                           iphasName,
                                                                           1,
                                                                           'sys',
                                                                           'sys',
                                                                           idPNMain,
                                                                           'n',
                                                                           'n'))
                f.write("INSERT INTO `PNMain_tbAngDiam`(`idPNMain`,`idtbAngDiam`) VALUES (%d,%d);\n" % (idPNMain,idtbAngDiamStart + i))
                f.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`InUse`,`userRecord`,`idPNMain`,`tempflag`) ")
                f.write("VALUES (%d,%.0f,%d,'%s',%d,'%s');\n" % (idtbAngDiamStart + i,
                                                                 300.,
                                                                 1,
                                                                 'sys',
                                                                 idPNMain,
                                                                 'n'))

    with open(hashpnFile,'w') as hf:
        for id in ids:
            hf.write('hashpn fetch all '+str(id)+' -w force\n')
            hf.write('hashpn brew all '+str(id)+' -w\n')
            hf.write('echo "finished HASH ID %d" >> logfile_IPHAS.log\n' % id)

    print('pngsInHash = ',pngsInHash)


if __name__ == '__main__':
#    createHashInFile()
#    checkHashOutFileRadii()
#    findInsidePN()
#    groupNewPNeAndAverageCoordinates()
#    makeHashFile()
    addNewPNeToHASH()
