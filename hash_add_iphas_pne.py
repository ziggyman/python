import numpy as np

import csvFree
import csvData
from myUtils import hmsToDeg, dmsToDeg, angularDistance

inFileName = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort.csv'
hashInFileName = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_check_hash'
hashOutFileName = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_hash_output.csv'
hashAngDiamFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbAngDiam_240920.csv'
hashAngExtFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbAngExt_290920.csv'
hashPNMainFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_290920.csv'

minDistFile = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_hashDists.csv'
newPNeFile = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_new_grouped.csv'
sqlCommandsFile = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_new_grouped_add_to_HASH.sql'
iphasTable = '/Users/azuri/daten/uni/HKU/HASH/IPHAS_listALL2_MASPN_sort_new_grouped_table.sql'

minDist = 100.

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

def addNewPNeToHASH():
    inputData = csvFree.readCSVFile(newPNeFile)
    with open(sqlCommandsFile,'w') as f:
        with open(iphasTable,'w') as ft:
            f.write("USE `MainGPN`;\n")
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


if __name__ == '__main__':
#    createHashInFile()
#    checkHashOutFileRadii()
#    findInsidePN()
    groupNewPNeAndAverageCoordinates()
