import csv
import numpy as np
import os
import pathlib
from PyAstronomy import pyasl
import unicodedata as ud

import csvFree,csvData
from myUtils import degToDMS,degToHMS,getPNGName,raDecToLonLat,dmsToDeg,hmsToDeg,angularDistance

spectraPath = '/Users/azuri/spectra/AAO_bulge'
subPaths = ['N2',
            'N3',
            'N4',
            'N5',
            'N7',
            'S18',
            'S19',
            'S20',
            'S21',
            'S22',
            'S26',
            'S3']

notFoundFileName = os.path.join(spectraPath,'notFound.list')
hashSearchFileName = os.path.join(spectraPath,'hash_search.csv')
hashAngDiamFileName = os.path.join(spectraPath,'hash_tbAngDiam.csv')
sqlFileOut = os.path.join(spectraPath,'hash_add_new.sql')
pnMainFile = os.path.join(spectraPath,'hash_PNMain.csv')
tbCNamesFile = os.path.join(spectraPath,'hash_tbCNames.csv')
raDecFileName = os.path.join(spectraPath,'raDec.csv')

def getNewPNGName(pngNames,lon,lat):
    pngNameTemp = getPNGName(lon,lat)
    pngName = pngNameTemp
    if pngName in pngNames:
        pngName = pngNameTemp+'a'
        if pngName in pngNames:
            pngName = pngNameTemp+'b'
            if pngName in pngNames:
                pngName = pngNameTemp+'c'
                if pngName in pngNames:
                    pngName = pngNameTemp+'d'
                    if pngName in pngNames:
                        pngName = pngNameTemp+'e'
                        if pngName in pngNames:
                            pngName = pngNameTemp+'f'
                            if pngName in pngNames:
                                pngName = pngNameTemp+'g'
                                if pngName in pngNames:
                                    pngName = pngNameTemp+'h'
                                    if pngName in pngNames:
                                        pngName = pngNameTemp+'i'
                                        if pngName in pngNames:
                                            pngName = pngNameTemp+'j'
    return pngName

def list2csv(path):
    listFile = list(pathlib.Path(path).glob('*.lis*'))[0].name
    print('listFile = ',listFile)
    with open(os.path.join(path,listFile),'r') as f:
        lines = f.readlines()
    lines = [line.strip() for line in lines]
    lines[0] = lines[0].replace(':','')
    #print('lines = ',lines)
    for i in range(50,0,-1):
        s = ''
        for j in range(i):
            s += ' '
        lines = [line.replace(s,' ') for line in lines]
    lines = [line.replace('\t',',') for line in lines]
    lines = [line.replace(' ',',') for line in lines]
    #print('lines = ',lines)
    with open(os.path.join(path,listFile[:listFile.rfind('.')]+'.csv'),'w') as f:
        for line in lines:
            if line.find(',') > 0:
#                if (line == lines[0]) or ((line != lines[0]) and ((not ',Type,' in lines[0]) or ((',Type,' in lines[0]) and (',P,' in line)))):
                f.write(line+'\n')

def findHashIDs(path):
    notFound = []
    if os.path.exists(notFoundFileName):
        with open(notFoundFileName,'r') as f:
            notFound = f.readlines()
    fibreListName = list(pathlib.Path(path).glob('*.csv'))[0].name
    fibreList = csvFree.readCSVFile(os.path.join(path,fibreListName))

    fitsFiles = list(pathlib.Path(path).glob('*.fits'))
    fitsFiles = [fitsFile.name for fitsFile in fitsFiles]
    targetNames = [fitsFile[:fitsFile.find('_')] for fitsFile in fitsFiles]

    searchList = []
    if os.path.exists(hashSearchFileName):
        with open(hashSearchFileName,'r') as f:
            searchList = f.readlines()
    for targetName in targetNames:
        found = False
        for i in range(fibreList.size()):
            if targetName == 'PNG357.12+01.66':
                print('searching for PNG357.12+01.66')
                if 'PNG357.12' in fibreList.getData('Name',i):
                    print('PNG357.12 found in ',fibreList.getData('Name',i))
            if fibreList.getData('Name',i).replace('PNG','') == targetName.replace('PNG',''):
                found = True
                searchList.append('%s,%s,%s\n' % (targetName,
                                                  degToHMS(np.degrees(float(fibreList.getData('RA',i)))),
                                                  degToDMS(np.degrees(float(fibreList.getData('DEC',i)))),
                                                  ))
        if not found:
            print('ERROR: did not find target <'+targetName+'> in '+fibreListName)
            if targetName == 'PNG357.12+01.66':
                print('PNG357.12+01.66 not found')
                STOP
            notFound.append(targetName+' in '+fibreListName+'\n')

    """ double check that all files mentioned in the fibre list are there"""
    for i in range(fibreList.size()):
        found = False
        for targetName in targetNames:
            if fibreList.getData('Name',i).replace('PNG','') == targetName.replace('PNG',''):
                found = True
        if not found:
            if 'Parked' not in fibreList.getData('Name',i):
                if ('Type' not in fibreList.header) or (('Type' in fibreList.header) and (fibreList.getData('Type',i) == 'P')):
                    notFound.append(fibreList.getData('Name',i)+'_*.fits in '+path[path.rfind('/')+1:]+'\n')

    with open(hashSearchFileName,'w') as f:
        for line in searchList:
            f.write(line)
    with open(notFoundFileName,'w') as f:
        for line in notFound:
            f.write(line)

def map2list(mapFileName,listFileName):
    lines = []
    with open(mapFileName,'r') as f:
        lines = f.readlines()
    lines = [line.replace(',','\t') for line in lines]
    with open(listFileName,'w') as f:
        f.write('Fibre: Name: RA: DEC:\n')
        for line in lines:
            while line.count('\t') > 3:
                line = line[:line.rfind('\t')]
            print('line = '+line+': line.count(+) = ',line.count('+'))
            if line[:3] == '127':
                print('line.find(PNG357.12+01.66) = ',line.find('PNG357.12+01.66'))
            f.write(line+'\n')

def findClosestMatch(hashFoundFile,outFileName):
    hashFound = csvFree.readCSVFile(hashFoundFile)
    closestMatch = csvData.CSVData()
    closestMatch.header = hashFound.header

    hashFoundNames = np.array([str(n) for n in hashFound.getData('id')])
    print('hashFoundNames = ',hashFoundNames)
    nDoubles = 0
    for name in hashFoundNames:
        name = str(name)
        idx = np.where(hashFoundNames == name)[0]
        if len(idx) == 0:
            arr = bytes(name, 'utf-8')
            arr2 = bytes(hashFoundNames[0], 'utf-8')
            for byte in arr:
                print(byte, end=' ')
            print("\n")
            for byte in arr2:
                print(byte, end=' ')
            print("\n")
            if name == hashFoundNames[0]:
                print('both are the same')
            else:
                print('both are still not the same')
            if str(name) == str(hashFoundNames[0]):
                print('both are the same if using unicode')
            else:
                print('both are still not the same if using unicode')

        print('found name <'+name+'> ',len(idx),' times in indices ',idx)
        lineToAdd = ''
        if len(idx) == 1:
            lineToAdd = hashFound.getData(idx[0])
        else:
            dists = []
            for i in range(len(idx)):
                dist = hashFound.getData('dist[arcsec]',idx[i])
                if dist != '':
                    dists.append(float(dist))
            if len(dists) > 0:
                minId = np.where(dists == np.min(dists))[0][0]
            else:
                minId = 0
            print('idx[',minId,'] = ',idx[minId])
            lineToAdd = hashFound.getData(idx[minId])
        closestMatchNames = closestMatch.getData('id')
        if name not in closestMatchNames:
            print('lineToAdd = ',lineToAdd)
            closestMatch.append(lineToAdd)
        else:
            print('name <'+name+'> already in closestMatch')
            nDoubles += 1
    print('nDoubles = ',nDoubles)
    csvFree.writeCSVFile(closestMatch,outFileName)

def checkDistances(inFileName, angDiamFileName, outFileName):
    hashSearch = csvFree.readCSVFile(hashSearchFileName[:hashSearchFileName.rfind('.')]+'_with_header.csv')
    hashFound = csvFree.readCSVFile(inFileName)
    angDiams = csvFree.readCSVFile(angDiamFileName)
    csvOut = csvData.CSVData()
    csvOut.header = hashFound.header
    nNotFound = 0
    ra_3238 = hmsToDeg('17:59:45.20')
    dec_3238 = dmsToDeg('-33:21:13.00')
    toDelete = []
    for i in range(hashFound.size()):
        name = hashFound.getData('id',i)
        idPNMain = hashFound.getData('pndb',i)
        if idPNMain == '':
            csvOut.append(hashFound.getData(i))
        else:
            dist = float(hashFound.getData('dist[arcsec]',i))
            found = False
            for j in range(angDiams.size()):
                if angDiams.getData('idPNMain',j) == idPNMain:
                    if angDiams.getData('InUse',j) == '1':
                        found = True
                        if dist > float(angDiams.getData('MajDiam',j)):
                            csvOut.append([name,'',''])
                        else:
                            csvOut.append(hashFound.getData(i))
            if not found:
                nNotFound += 1
                print('Problem: did not find an angular diameter for <'+name+'>: idPNMain = ',idPNMain,', dist = ',dist)
                if dist > 50.:
                    csvOut.append([name,'',''])
                else:
                    csvOut.append(hashFound.getData(i))
        for j in range(hashSearch.size()):
            if hashSearch.getData('id',j) == csvOut.getData('id',csvOut.size()-1):
                ra = hmsToDeg(hashSearch.getData('ra',j))
                dec = dmsToDeg(hashSearch.getData('dec',j))
                angDist = angularDistance(ra,dec,ra_3238,dec_3238)
                print('ra = ',ra,', dec = ',dec,': angDist = ',angDist)
                angDistPyAsl = pyasl.getAngDist(ra, dec, ra_3238, dec_3238) * 3600.
                if angDist < 800:
                    if csvOut.getData('pndb',csvOut.size()-1) != '3238':
                        toDelete.append([csvOut.getData('pndb',csvOut.size()-1), ra, dec, angDist, angDistPyAsl])
                        csvOut.setData('pndb',csvOut.size()-1,'3238')
    csvFree.writeCSVFile(csvOut, outFileName)
    for i in toDelete:
        print('toDelete : ',i)
    #STOP

def addNewCandidates(inFileName, sqlFileName):
    ids = ''
    idsNew = ''
    nIds = 0
    idPNMainStart = 0
    idtbCNamesStart = 0
    newNames = []
    rows = csv.DictReader(open(pnMainFile))
    pngNames = []
    for row in rows:
        pngNames.append(row['PNG'])
        idPNMainStart = int(row['idPNMain'])+1
    rows = csv.DictReader(open(tbCNamesFile))
    for row in rows:
        idtbCNamesStart = int(row['idtbCNames'])+1
    print('idPNMainStart = ',idPNMainStart)
    print('idtbCNamesStart = ',idtbCNamesStart)
    nNew = 0
    nFound = 0
    with open(sqlFileName,'w') as f:
        rows = csv.DictReader(open(inFileName))
        for row in rows:
            if row['pndb'] == '':
                nNew += 1
                print(row['id']+": new PN")
                rowsRaDec = csv.DictReader(open(raDecFileName))
                ra = ''
                dec = ''
                for rowRaDec in rowsRaDec:
                    if row['id'] == rowRaDec['id']:
                        ra = rowRaDec['RA']
                        dec = rowRaDec['DEC']
                lon,lat = raDecToLonLat(hmsToDeg(ra),dmsToDeg(dec))
                pngName = getNewPNGName(pngNames,lon,lat)
                print('pngName = <'+pngName+'>')
                pngNames.append(pngName)
                if True:
                    newNames.append(row['id'])
                    f.write("USE `MainGPN`;\n")
                    f.write("INSERT INTO `PNMain`(`idPNMain`,`PNG`,`refPNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDECJ2000`,`Glon`,`Glat`,`refCoord`,`Catalogue`,`refCatalogue`,`userRecord`,`domain`,`refDomain`,`PNstat`,`refPNstat`,`refSimbadID`,`show`)")
                    f.write("VALUES (%d,'%s','%s','%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n"
                            % (idPNMainStart,
                            pngName,
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
                            'y'
                            ))
                    f.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`)")
                    f.write("VALUES (%d,'%s','%s',%d,'%s','%s',%d);\n"
                            % (idtbCNamesStart,
                            row['id'],
                            'ziggy',
                            1,
                            'ziggy',
                            'ziggy',
                            idPNMainStart
                            ))
                    f.write("INSERT INTO `PNMain_tbCNames`(`idPNMain`,`idtbCNames`)")
                    f.write("VALUES (%d,%d);\n"
                            % (idPNMainStart,
                            idtbCNamesStart
                            ))
                    idtbCNamesStart += 1
                idPNMainStart += 1
                idsNew += str(idPNMainStart)+','
                nIds += 1
            else:
                nFound += 1
                print(row['id']+": "+row['pndb'])
                nameFound = False
                tbCNames = csv.DictReader(open(tbCNamesFile))
                for tbCName in tbCNames:
                    if tbCName['Name'].replace(' ','') == row['id'].replace(' ',''):
                        nameFound = True
                if not nameFound:
                    if row['id'][:2] != 'PN':
                        f.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`)")
                        f.write("VALUES (%d,'%s','%s',%d,'%s','%s',%d);\n"
                                % (idtbCNamesStart,
                                row['id'],
                                'ziggy',
                                0,
                                'ziggy',
                                'ziggy',
                                int(row['pndb']),
                                ))
#                        f.write("INSERT INTO `PNMain_tbCNames`(`idtbCNames`,`idPNMain`)")
#                        f.write("VALUES (%d,%d);\n"
#                                % (idtbCNamesStart,
#                                int(row['pndb']),
#                                ))
                        idtbCNamesStart += 1
                ids += row['pndb']+','

    print('ids = ',nIds,': ',ids)
    print('idsNew = ',idsNew,' newNames = ',newNames)
    print('nNew = ',nNew,', nFound = ',nFound)

def getIds(filename):
    lis = csvFree.readCSVFile(filename)
    ids = lis.getData('pndb')
    idsOut = []
    for id in ids:
        if id not in idsOut:
            idsOut.append(id)
    for id in idsOut:
        print(id, end=',')
    print('\n')

if __name__ == '__main__':
    hashFoundFile = os.path.join(spectraPath,'hash_found_all.csv')
    hashFoundClosestFile = hashFoundFile[:hashFoundFile.rfind('.')]+'_closest.csv'
    hashFoundDistCheckFile = hashFoundClosestFile[:hashFoundClosestFile.rfind('.')]+'_distances_checked.csv'
    #map2list(os.path.join(spectraPath,'N2/map.txt'),os.path.join(spectraPath,'N2/map.list'))
    if False:
        if os.path.exists(notFoundFileName):
            os.remove(notFoundFileName)
        if os.path.exists(hashSearchFileName):
            os.remove(hashSearchFileName)
        for sub in subPaths:
            list2csv(os.path.join(spectraPath,sub))
            findHashIDs(os.path.join(spectraPath,sub))
    if False:
        findClosestMatch(hashFoundFile,hashFoundClosestFile)
        checkDistances(hashFoundClosestFile,hashAngDiamFileName,hashFoundDistCheckFile)
        addNewCandidates(hashFoundDistCheckFile, sqlFileOut)
    getIds(hashFoundDistCheckFile)
