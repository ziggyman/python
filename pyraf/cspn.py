import csv
import numpy as np
from matplotlib import pyplot as plt
import csvFree,csvData
from myUtils import angularDistancePyAsl,hmsToDeg,dmsToDeg,raDecToLonLat

inputSGFile = '/Users/azuri/daten/uni/HKU/publications/CSPN/gonzalez-santamaria2021_tablea1.dat'
inputHashFile = '/Users/azuri/daten/uni/HKU/publications/CSPN/hash_tbCSCoords_240122.csv'
inputHashDiamFile = '/Users/azuri/daten/uni/HKU/publications/CSPN/hash_tbAngDiam_240122.csv'
inputHashPNMainFile = '/Users/azuri/daten/uni/HKU/publications/CSPN/hash_PNMain_240122.csv'

csvDiam = csvFree.readCSVFile(inputHashDiamFile)

outputSQLFile = '/Users/azuri/daten/uni/HKU/publications/CSPN/hash_tbCSCoords_240122.sql'

def fixHashFile():
    with open(outputSQLFile,'w') as w:
        csvHash = csvData.CSVData()
        with open(inputHashFile,'r') as f:
            hashData = csv.DictReader(f)
            print('hashData.fieldnames = ',hashData.fieldnames)
            csvHash.header = hashData.fieldnames
        with open(inputHashFile,'r') as f:
            hashData = csv.DictReader(f)
            print('hashData.fieldnames = ',hashData.fieldnames)
            nRows = 0
            for row in hashData:
                nRows += 1
                csvHash.append([row[x] for x in hashData.fieldnames])
#                if row['CS_DRAJ2000'] == 'NULL':
                print("row['CS_DECJ2000'].split(':')[0][:2] = <"+row['CS_DECJ2000'].split(':')[0][:2]+">")
                if row['CS_DECJ2000'].split(':')[0][:2] == '-0':
                    print('row = ',row)
                    dra = hmsToDeg(csvHash.getData('CS_RAJ2000',csvHash.size()-1))
                    ddec = dmsToDeg(csvHash.getData('CS_DECJ2000',csvHash.size()-1))
                    lon, lat = raDecToLonLat(dra,ddec)
                    csvHash.setData('CS_DRAJ2000',
                                    csvHash.size()-1,
                                    str(dra))
                    csvHash.setData('CS_DDECJ2000',
                                    csvHash.size()-1,
                                    str(ddec))
                    csvHash.setData('CS_Glon',
                                    csvHash.size()-1,
                                    str(lon))
                    csvHash.setData('CS_Glat',
                                    csvHash.size()-1,
                                    str(lat))
                    w.write("UPDATE MainGPN.tbCSCoords SET CS_DRAJ2000 = "+str(dra)+", CS_DDECJ2000 = "+str(ddec)+", CS_Glon = "+str(lon)+", CS_Glat = "+str(lat)+", CSstat = 'p' WHERE idtbCSCoords = "+csvHash.getData('idtbCSCoords',csvHash.size()-1)+";\n")

    return csvHash


def plotHistogram():
    csvSG = csvData.CSVData()
    csvSG.header = ['PNG','Name','DRA','DDec','flag']
    with open(inputSGFile,'r') as f:
        linesSG = f.readlines()
    for iLine in range(len(linesSG)):
        linesSG[iLine] = linesSG[iLine].strip().split()
        print('linesSG[',iLine,'] = ',linesSG[iLine])
        csvSG.append([linesSG[iLine][1][1:],linesSG[iLine][2],linesSG[iLine][5],linesSG[iLine][6],linesSG[iLine][4]])

    csvPNMain = csvFree.readCSVFile(inputHashPNMainFile)

    csvHash = fixHashFile()
    csvHash.addColumn('PNG')
    for iRow in range(csvHash.size()):
        csvHash.setData('PNG',
                        iRow,
                        csvPNMain.getData('PNG',
                                          csvPNMain.find('idPNMain',
                                                         csvHash.getData('idPNMain',
                                                                         iRow
                                                                        )
                                                        )[0]
                                         )
                       )

    print('csvSG.size() = ',csvSG.size(),', csvHash.size() = ',csvHash.size(),', csvPNMain.size() = ',csvPNMain.size())

    print('csvPNMain.header = ',csvPNMain.header)

    distances = []
    print('csvSG.header = ',csvSG.header)
    for i in range(csvHash.size()):
        if csvHash.getData('InUse',i) == '1':
            posSG = csvSG.find('PNG',csvHash.getData('PNG',i))
            print('i = ',i,': posSG = ',posSG)
            if len(posSG) > 1:
                print('Problem: len(posSG) = ',len(posSG),' > 1')
                STOP
            if posSG[0] >= 0:
                print("csvHash.getData('CS_DRAJ2000',",i,") = ",csvHash.getData('CS_DRAJ2000',i))
                print("csvHash.getData('CS_DDECJ2000',",i,") = ",csvHash.getData('CS_DDECJ2000',i))
                print("csvSG.getData('DRA',",posSG[0],") = ",csvSG.getData('DRA',posSG[0]))
                print("csvSG.getData('DDec',",posSG[0],") = ",csvSG.getData('DDec',posSG[0]))
                distances.append([angularDistancePyAsl(float(csvHash.getData('CS_DRAJ2000',i)),
                                                    float(csvHash.getData('CS_DDECJ2000',i)),
                                                    float(csvSG.getData('DRA',posSG[0])),
                                                    float(csvSG.getData('DDec',posSG[0]))) * 3600.,
                                  csvHash.getData('idPNMain',i),
                                  csvDiam.getData('MajDiam',csvDiam.find('idPNMain',csvHash.getData('idPNMain',i))[0]),
                                  csvHash.getData('refCSstat',i),
                                  float(csvSG.getData('DRA',posSG[0])),
                                  float(csvSG.getData('DDec',posSG[0])),
                                  float(csvHash.getData('CS_DRAJ2000',i)),
                                  float(csvHash.getData('CS_DDECJ2000',i)),
                                 ]
                                )
                if distances[len(distances)-1][0] > 5:
                    print('distances[',len(distances)-1,'] = ',distances[len(distances)-1])

    distances = np.array(distances)
    print('distances = ',len(distances),': ',distances)
#    print('[d[0] for d in distances] = ',[d[0] for d in distances])
#    print('np.array( [float(d[0]) for d in distances])>5. = ',np.array( [float(d[0]) for d in distances])>5.)
    
    largeDistances = distances[np.array([float(d[0]) for d in distances]) > 5.]
    print('largeDistances = ',len(largeDistances),': ',largeDistances)
    histVals = plt.hist(np.sort(np.array([float(distance[0]) for distance in distances])),bins=40,range=[0.,40.])
    print('histVals = ',histVals)
    plt.xlabel('distance in "')
    plt.ylabel('number of CSPN')
    plt.show()

    return csvHash, csvSG

def checkForStellarPNe(csvHash, csvSG):
    stellarPNe = []
    for i in range(csvSG.size()):
        hashPos = csvHash.find('PNG',csvSG.getData('PNG',i))
        if len(hashPos) > 1:
            print('ERROR: more than one PNG found for PNG ',csvSG.getData('PNG',i))
            print('idPNMains = ',csvHash.getData('idPNMain',hashPos))
            STOP
        if hashPos[0] < 0:
            print('ERROR: PNG <'+csvSG.getData('PNG',i)+'> not found')
            STOP
        idPNMain = csvHash.getData('idPNMain',hashPos[0])
        diamPos = csvDiam.find('idPNMain',idPNMain)
        if diamPos[0] > -1:
            majDiam = csvDiam.getData('MajDiam',diamPos)
            inUse = csvDiam.getData('InUse',csvDiam.find('idPNMain',idPNMain))
            for j in range(len(inUse)):
                if inUse[j] == '1':
                    majDiam = majDiam[j]
                    print('PNG ',csvSG.getData('PNG',i),': majDiam = ',majDiam)
                    if float(majDiam) < 5.:
                        stellarPNe.append([majDiam,idPNMain])
    print('stellarPNe = ',len(stellarPNe),': ',stellarPNe)

def fixPNG(csvData):
    with open('/Users/azuri/daten/uni/HKU/publications/CSPN/fixPNG.sql','w') as f:
        for i in range(csvData.size()):
            png = csvData.getData('PNG',i)
            oldPNG = png
            print('original png = <'+png+'>')
            if png.find('+-') >= 0:
                if len(png)-png.rfind('-')-1 < 4:
                    png = png.replace('+-','-0')
                else:
                    png = png.replace('+-','-')
                    print('new png = <'+png+'>')
                    #STOP
                print('new png = <'+png+'>')
            if png.find('.') != 3:
                if png.find('.') == 1:
                    png = '00'+png
                else:
                    png = '0'+png
                print('new png = <'+png+'>')
            pmPos = png.find('+')
            if pmPos < 0:
                pmPos = png.find('-')
            if pmPos >= 0:
                prefix = png[:pmPos+1]
                latLen = png.rfind('.')-pmPos-1
                print('png = ',png,': latLen = ',latLen)
                if latLen < 2:
                    png = prefix+'0'+png[pmPos+1:]
                    print('new png = <'+png+'>')
                if pmPos < 5:
                    if pmPos < 4:
                        png = '00'+png
                        print('new png = <'+png+'>')
                        STOP
                    if pmPos < 5:
                        png = '0'+png
            if oldPNG != png:
                if csvData.find('PNG',png)[0] >= 0:
                    if csvData.find('PNG',png+'a')[0] < 0:
                        png = png+'a'
                    else:
                        if csvData.find('PNG',png+'b')[0] < 0:
                            png = png+'b'
                        else:
                            if csvData.find('PNG',png+'c')[0] < 0:
                                png = png+'c'
                            else:
                                print('problem: png <'+png+'c'+'> already exists')
                                STOP
                f.write("UPDATE MainGPN.PNMain SET PNG = '"+png+"' WHERE idPNMain = "+csvData.getData('idPNMain',i)+";\n")
                csvData.setData('PNG',i,png)
    return csvData


if __name__ == '__main__':
    #fixHashFile()
    csvHash, csvSG = plotHistogram()
    csvHash = csvFree.readCSVFile(inputHashPNMainFile)
    csvHash = fixPNG(csvHash)
    checkForStellarPNe(csvHash, csvSG)
