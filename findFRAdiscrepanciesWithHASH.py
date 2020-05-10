import numpy as np
from csvFree import readCSVFile,writeCSVFile
import csvData
from myUtils import dmsToDeg, hmsToDeg, angularDistancePyAsl, degToArcsec

allHASHobjectsFile = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/HASH_all_objects.csv'
allFRAobjectsFile = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/FRA_all_with_header.csv'
hashOutFile = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/FRAall_HASHoutput.csv'
hashCommonNamesFile = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/hash-commonNames.csv'
discrepanciesAFile = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/TableI I_ 120objets_25032020_discrepancies.csv'
discrepanciesBFile = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/TableI_568PNG_25032020_discrepancies.csv'

disrepanciesOutFileName = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/April2020/discrepancies.csv'

maxAngularDistance = 15.

def combineDiscrepancies():
    discrepanciesA = readCSVFile(discrepanciesAFile)
    discrepanciesB = readCSVFile(discrepanciesBFile)

    csvOut = csvData.CSVData()
    csvOut.header = ['Name','RA','DEC']

    for c in [discrepanciesA,discrepanciesB]:
        for i in range(c.size()):
            lineOut = [c.getData('NOM',i),c.getData('AD:(J2000)',i),c.getData('DEC (J2000)',i)]
            csvOut.append(lineOut)
    print('csvOut.size() = ',csvOut.size())
    return csvOut

def getName(csv,keyWord,line):
    name = csv.getData(keyWord,line)
    if name[0:2] == 'LD':
        name = 'LDu'+name[2:]
    name = name.replace(' ','')
    return name

def getCommonNames(csv,hashIDKeyWord,hashID):
    hashIDs = csv.getData(hashIDKeyWord)
    commonNames = []
    for i in range(len(hashIDs)):
        if hashIDs[i] == hashID:
            commonNames.append(csv.getData('Name',i).replace(' ',''))
    return commonNames

def checkCombinedDiscrepancies(combinedDiscrepancies, allHASHobjects, hashCommonNames):
    csvOut = csvData.CSVData()
    csvOut.header = ['Name','idPNMain','HASH common names', 'RA FRA', 'RA HASH', 'DEC FRA', 'DEC HASH', 'angular distance [arcsec]']
    csvOut.data = []
    emptyData = ['','','', '', '', '', '', '']
    for i in range(combinedDiscrepancies.size()):
        ra = combinedDiscrepancies.getData('RA',i)
        dec = combinedDiscrepancies.getData('DEC',i)
        raDeg = hmsToDeg(ra)
        decDeg = dmsToDeg(dec)
        name = getName(combinedDiscrepancies,'Name',i)
        print('name = <'+name+'>')
        found = False
        for j in range(allHASHobjects.size()):
            if name == allHASHobjects.getData('Name',j).replace(' ',''):
                hashID = allHASHobjects.getData('idPNMain',j)
                angDist = degToArcsec(angularDistancePyAsl(raDeg,decDeg,hmsToDeg(allHASHobjects.getData('RAJ2000',j)),dmsToDeg(allHASHobjects.getData('DECJ2000',j))))
                commonNames = getCommonNames(hashCommonNames,'idPNMain',hashID)
                print('Name = '+name+': HASH ID = ',hashID,': commonNames = ',commonNames,': ra = ',ra,', RA = ',allHASHobjects.getData('RAJ2000',j),', dec = ',dec,', DEC = ',allHASHobjects.getData('DECJ2000',j),', angDist = ',angDist)
                found = True
        if not found:
            print('ERROR: object with name <'+name+'> not found in HASH')
            for j in range(allHASHobjects.size()):
                angDist = degToArcsec(angularDistancePyAsl(raDeg,decDeg,hmsToDeg(allHASHobjects.getData('RAJ2000',j)),dmsToDeg(allHASHobjects.getData('DECJ2000',j))))
                if angDist < maxAngularDistance:
                    hashID = allHASHobjects.getData('idPNMain',j)
                    commonNames = getCommonNames(hashCommonNames,'idPNMain',hashID)
                    if name in commonNames:
                        print('Name = '+name+': HASH ID = ',hashID,': commonNames = ',commonNames,': ra = ',ra,', RA = ',allHASHobjects.getData('RAJ2000',j),', dec = ',dec,', DEC = ',allHASHobjects.getData('DECJ2000',j),', angDist = ',angDist)
                        found = True
                    else:
                        print('Name = <'+name+'>: object found within ',maxAngularDistance,' arcsec: angDist = ',angDist,': ra = ',ra,', RA = ',allHASHobjects.getData('RAJ2000',j),', dec = ',dec,', DEC = ',allHASHobjects.getData('DECJ2000',j),', HASH name = <'+allHASHobjects.getData('Name',j))
                        csvOut.append(emptyData)
                        csvOut.setData('Name',csvOut.size()-1,name)
                        csvOut.setData('idPNMain',csvOut.size()-1,hashID)
                        cNames = commonNames[0]
                        for k in np.arange(1,len(commonNames),1):
                            cNames += ';'+commonNames[k]
                        csvOut.setData('HASH common names',csvOut.size()-1,cNames)
                        csvOut.setData('RA FRA',csvOut.size()-1,ra)
                        csvOut.setData('RA HASH',csvOut.size()-1,allHASHobjects.getData('RAJ2000',j))
                        csvOut.setData('DEC FRA',csvOut.size()-1,dec)
                        csvOut.setData('DEC HASH',csvOut.size()-1,allHASHobjects.getData('DECJ2000',j))
                        csvOut.setData('angular distance [arcsec]',csvOut.size()-1,str(angDist))
    writeCSVFile(csvOut,disrepanciesOutFileName)

if __name__ == "__main__":
    allHASHobjects = readCSVFile(allHASHobjectsFile)
    allFRAobjects = readCSVFile(allFRAobjectsFile)
    hashOut = readCSVFile(hashOutFile)
    hashCommonNames = readCSVFile(hashCommonNamesFile)

    combinedDiscrepancies = combineDiscrepancies()
    checkCombinedDiscrepancies(combinedDiscrepancies, allHASHobjects, hashCommonNames)

