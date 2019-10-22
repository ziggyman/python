import numpy as np
import os
from astropy.io import fits

def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

def readHASHSelection(fName='/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/FrenchAmateursCatalogue.txt'):
    lines = readFileToArr(fName)
    lines[1] = lines[1].split(' ')
    catalogueKeys = [lines[0]]
    for key in lines[1]:
        catalogueKeys.append(key)
    print('catalogueKeys = ',catalogueKeys)
    catLines = []
    for iLine in np.arange(2,len(lines),2):
        strLine = lines[iLine+1].split('\t')
        print('iLine = ',iLine,': strLine = ',strLine)
        catLine = {catalogueKeys[0]:lines[iLine],catalogueKeys[1]:strLine[0],catalogueKeys[2]:strLine[1],catalogueKeys[3]:strLine[2]}
        print('catLine = ',catLine)
        catLines.append(catLine)
    print('lines = ',lines)
    print('lines[0] = ',lines[0])
    print('lines[1] = ',lines[1])
    print('catLines = ',catLines)
    return catLines

def readLeDu(fName='/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/LeDu2017_table2.csv'):
    lines = readFileToArr(fName)
    keys = [key.strip("'") for key in lines[0].split(',')]
    catLines = []
    print('keys = ',len(keys),': ',keys)
    for lineI in lines[1:]:
        line = []
        lineII = lineI[1:]
        while lineII.find("'") >= 0:
            print('lineII = ',lineII)
            if lineII[:3] == 'ULL':
                print('lineII[:3] = ',lineII[:3],' == ULL')
                val = 'NULL'
                if len(lineII) > 3:
                    lineII = lineII[len(val)+1:]
            else:
                val = lineII[:lineII.find("'")]
                if len(lineII) > 3:
                    lineII = lineII[len(val) + 3:]
            print('val = ',val)
            line.append(val)
        if lineII[:3] == 'ULL':
            line.append('NULL')
        if len(line) != 0:
            print('line = ',len(line),': ',line)
            if len(keys) != len(line):
                print('ERROR: keys = ',len(keys),': ',keys)
                print('ERROR: lineI = ',len(lineI),': ',lineI)
                print('ERROR: line = ',len(line),': ',line)
                STOP
            myDict = {}
            for iKey in np.arange(0,len(keys),1):
                myDict.update({keys[iKey]:line[iKey]})
            print('myDict = ',myDict)
            catLines.append(myDict)
    print('catLines = ',len(catLines),': ',catLines)
    return catLines

#frenchAmateursCatalogue = readHASHSelection()
#leDuCatalogue = readLeDu()
#LAstrCatalogue = readHASHSelection('/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/2017LAstr131b46L.txt')

def renameFRAFits(fName):
    lines = readFileToArr(fName)
    linesOut = []
    for line in lines:
        fitsNameIn = line[line.rfind('/')+1:]
        fitsNameInTemp = fitsNameIn
        fitsNameOut = line[0:line.rfind('/ENVO')+1]+'all/'
        if fitsNameIn[0] == '_':
            fitsNameInTemp = fitsNameIn[1:]
        fitsNameOut += fitsNameInTemp[0:fitsNameInTemp.find('_')+1].upper()+'HP'
#            date = fitsNameInTemp[fitsNameInTemp.find('_')+1:fitsNameInTemp.find('_')+9]
        date = fitsNameInTemp[fitsNameInTemp.find('_20')+1:fitsNameInTemp.find('_20')+9]
        fitsNameOut += date[6:]+date[4:6]+date[2:4]+'.fits'
        linesOut.append(fitsNameOut)
        print('fitsNameIn = <'+fitsNameIn+'>: date = <'+date+'>: fitsNameOut = <'+fitsNameOut+'>')
        hdu_list = fits.open(line)
        hdu_list.info()
        image_data = fits.getdata(line)
        print(image_data.shape)

renameFRAFits('/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/allFits.list')

