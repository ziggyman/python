import numpy as np
import os

def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

def readFrenchAmateursCatalogue(fName='/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/FrenchAmateursCatalogue.txt'):
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

#frenchAmateursCatalogue = readFrenchAmateursCatalogue()
leDuCatalogue = readLeDu()
