import numpy as np
from PyAstronomy import pyasl
from myUtils import lonLatToRaDec,getHeader
import csvFree
import csvData

def fixReesTable(fNameIn, fNameOut=None):
    newLines = []
    with open(fNameIn, 'r') as f:
        lines = f.readlines()
        print('lines = ',lines)
        for line in lines:
            line = line.replace(' ± ','±')
            line = line.replace(' ±','±')
            lineA = ''
            lineB = ''
            if line[0:len(line)-5].find('NTT') > 0:
                str = 'NTT'
            else:
                str = 'HST'
            pos = line.find(str)
            lineA = line[0:pos+3]+'\n'
            lineB = line[pos+4:]
            newLines.append(lineA.replace(' ',','))
            newLines.append(lineB.replace(' ',','))
    if fNameOut:
        with open(fNameOut,'w') as f:
            f.write('PNG,Morphology,PA,GPA,Telescope')
            for line in newLines:
                f.write(line)

def fixWeidmannTable(fNameIn, fNameOut=None):
    newLines = []
    with open(fNameIn, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.replace('−','-')
            if line.find(',') > 0:
                newLine = line
            else:
                keys = line.split(' ')
                for iKey in np.arange(0,len(keys),1):
                    if keys[iKey][0] != '.':
                        if iKey == 0:
                            newLine = keys[iKey]
                            if keys[iKey + 1][0] != '.':
                                newLine += ' '+keys[iKey+1]
                            iKey += 1
                        elif iKey > 1:
                            newLine += ','+keys[iKey]
            newLines.append(newLine)
    if fNameOut:
        with open(fNameOut, 'w') as f:
            for line in newLines:
                f.write(line)

def crossMatchRees(fNameInA, fNameInB, fNameOutMatch=None, fNameOutNoMatch=None):
    with open(fNameInA,'r') as fA:
        linesA = fA.readlines()
    with open(fNameInB,'r') as fB:
        linesB = fB.readlines()
    match = []
    noMatch = []
    for lineA in linesA:
        lineATemp = lineA[lineA.find(',')+1:]
        pngNameA = lineATemp[:lineATemp.find(',')]
        print('pngNameA = <'+pngNameA+'>')
        found = False
        newLine = lineA
        for lineB in linesB:
            pngNameB = lineB[:lineB.find(',')]
#            print('pngNameB = <'+pngNameB+'>')
            if pngNameA == pngNameB:
                found = True
                newLine = newLine[:newLine.find('\n')] + ',' + lineB
        if found:
            match.append(newLine)
        else:
            noMatch.append(newLine)
    if fNameOutMatch:
        with open(fNameOutMatch,'w') as f:
#            f.write(linesA[0][:linesA[0].find('\n')]+','+linesB[0])
            for line in match:
                f.write(line)
    if fNameOutNoMatch:
        with open(fNameOutNoMatch,'w') as f:
            f.write(linesA[0])
            for line in noMatch:
                f.write(line)

def crossMatchWeidmann(fNameInW, fNameInHASH, fNameOutMatch=None, fNameOutNoMatch=None):
    with open(fNameInW, 'r') as f:
        linesW = f.readlines()
    with open(fNameInHASH, 'r') as f:
        linesHASH = f.readlines()
    for iLine in np.arange(0,len(linesHASH),1):
        linesHASH[iLine] = linesHASH[iLine].rstrip('\n').split(',')
    match = []
    noMatch = []
    headerW = []
    nFound = 0
    for lineW in linesW:
        lineW = lineW.rstrip('\n')
        itemsW = lineW.split(',')
        if len(headerW) == 0:
            headerW = itemsW
        else:
            nameW = itemsW[0]
#            print('read name <'+nameW+'> from Weidmann file')
            found = False
            for iLine in np.arange(0,len(linesHASH),1):
#                print('nameHASH = <'+linesHASH[iLine][2].strip('"')+'>')
                if (linesHASH[iLine][2].strip('"') == nameW) or ('G'+linesHASH[iLine][1] == nameW):
                    print('nameW = <'+nameW+'> found in HASH as '+linesHASH[iLine][2].strip('"'))
                    found = True
                    nFound += 1
            if not found:
                l = float(itemsW[1])
                b = float(itemsW[2])
                ra1,dec1 = lonLatToRaDec(l, b)
                for iLine in np.arange(1,len(linesHASH),1):
                    ra2 = float(linesHASH[iLine][6])
                    dec2 = float(linesHASH[iLine][7])
                    dist = pyasl.getAngDist(ra1, dec1, ra2, dec2) * 3600.
#                    print('dist = ',dist)
                    if dist < 5.:
                        found = True
                        nFound += 1
                        print('nameW = <'+nameW+'> identified as ',linesHASH[iLine][2])
            if not found:
                print('nameW = <'+nameW+'> not found in HASH: l=',itemsW[1],', b=',itemsW[2])

#                print('nameW = <'+nameW+'> NOT found in HASH')
    print('found ',nFound,' names in both catalogues out of ',len(linesW),' PN in Weidmann table')

def crossmatch_hash_old_and_new():
    hash_old_filename = '/Users/azuri/daten/uni/HKU/PN_orientation_angles/HASH_bipolar+elliptical_true+likely_PNe+004.2-05.9+005.9-02.6+001.4+05.3+357.5+03.2+357.6-03.3.csv'
    hash_new_filename = '/Users/azuri/daten/uni/HKU/PN_orientation_angles/HASH_bipolar+elliptical_true+likely_PNe+004.2-05.9+005.9-02.6+001.4+05.3+357.5+03.2+357.6-03.3_new.csv'
    hash_old = csvFree.readCSVFile(hash_old_filename)
    hash_new = csvFree.readCSVFile(hash_new_filename)

    hash_diff = csvData.CSVData()
    hash_diff.header = hash_new.header
    for i in range(hash_new.size()):
        idPNMain = hash_new.getData('idPNMain',i)
        if hash_old.find('idPNMain',idPNMain)[0] < 0:
            hash_diff.append(hash_new.getData(i))
    print('hash_diff.size() = ',hash_diff.size())
    csvFree.writeCSVFile(hash_diff,'/Users/azuri/daten/uni/HKU/PN_orientation_angles/HASH_bipolar+elliptical_true+likely_PNe+004.2-05.9+005.9-02.6+001.4+05.3+357.5+03.2+357.6-03.3_diff.csv')

    hash_remove = csvData.CSVData()
    hash_remove.header = hash_old.header
    for i in range(hash_old.size()):
        idPNMain = hash_old.getData('idPNMain',i)
        if hash_new.find('idPNMain',idPNMain)[0] < 0:
            hash_remove.append(hash_old.getData(i))
    print('hash_remove.size() = ',hash_remove.size())
    csvFree.writeCSVFile(hash_remove,'/Users/azuri/daten/uni/HKU/PN_orientation_angles/HASH_bipolar+elliptical_true+likely_PNe+004.2-05.9+005.9-02.6+001.4+05.3+357.5+03.2+357.6-03.3_remove.csv')

def identifyNTTobjects():
    fitslistFileName = '/Users/azuri/daten/uni/HKU/PN_orientation_angles/NTT/fits.list'
    with open(fitslistFileName,'r') as f:
        lines = f.readlines()
    hashFull = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/PN_orientation_angles/hash_PNMain_full_260922.csv')
    with open(fitslistFileName[:fitslistFileName.rfind('.')]+'_idPN.text','w') as f:
        for line in lines:
            header = getHeader(line.strip(),0)
            object = header['OBJECT']
            if object == 'PN003.3-1.6':
                object = 'PN003.3-01.6'
            if 'PN' in object:
                idx = hashFull.find('PNG',object.replace('PN','').lower())
                if (len(idx) > 1) or (idx[0] == -1):
                    idx = hashFull.find('PNG',object.replace('PN','').replace('A',''))
                    if (len(idx) > 1) or (idx[0] == -1):
                        print(line.strip()+': header[OBJECT] = ',object)
                        print('idx = ',idx)
                idPNMain = hashFull.getData('idPNMain',idx[0])
                f.write(line.strip()+' '+object+' '+idPNMain+'\n')

#identifyNTTobjects()
#STOP
crossmatch_hash_old_and_new()

#fixReesTable('/Users/azuri/daten/uni/HKU/PN alignment/Rees_Zijlstra_table.txt','/Users/azuri/daten/uni/HKU/PN alignment/Rees_Zijlstra_table.csv')
#fixWeidmannTable('/Users/azuri/daten/uni/HKU/PN alignment/Weidmann+Diaz_table.txt','/Users/azuri/daten/uni/HKU/PN alignment/Weidmann+Diaz_table.csv')
#crossMatchRees('/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe.csv','/Users/azuri/daten/uni/HKU/PN alignment/Rees_Zijlstra_table.csv','/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe+Rees_Zijlstra.csv','/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe-Rees_Zijlstra.csv')
#crossMatchWeidmann('/Users/azuri/daten/uni/HKU/PN alignment/Weidmann+Diaz_table.csv', '/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe.csv','/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe+Weidmann.csv','/Users/azuri/daten/uni/HKU/PN alignment/HASH_bipolar+elliptical_true_PNe-Weidmann.csv')
