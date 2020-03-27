import os
import numpy as np
import difflib

import hammer as ham
import csvFree, csvData

path = '/Users/azuri/daten/uni/HKU/PN alignment'
rzFile = os.path.join(path, 'Rees_Zijlstra_table.csv')
myDataFile = os.path.join(path, 'PN-alignments.csv')
#myHashFile = os.path.join(path, 'HASH_bipolar+elliptical_true_PNe.csv')
myHashFile = os.path.join(path, 'HASH_true_PNe+004.2-05.9+005.9-02.6.csv')
dataFileOut = os.path.join(path, 'Rees_Zijlstra_table_with_HASH-ID.csv')

csvOut = csvData.CSVData()
csvOut.header = ['HASH ID','EPA','flag','source','GPA','DRAJ2000','DDECJ2000','l','b','GPA2','comment']

rzData = csvFree.readCSVFile(rzFile)
myData = csvFree.readCSVFile(myDataFile)
myHashData = csvFree.readCSVFile(myHashFile)

rzPNGs = rzData.getData('PNG')
rzMorph = rzData.getData('Morphology')
myPNGs = myHashData.getData('PNG')
print('rzPNGs = ',rzPNGs)
print('myPNGs = ',myPNGs)

print('repr(000.2-01.9) = <'+repr('000.2-01.9')+'>')
print('repr(000.2−01.9) = <'+repr('000.2−01.9')+'>')
print('difflib.ndiff(000.2-01.9,000.2−01.9) = ',''.join(difflib.ndiff('000.2-01.9','000.2−01.9')))
for line in difflib.ndiff('000.2-01.9','000.2−01.9'):
    print(line)
if '000.2-01.9' in rzPNGs:
    print('000.2-01.9 is in rzPNGs')
else:
    print('000.2-01.9 is not in rzPNGs')

if '000.2−01.9' in rzPNGs:
    print('000.2−01.9 is in rzPNGs')
else:
    print('000.2−01.9 is not in rzPNGs')

if '000.2-01.9' in myPNGs:
    print('000.2-01.9 is there')
else:
    print('000.2-01.9 is NOT there')
#STOP
nFound = 0
nNotFound = 0
for i in np.arange(0,len(rzPNGs),1):
    lineOut = []
    p = rzPNGs[i]
    found = False
    if rzMorph[i] in ['Bipolar']:
        for j in np.arange(0,len(myPNGs),1):
            if p == myPNGs[j]:
                print('PNG '+p+' found in my data')
                found = True
                nFound += 1
                lineOut.append(myHashData.getData('idPNMain',j))
                lineOut.append(rzData.getData('PA',i))
                lineOut.append('1')
                lineOut.append(rzData.getData('Telescope',i))
                lineOut.append(rzData.getData('GPA',i).split('±')[0])
                lineOut.append(myHashData.getData('DRAJ2000',j))
                lineOut.append(myHashData.getData('DDECJ2000',j))
                lineOut.append(myHashData.getData('Glon',j))
                lineOut.append(myHashData.getData('Glat',j))
                lineOut.append(rzData.getData('GPA',i))
                lineOut.append('')
                csvOut.append(lineOut)


        if not found:
            print('PNG '+p+' not found in my data')
            nNotFound += 1
print('found ',nFound,' PNe in my data')
print('did not find ',nNotFound,' PNe in my data')

csvFree.writeCSVFile(csvOut,dataFileOut)
