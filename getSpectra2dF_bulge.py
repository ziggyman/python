import numpy as np
import csvFree,csvData
from myUtils import raDecToLonLat

csvFitsFiles = '/Users/azuri/daten/uni/HKU/HASH/hash_FitsFiles_181121.csv'

fitsFilesTable = csvFree.readCSVFile(csvFitsFiles)

setnames = fitsFilesTable.getData('setname')
outsideBulge = []
for i in range(len(setnames)):
    if setnames[i] == 'AAOmega_bulge':
        dra = float(fitsFilesTable.getData('DRAJ2000',i))
        ddec = float(fitsFilesTable.getData('DDECJ2000',i))
        l,b = raDecToLonLat(dra,ddec)
        if abs(l) > 10. or abs(b) > 10.:
            outsideBulge.append(i)

print('outsideBulge = ',len(outsideBulge),': ',outsideBulge)
STOP

instrument = fitsFilesTable.getData('instrument')

for i in np.arange(len(instrument)-1,-1,-1):
    inst = instrument[i]
    if '2dF' not in inst:
        fitsFilesTable.removeRow(i)

print('found ',fitsFilesTable.size(),' spectra from 2dF')
instrument = fitsFilesTable.getData('instrument')
print('instrument = ',instrument)

ras = [float(x) for x in fitsFilesTable.getData('DRAJ2000')]
decs = [float(x) for x in fitsFilesTable.getData('DDECJ2000')]

lsAndBs = []
for i in range(len(ras)):
    lsAndBs.append(raDecToLonLat(ras[i],decs[i]))

print('lsAndBs = ',lsAndBs)

for i in np.arange(len(lsAndBs)-1,-1,-1):
    l,b = [lsAndBs[i][0],lsAndBs[i][1]]
    if abs(l) > 10. or abs(b) > 10.:
        fitsFilesTable.removeRow(i)

print('setNames = ',fitsFilesTable.getData('setname'))

    
