import os
import numpy as np
import csvFree,csvData
from myUtils import fix_UsrComments

listA = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/Master_WR_List.txt')
listB = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/WR-observing-list-ATAC10B.txt')
comments = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/hash_WR.csv')
hashPNMain = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/hash_PNMain_250624.csv')
hashCNames = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/hash_CNames_250624.csv')
fitsFiles = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/hash_FitsFiles_250624.csv')
spectraInfo = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/hash_spectraInfo_250624.csv')
fix_UsrComments('/Users/azuri/daten/uni/HKU/interns_projects/oriol/hash_userComments_250624.csv')
userComments = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/hash_userComments_250624.csv')
CSCoords = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/hash_CSCoords_250624.csv')

wrData = csvData.CSVData()
wrData.header = ['idPNMain','PNG','RA','DEC','classification','notes']

hashPNMain_idPNMains = np.array(hashPNMain.getData('idPNMain'))
print('hashPNMains_idPNMains = ',hashPNMain_idPNMains)
#for idPNMain in hashPNMain_idPNMains:
#    print(idPNMain)

for i in range(comments.size()):
    idPNMain = comments.getData('idPNMain',i)
    print('searching for idPNMain = ',idPNMain)
#    print('int(idPNMain) = ',int(idPNMain))
    idx = np.where(hashPNMain_idPNMains == idPNMain)[0]
    print('idx = ',idx)
    PNG = hashPNMain.getData('PNG',idx)
    RA = hashPNMain.getData('RAJ2000',idx)
    DEC = hashPNMain.getData('DECJ2000',idx)
    classification = ''
    note = comments.getData('comment',i)
    if idPNMain not in wrData.getData('idPNMain'):
        wrData.append([idPNMain,PNG,RA,DEC,classification,note])
    else:
        idx = np.where(np.array(wrData.getData('idPNMain')) == idPNMain)[0]
        print('idPNMain already in wrData at position ',idx)
        wrData.setData('notes',idx,wrData.getData('notes',idx)+' '+note)
print('wrData.size() = ',wrData.size())

hash_PNGs = np.array(hashPNMain.getData('PNG'))
print('hash_PNGs = ',hash_PNGs)
for i in range(listA.size()):
    png = listA.getData('PNG',i)
    png = png[png.find('G')+1:]
    print('listA: i=',i,': searching for PNG ',png)
    idx = np.where(hash_PNGs == png)[0]
    print('idx = ',idx)
    idPNMain = hashPNMain.getData('idPNMain',idx)
    RA = listA.getData('RA',i).replace(' ',':')
    DEC = listA.getData('DEC',i).replace(' ',':')
    classification = listA.getData('classification',i)
    note = ''
    if idPNMain not in wrData.getData('idPNMain'):
        wrData.append([idPNMain,png,RA,DEC,classification,note])
    else:
        idx = np.where(np.array(wrData.getData('idPNMain')) == idPNMain)[0]
        print('found at position ',idx)
        wrData.setData('classification',idx,classification)
print('wrData.size() = ',wrData.size())

hash_names = np.array([name.lower().replace(' ','') for name in hashCNames.getData('Name')])
print('hash_names = ',hash_names)
print('listB')
print('listB.size() = ',listB.size())
for i in range(listB.size()):
    name = listB.getData('Name',i).lower()
    print('searching for name ',name)
    idx = np.where(hash_names == name)[0]
    print('idx = ',idx)
    if len(idx) == 0:
        print('could not find name <'+name+'> in hashCNames')
        idx = np.where(hash_names == name[:3]+'J'+name[3:])[0]
        if len(idx) == 0:
            idx = np.where(hash_names == 'pn'+name)[0]
            if len(idx) == 0:
                print('RA = ',listB.getData('RA',i),', DEC = ',listB.getData('DEC',i))
                STOP
    idPNMain = hashCNames.getData('idPNMain',idx)
    idx = np.where(hashPNMain_idPNMains == idPNMain)[0]
    png = hashPNMain.getData('PNG',idx)
    RA = listB.getData('RA',i)
    DEC = listB.getData('DEC',i)
    classification = ''
    note = ''
    if idPNMain not in wrData.getData('idPNMain'):
        wrData.append([idPNMain,png,RA,DEC,classification,note])
    else:
        idx = np.where(np.array(wrData.getData('idPNMain')) == idPNMain)[0]
        print('found at position ',idx)
        wrData.setData('classification',idx,classification)
print('wrData.size() = ',wrData.size())

comments_idPNMains = comments.getData('idPNMain')
nCoords = 0
for i in range(wrData.size()):
    idPNMain = wrData.getData('idPNMain',i)
    for j in range(CSCoords.size()):
        if CSCoords.getData('idPNMain',j) == idPNMain:
            wrData.setData('RA',i,CSCoords.getData('CS_RAJ2000',j))
            wrData.setData('DEC',i,CSCoords.getData('CS_DECJ2000',j))
            nCoords += 1
    if idPNMain not in comments_idPNMains:
        for j in range(userComments.size()):
            if userComments.getData('idPNMain',j) == idPNMain:
                if userComments.getData('comment',j) not in wrData.getData('notes',i):
                    wrData.setData('notes',i,wrData.getData('notes',i)+'; '+userComments.getData('comment',j))

csvFree.writeCSVFile(wrData,'/Users/azuri/daten/uni/HKU/interns_projects/oriol/all_WR_CSPN.csv')

hash_wr = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/all_WR_CSPN.csv')
print('nCoords = ',nCoords)
with open('/Users/azuri/daten/uni/HKU/interns_projects/oriol/simbad_query.txt','w') as f:
    for i in range(hash_wr.size()):
        f.write('query coo %s %s radius=2s frame=FK5 epoch=2000\n' % (hash_wr.getData('RA',i), hash_wr.getData('DEC',i)))

wrDataOld = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/all_WR_CSPN.csv.bak')
idPNMainsOld = wrDataOld.getData('idPNMain')
hash_wr_temp = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/interns_projects/oriol/all_WR_CSPN.csv')
for i in range(wrData.size()):
    idPNMain = wrData.getData('idPNMain',i)
    if idPNMain not in idPNMainsOld:
        print('idPNMain ',idPNMain,' is new to list')
        hash_wr_temp.removeRow(np.where(np.array(hash_wr_temp.getData('idPNMain')) == idPNMain)[0])
idPNMainsTemp = hash_wr_temp.getData('idPNMain')
for i in range(wrData.size()):
    idPNMain = wrData.getData('idPNMain',i)
    if idPNMain not in idPNMainsTemp:
        hash_wr_temp.append(hash_wr.getData(i))

idPNMainsTemp = np.array(hash_wr_temp.getData('idPNMain'))
for idPNMain in ['14820','2915','1185','1266','1201','1250','2379','2503','190','125','184','207','2907','3060','2093','2503','1002','1185','1321','4246','205']:
    idx = np.where(idPNMainsTemp == idPNMain)[0]
    if len(idx) > 0:
        hash_wr_temp.removeRow(idx)
        print('removed idPNMain ',idPNMain)
        idPNMainsTemp = np.array(hash_wr_temp.getData('idPNMain'))
csvFree.writeCSVFile(hash_wr_temp,'/Users/azuri/daten/uni/HKU/interns_projects/oriol/all_WR_CSPN.csv')

idPNMainsTemp = np.array(hash_wr_temp.getData('idPNMain'))
for idPNMain in idPNMainsTemp:
    if len(np.where(idPNMainsTemp == idPNMain)[0]) > 1:
        print('found idPNMain ',idPNMain,' more than once')

setNames = np.array(spectraInfo.getData('name'))
fitsFiles_idPNMains = np.array(fitsFiles.getData('idPNMain'))
with open('/Users/azuri/daten/uni/HKU/interns_projects/oriol/copySpectra.sh','w') as f:
    for i in range(hash_wr_temp.size()):
        idPNMain = hash_wr_temp.getData('idPNMain',i)
        idx = np.where(fitsFiles_idPNMains == idPNMain)[0]
        print('idPNMain = ',idPNMain,': idx = ',idx)
        for id in idx:
            setName = fitsFiles.getData('setname',id)
            if 'HLA' not in setName:
                fName = fitsFiles.getData('fileName',id)
                path = spectraInfo.getData('path',np.where(setNames == setName))[0]
                f.write('cp -p /data/mashtun/'+path+'/'+fName+' wr_spectra/\n')
        if len(idx) == 0:
            print('did find no spectra for idPNMain ',idPNMain)
