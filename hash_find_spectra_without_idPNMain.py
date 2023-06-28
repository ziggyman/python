import csvFree,csvData

spectraInfo = csvFree.readCSVFile('/Users/azuri//daten/uni/HKU/HASH/hash_spectraInfo_120623.csv')
fitsFiles = csvFree.readCSVFile('/Users/azuri//daten/uni/HKU/HASH/hash_fitsFiles_120623.csv')

found = fitsFiles.find('idPNMain','NULL')
print(len(found),': ',found)

fNames = []
setNames = []
with open('/Users/azuri//daten/uni/HKU/HASH/copy_spectra_without_idPNMain.sh','w') as f:
    f.write('mkdir spectra_without_idPNMain\n')
    for i in range(len(found)):
        path = spectraInfo.getData('path',spectraInfo.find('name',fitsFiles.getData('setname',found[i]))[0])
        print('setname = ',fitsFiles.getData('setname',found[i]),', fileName = ',fitsFiles.getData('fileName',found[i]))
        if not fitsFiles.getData('setname',found[i]) in setNames:
            f.write('mkdir spectra_without_idPNMain/'+fitsFiles.getData('setname',found[i])+'\n')
            setNames.append(fitsFiles.getData('setname',found[i]))
        f.write('cp /data/mashtun/'+path+fitsFiles.getData('fileName',found[i])+' spectra_without_idPNMain/'+fitsFiles.getData('setname',found[i])+'/\n')
