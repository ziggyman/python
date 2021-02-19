import csv
import os

path = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/spectra'
fitsFileName = os.path.join(path,'hash_fitsfiles_200121.csv')
spectraInfoFileName = os.path.join(path,'hash_spectraInfo_200121.csv')
listFileName = os.path.join(path,'HASH_spectres.csv')

outFileName = os.path.join(path,'copySpectra.sh')
outFileList = os.path.join(path,'spectraFiles.txt')

nSpectra = 0
with open(outFileName,'w') as f:
    with open(outFileList,'w') as fl:
        with open(listFileName,'r') as listFile:
            spectraList = csv.DictReader(listFile)
            for item in spectraList:
                id = item['Id HASH']
                print('searching for id ',id,', name = ',item['Name'])
                with open(fitsFileName,'r') as fitsFile:
                    fitsList = csv.DictReader(fitsFile)
                    for fits in fitsList:
                        if fits['idPNMain'] == id:
                            with open(spectraInfoFileName,'r') as spectraInfoFile:
                                spectraInfoList = csv.DictReader(spectraInfoFile)
                                setName = fits['setname']
                                fitsName = fits['fileName']
                                print('id ',id,' found in setName = ',setName)
                                for specInfo in spectraInfoList:
                                    if specInfo['name'] == setName:
                                        fitsPath = os.path.join('/data/mashtun/',specInfo['path'])
                                        f.write('cp '+os.path.join(fitsPath,fitsName)+' FRA_spectra/\n')
                                        fl.write(item['Name']+': '+fitsName+'\n')
                                        nSpectra += 1
print('found ',nSpectra,' files')
