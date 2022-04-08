import csvFree,csvData

fitsFilesFile = '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/vanessa/hash_FitsFiles.csv'
spectraInfoFile = '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/vanessa/hash_spectraInfo.csv'

csvFiles = csvFree.readCSVFile(fitsFilesFile)
csvInfo = csvFree.readCSVFile(spectraInfoFile)

#idFiles = ['/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/vanessa/getSpectraID_bipolar.dat',
#           '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/vanessa/getSpectraID_elliptical.dat',
#           '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/vanessa/getSpectraID_round.dat',
#           ]
#idFiles = ['/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/burton/bipolar.txt',
#           '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/burton/elliptical.txt',
#           '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/burton/round.txt',
#           ]
idFiles = ['/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/helen/round.txt',
#           '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/helen/bipolar.txt',
#           '/Users/azuri/daten/uni/HKU/HASH/HASH_spectra/helen/elliptical.txt',
           ]

for idFile in idFiles:
    prefix = idFile[idFile.rfind('/')+1:idFile.rfind('.')]
    print('prefix = <'+prefix+'>')
    with open(idFile,'r') as f:
        ids = f.readlines()
    for i in range(len(ids)):
        ids[i] = ids[i].strip()
    print('ids = ',ids)
    with open(idFile[:idFile.rfind('/')+1]+prefix+'_fetch.bat','w') as w:
        w.write('mkdir helen\n')
        w.write('mkdir helen/'+prefix+'\n')
        for id in ids:
            specNames = csvFiles.getData('fileName',csvFiles.find('idPNMain',id))
            print('specNames = ',specNames)
            w.write('mkdir helen/'+prefix+'/'+str(id)+'\n')
            for iSpec in range(len(specNames)):
                w.write('cp /data/mashtun/'+csvInfo.getData('path',csvInfo.find('name',csvFiles.getData('setname',csvFiles.find('fileName',specNames[iSpec])[0]))[0])+specNames[iSpec]+' helen/'+prefix+'/'+str(id)+'/\n')
