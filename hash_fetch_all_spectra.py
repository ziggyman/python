import csvFree,csvData
fName  = '/Users/azuri/daten/uni/HKU/interns_projects/alex/hash_fitsfiles_25042023.csv'
fNameInfo = '/Users/azuri/daten/uni/HKU/interns_projects/alex/hash_spectrainfo_25042023.csv'
csvFitsFiles = csvFree.readCSVFile(fName)
csvSpectraInfo = csvFree.readCSVFile(fNameInfo)

allDirs = []
for i in range(csvFitsFiles.size()):
    setname = csvFitsFiles.getData('setname',i)
    if setname not in allDirs:
        allDirs.append(setname)

print('allDirs = ',allDirs)

with open('/Users/azuri/daten/uni/HKU/interns_projects/alex/fetch_spectra.sh','w') as f:
    for d in allDirs:
        path = csvSpectraInfo.getData('path',csvSpectraInfo.find('name',d)[0])
        f.write('scp -pr tcooper@corona:/data/mashtun/'+path+' /Users/azuri/daten/uni/HKU/interns_projects/alex/allHASHspectra/\n')
