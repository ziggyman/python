import CSVFree,CSVData
fName  = '/Users/azuri/daten/uni/HKU/interns_projects/alex/hash_fitsfiles_25042023.csv'
csvFitsFiles = CSVFree.readCSVFile(fName)

allDirs = []
for i in range(csvFitsFiles.size()):
    setname = csvFitsFiles.getData('setname',i)
    if setname not in allDirs:
        allDirs.append(setname)

print('allDirs = ',allDirs)
