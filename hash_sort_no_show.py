import csvData,csvFree

inFileName = '/Users/azuri/daten/uni/HKU/HASH/hash-no-show.csv'

inData = csvFree.readCSVFile(inFileName)
print('inData.size() = ',inData.size())

pnStat = inData.getData('PNstat')
print('pnStat = ',pnStat)

tlps = csvData.CSVData()
tlps.header = inData.header

newCandidates = csvData.CSVData()
newCandidates.header = inData.header

others = csvData.CSVData()
others.header = inData.header

for i in range(inData.size()):
    pnStat = inData.getData('PNstat',i)
    if pnStat in ['T','L','P']:
        tlps.append(inData.getData(i))
    elif pnStat == 'c':
        newCandidates.append(inData.getData(i))
    else:
        others.append(inData.getData(i))

csvFree.writeCSVFile(tlps,inFileName[:-4]+'_TLPs.csv')
csvFree.writeCSVFile(newCandidates,inFileName[:-4]+'_newCandidates.csv')
csvFree.writeCSVFile(others,inFileName[:-4]+'_others.csv')
