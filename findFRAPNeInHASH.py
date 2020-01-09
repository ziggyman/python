import csvData
import csvFree

longList = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/Jan2020/TableI_507PNG_08122019_JaunePasDansHASH.csv'
shortList = '/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/Jan2020/TableI_105objets_20112019_JaunePasDansHASH.csv'

keepColumns= ['NOM','AD:(J2000)','DEC (J2000)']
print('keepColumns = ',keepColumns)
for iList in [longList, shortList]:
    csv = csvFree.readCSVFile(iList,',',False)
    if csv.header[0] == '\ufeffNOM':
        csv.renameColumn('\ufeffNOM','NOM')
        print('renamed column "\ufeffNOM" to <'+csv.header[0]+'>')
    if csv.header[len(csv.header)-1] == 'HASH ID\r':
        csv.renameColumn('HASH ID\r','HASH ID')
        print('renamed column "HASH ID\r" to <'+csv.header[len(csv.header)-1]+'>')
    print('csv = ',csv)
    print('csv.header = ',csv.header)
    for keyword in csv.header:
        print('checking keyword <'+keyword+'>')
        if keyword not in keepColumns:
            print('keyword <'+keyword+'> not found in keepColumns => removing column')
            csv.removeColumn(keyword)
    print('csv.header = ',csv.header)
    print('csv.data = ',csv.data)

