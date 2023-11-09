import csvFree,csvData

all = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/CFP-GE_Oct09.txt','\t',False)
print('all = ',all.header)
print(all.data)
names_all = all.getData('Name')
names_all = [name[:name.find(',')].lower() for name in names_all]
print('names_all = ',names_all)

new = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/possible_authors.csv')
print('new = ',new.header)
print(new.data)
names_new = new.getData('name')
names_new = [name[name.rfind(' ')+1:].lower() for name in names_new]
print('names_new = ',names_new)

iaus384 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/galaxies special issue guest editor/CFP-GE_Oct09.txt','\t',False)
names_iaus384 = iaus384.getData('Name')
names_iaus384 = [name[name.rfind(' ')+1:].lower() for name in names_iaus384]
print('names_iaus384 = ',names_iaus384)
