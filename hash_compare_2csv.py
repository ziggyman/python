import csvData
import csvFree

a = '/Users/azuri/daten/uni/HKU/HASH/barlow1_sql.csv'
b = '/Users/azuri/daten/uni/HKU/HASH/barlow1_hash.csv'

csvA = csvFree.readCSVFile(a)
csvB = csvFree.readCSVFile(b)

idsA = csvA.getData('idPNMain')
idsB = csvB.getData('idPNMain')

print('idsA = ',idsA)
print('idsB = ',idsB)

for idA in idsA:
    if idA in idsB:
        print('found id '+idA)
    else:
        print('could not find id '+idA)

for i in range(csvA.size()):
    if csvA.getData('idPNMain',i) == '4506':
        print(csvA.getData(i))