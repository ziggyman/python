import csvFree,csvData

tablesInfo = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_MainGPN.tablesInfo.csv')

tbIRFlux = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_MainGPN.tbIRFlux.csv')
tbRadioCont = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_MainGPN.tbRadioCont.csv')
varTypes = []
IRvarNames = []
RCvarNames = []
for i in range(tablesInfo.size()):
    if tablesInfo.getData('varType',i) not in varTypes:
        varTypes.append(tablesInfo.getData('varType',i))
    if tablesInfo.getData('varTable',i) == 'tbIRFlux':
        data = tbIRFlux.getData(tablesInfo.getData('varColumn',i))
        IRvarNames.append(tablesInfo.getData('varColumn',i))
        print(tablesInfo.getData('varName',i),': data = ',data)
        if tablesInfo.getData('varType',i) == 'REAL':
            for j in range(len(data)):
                if data[j] != 'NULL':
                    a = float(data[j])
                    print(tablesInfo.getData('varColumn',i),': data[',j,'] = <',data[j],'>: a = ',a)

    if tablesInfo.getData('varTable',i) == 'tbRadioCont':
        data = tbRadioCont.getData(tablesInfo.getData('varColumn',i))
        RCvarNames.append(tablesInfo.getData('varColumn',i))
        print(tablesInfo.getData('varName',i),': data = ',data)
        if tablesInfo.getData('varType',i) == 'REAL':
            for j in range(len(data)):
                if data[j] != 'NULL':
                    a = float(data[j])
                    print(tablesInfo.getData('varColumn',i),': data[',j,'] = <',data[j],'>: a = ',a)


print('varTypes = ',varTypes)
print('IRvarNames = ',IRvarNames)
print('RCvarNames = ',RCvarNames)

fullView = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_MainGPN.GPNFullView.csv')
fullViewCols = []
for i in range(len(fullView.header)):
    if ('_Flux_' in fullView.header[i]) or ('_Mag_' in fullView.header[i]):
        fullViewCols.append(fullView.header[i])
        data = fullView.getData(fullView.header[i])
        for j in range(len(data)):
            if data[j] != 'NULL':
                a = float(data[j])
print('checked columns ',fullViewCols)
