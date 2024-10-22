import csvFree,csvData

table = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/publications/clusters/hash_CSCoords_081024.csv')
PNMain = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/publications/clusters/hash_PNMain_081024.csv')

idPNMains = []
for i in range(table.size()):
    idPNMain = table.getData('idPNMain',i)
    for j in range(PNMain.size()):
        if idPNMain == PNMain.getData('idPNMain',j):
            if PNMain.getData('PNstat',j) in ['T','L','P']:
                if idPNMain not in idPNMains:
                    idPNMains.append(idPNMain)
print('found ',len(idPNMains),' CSPN')
