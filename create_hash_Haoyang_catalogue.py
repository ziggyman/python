import csvFree,csvData

ids = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/hash_Haoyang_idPNMains.csv').getData('idPNMain')

i = 1
with open('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/hash_Haoyang.sql','w') as f:
    for id in ids:
        f.write("INSERT INTO `Haoyang`(`idHaoyang`,`idPNMain`,`mapFlag`) VALUES (%d,%d,'y');\n" % (i,int(id)))
        i += 1
