import csvFree,csvData

cat = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/hash_LiParker1.csv')
catIDs = cat.getData('idPNMain')

with open('/Users/azuri/daten/uni/HKU/observing/hash_id_Yushan.txt','r') as f:
    ids = f.readlines()

with open('/Users/azuri/daten/uni/HKU/observing/hash_id_Yushan.sql','w') as f:
    iStart = 32
    for id in ids:
        id = id.strip()
        if id not in catIDs:
            f.write("INSERT INTO `Li-Parker-1`(`idLi-Parker-1`,`idPNMain`,`mapFlag`) VALUES (%d,%d,'y');\n" % (iStart, int(id)))
            iStart += 1
