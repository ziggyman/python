import csv

hashFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_PNSpectra_Sources_FitsFiles_290920.csv'
idsFileNameOut = '/Users/azuri/daten/uni/HKU/HASH/hash_saao_spectra_ids.dat'

ids = []
with open(hashFileName) as csvFile:
    hashData = csv.DictReader(csvFile)

    with open(idsFileNameOut,'w') as f:
        for row in hashData:
            setname = row['setname']
            if (setname[0:4] == 'SAAO') or (setname == 'FRA'):
                id = row['idPNMain']
                if id != 'NULL':
                    foundId = False
                    foundSet = False
                    for i in range(len(ids)):
                        if id == ids[i][0]:
                            foundId = True
                            if setname == ids[i][1]:
                                foundSet = True
                            if foundId and not foundSet:
                                f.write('%s,' % (id))
                    ids.append([id,setname])
