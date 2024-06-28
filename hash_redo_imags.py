import numpy as np
import csvFree,csvData

pnMain = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_210624.csv')
idPNMains = pnMain.getData('idPNMain')
with open('/Users/azuri/daten/uni/HKU/HASH/redo_images.txt','w') as f:
    ids = []
    new = True
    for i in np.arange(len(idPNMains) - 1, -1, -1):
        ids.append(idPNMains[i])
        if new:
            f.write('hashpn fetch all %d' % (int(idPNMains[i])))
            new = False
        else:
            f.write(',%d' % (int(idPNMains[i])))
        if i % 100 == 0:
            f.write(' -w force\n')
            f.write('hashpn brew all %d' % (int(ids[0])))
            for id in ids[1:]:
                f.write(',%d' % (int(id)))
            f.write(' -w\n')
            ids = []
            new = True
