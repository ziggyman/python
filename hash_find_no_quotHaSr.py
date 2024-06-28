import numpy as np
import csvFree,csvData

hash_shs = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_SHS.csv')
hash_quot_HaSr = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/hash_quot_HaSr.csv')

idPNMains_quot_HaSr = np.array(hash_quot_HaSr.getData('idPNMain'))
idPNMain_wha = []
idPNMain_wsr = []
idPNMain_has_quot_HaSr = hash_quot_HaSr.getData('idPNMain')
idPNMain_no_quot_HaSr = []
for i in range(hash_shs.size()):
    idPNMain = hash_shs.getData('idPNMain',i)
    if hash_shs.getData('filename',i).endswith('wsr.fits'):
        idPNMain_wsr.append(idPNMain)
    elif hash_shs.getData('filename',i).endswith('wha.fits'):
        idPNMain_wha.append(idPNMain)

for id in idPNMain_wha:
    if id in idPNMain_wsr:
        if id not in idPNMain_has_quot_HaSr:
            idPNMain_no_quot_HaSr.append(id)
print('idPNMains_no_quot_HaSr = ',len(idPNMain_no_quot_HaSr),': ',idPNMain_no_quot_HaSr)
