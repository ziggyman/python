import csvData
import csvFree
import numpy as np

csvA = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/ziggy_Calern_PN_candidates_May2019.csv')
csvA.append(csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/ziggy_Calern_PN_candidates_May2019_II.csv'))
csvA.append(csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/ziggy_Calern_PN_candidates_May2019_III.csv'))

sorted = np.sort(csvFree.convertStringVectorToUnsignedVector(csvA.getData('idPNMain')))

print('sorted = ',sorted)

fNameOut = '/Users/azuri/daten/uni/HKU/HASH/ziggy_Calern_PN_candidates_May2019_out.csv'

csvOut = csvData.CSVData()
csvOut.header = csvA.header

for i in np.arange(0,len(sorted),1):
    found = csvA.find('idPNMain',str(sorted[i]))
    print('i = ',i,': found = ',found)
    csvOut.append(csvA.getData(found))

print('csvOut.size() = ',csvOut.size())

csvFree.writeCSVFile(csvOut,fNameOut)
