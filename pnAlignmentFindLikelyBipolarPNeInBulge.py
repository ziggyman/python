import os
import csvFree,csvData

path = '/Users/azuri/daten/uni/HKU/PN alignment'
fNameIn = 'HASH_bipolar_likely_PN.csv'
csv = csvFree.readCSVFile(os.path.join(path,fNameIn))

csvOut = csvData.CSVData()
csvOut.header = csv.header
for i in range(csv.size()):
    if csv.getData('PNstat',i) == 'L':
        lon = float(csv.getData('Glon',i))
        lat = float(csv.getData('Glat',i))
        if ((lon < 10.) or (lon > 350.)) and (((lat > 0.) and (lat < 10.)) or ((lat < 0.) and (lat) > -10.)):
            print('lon = ',lon,', lat = ',lat)
            csvOut.append(csv.getData(i))
csvFree.writeCSVFile(csvOut,os.path.join(path,fNameIn[:-4]+'_bulge.csv'))
