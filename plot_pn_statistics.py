import csvData
import csvFree
import numpy as np
import matplotlib.pyplot as plt

files = ['/Users/azuri/daten/uni/HKU/HASH/PN_statistics/northern PN.txt',
         '/Users/azuri/daten/uni/HKU/HASH/PN_statistics/southPN.txt']

# read one of our files and return the data for each column
def readFile(file):
    dataArr = []
    with open(file,'r') as f:
        lines = f.readlines()
        i = 0
        for line in lines:
            fields = line.split()
            dataArr.append(fields)
            print('line ',i,' = <',fields,'>')
            i += 1
    return dataArr

# read all our files and return one comma-separated-values object
def readFiles():
    csv = csvData.CSVData()
    for iFile in np.arange(0,len(files),1):
        dataTemp = readFile(files[iFile])
        if iFile == 0:
            csv.header = dataTemp[0]
            csv.data = dataTemp[1:]
        else:
            csv.append(dataTemp[1:])
    return csv

#csv = readFiles()
csvFile = '/Users/azuri/daten/uni/HKU/HASH/PN_statistics/PNdata.csv'
#csvFree.writeCSVFile(csv,csvFile)

pnData = csvFree.readCSVFile(csvFile)
print('pnData.header = ',pnData.header)
print('pnData.data = ',pnData.data)
print("pnData.getData('ID') = ",pnData.getData('ID'))

# read data from csv object and convert them to a numpy float array
NII = np.array([float(a) for a in pnData.getData('[NII]-6583')])
#NII = np.array([float(a) for a in pnData.getData('[NII]-6548')])
Halpha = np.array([float(a) for a in pnData.getData('Halpha')])
Hbeta = np.array([float(a) for a in pnData.getData('Hbeta')])
SII = np.array([float(a) for a in pnData.getData('[SII]-6716')])
OIII = np.array([float(a) for a in pnData.getData('[OIII]-5007')])
#OIII_4363 = np.array([float(a) for a in pnData.getData('[OIII]-4363')])

# set all values equal to zero to 0.001 to avoid inf
Halpha[np.where(Halpha == 0.)] = 0.001
Hbeta[np.where(Hbeta == 0.)] = 0.001

# ionization stages from [NII] and [OIII]
NII_abundance = NII / Halpha
OIII_abumdance = OIII / Hbeta
NII[np.where(NII == 0.)] = 0.001

OIII_over_NII = OIII_abumdance / NII_abundance

plt.hist(OIII_over_NII,20,range=[0.,5.])
plt.xlabel('[OIII] / [NII}')
plt.ylabel('counts')
plt.show()

# shock excitation
SII_over_Halpha = SII / Halpha

plt.hist(SII_over_Halpha,20,range=[0.,2.])
plt.xlabel('[SII] / H_alpha')
plt.ylabel('counts')
plt.show()

# temperature
#OIII_4363