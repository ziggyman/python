import csvFree,csvData
import numpy as np
import matplotlib.pyplot as plt
from myUtils import angularDistancePyAsl,lonLatToRaDec,hmsToDeg,dmsToDeg,degToDMS,degToHMS,readmeToCSV

dataFileName = '/Users/azuri/daten/uni/HKU/HASH/allGalacticTLP_090822.csv'

data = csvFree.readCSVFile(dataFileName)

print('header = ',data.header)

noCSPN = data.find('CS_Glon','')
print(noCSPN)
print('data.size() = ',data.size())
for i in np.arange(len(noCSPN),0,-1):
    data.removeRow(i)
print('data.size() = ',data.size())

nT = len(data.find('PNstat','T'))
nL = len(data.find('PNstat','L'))
nP = len(data.find('PNstat','P'))

print('nT = ',nT)
print('nL = ',nL)
print('nP = ',nP)

nE = len(data.find('mainClass','E'))
nR = len(data.find('mainClass','R'))
nB = len(data.find('mainClass','B'))
nI = len(data.find('mainClass','I'))
nA = len(data.find('mainClass','A'))
nS = len(data.find('mainClass','S'))

print('nE+nR+nB+nI+nA+nS = ',nE+nR+nB+nI+nA+nS)
noMainClass = data.size() - (nE+nR+nB+nI+nA+nS)

def plot():
    labels = '%d T' % (nT), '%d L' % (nL), '%d P' % (nP)
    sizes = [nT,nL,nP]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=False, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('/Users/azuri/daten/uni/HKU/HASH/allGalacticTLP_090822_TLP.pdf', bbox_inches='tight')
    plt.show()


    labels = '%d E' % (nE), '%d R' % (nR), '%d B' % (nB), '%d A' % (nA), '%d S' % (nS), '%d I' % (nI), '%d not classified' % (noMainClass)
    sizes = [nE,nR,nB,nA,nS,nI,noMainClass]
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=False, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('/Users/azuri/daten/uni/HKU/HASH/allGalacticTLP_090822_ERBIAS.pdf', bbox_inches='tight')
    plt.show()

def getIDPNMainFromPNG(pngs, hashMainGPN):
    ids = []
    for png in pngs:
        idPos = hashMainGPN.find('PNG',png)
        if idPos[0] < 0:
            print('ERROR: could not find PNG '+png)
            STOP
        ids.append(hashMainGPN.getData('idPNMain',idPos[0]))
    return ids


if False:
    gonzalesFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51/tablea1.csv'
    hashMainGPNFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_MainGPN_200622.csv'
    hashMainGPN = csvFree.readCSVFile(hashMainGPNFileName)

    gonzales = csvFree.readCSVFile(gonzalesFileName)
    pngs = gonzales.getData('PNG')
    idPNMainsGonzales = getIDPNMainFromPNG(pngs, hashMainGPN)
    gonzales.addColumn('idPNMain',idPNMainsGonzales)

    data = csvFree.readCSVFile(dataFileName)
    CS_Glons = []
    CS_Glats = []
    hashIDs = []
    mainClasses = []
    notFound = []
    for i in range(len(idPNMainsGonzales)):
        idx = data.find('idPNMain',idPNMainsGonzales[i])[0]
        if idx >= 0:
            CS_Glons.append(data.getData('CS_Glon', idx))
            CS_Glats.append(data.getData('CS_Glat', idx))
            hashIDs.append(data.getData('idPNMain', idx))
        else:
            print('ERROR: idPNMain ',idPNMainsGonzales[i],' not found')
            notFound.append(idPNMainsGonzales[i])
    print('notFound = ',len(notFound),': ',notFound)
    #STOP

    print('len(CS_Glons) = ',len(CS_Glons))
    print('len(CS_Glats) = ',len(CS_Glats))

    hashDRADECs = []
    for i in range(len(CS_Glons)):
        if CS_Glons[i] != '':
            hashDRADECs.append(lonLatToRaDec(float(CS_Glons[i]),float(CS_Glats[i])))
            hashDRADECs[len(hashDRADECs)-1].append(hashIDs[i])

    print('hashDRADECs = ',len(hashDRADECs),': ',hashDRADECs)

    distances = []
    for i in range(len(hashDRADECs)):
        idx = gonzales.find('idPNMain',hashDRADECs[i][2])[0]
        distances.append(angularDistancePyAsl(hashDRADECs[i][0],
                                            hashDRADECs[i][1],
                                            float(gonzales.getData('RAdeg',idx)),
                                            float(gonzales.getData('DEdeg',idx))))
        distances[len(distances)-1] = [distances[len(distances)-1]*3600.,
                                    hashDRADECs[i][2],
                                    degToHMS(float(hashDRADECs[i][0])),
                                    degToDMS(float(hashDRADECs[i][1])),
                                    degToHMS(float(gonzales.getData('RAdeg',idx))),
                                    degToDMS(float(gonzales.getData('DEdeg',idx))),
                                    ]
    print('distances = ',distances)
    possibleWrongCSinGonzales = []
    for i in range(len(distances)):
        if distances[i][0] > 2.:
            print(distances[i])
            possibleWrongCSinGonzales.append(distances[i])
    print('possibleWrongCSinGonzales = ',len(possibleWrongCSinGonzales),': ',possibleWrongCSinGonzales)
    xmax = 381
    plt.hist([i[0] for i in distances],bins=int(xmax/2.5))
    plt.xlim([0,xmax])
    plt.ylim([0,50])
    plt.xlabel('Distance HASH to GS2021 ["]')
    plt.ylabel('Number of CSPN')
    plt.savefig('/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51//allGalacticTLP_distance_histogram_y0-50_x0-380.pdf', bbox_inches='tight')
    plt.show()

    xmax=50
    plt.hist([i[0] for i in distances],bins=int(381/2.5))
    plt.xlim([0,xmax])
    plt.ylim([0,50])
    plt.xlabel('Distance HASH to GS2021 ["]')
    plt.ylabel('Number of CSPN')
    plt.savefig('/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51//allGalacticTLP_distance_histogram_y0-50_x0-50.pdf', bbox_inches='tight')
    plt.show()

chornayReadMeFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A110_Chorney2021/J_A+A_656_A110/ReadMe'
chornayDataFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A110_Chorney2021/J_A+A_656_A110/tablea1.dat'
chornayCSVFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A110_Chorney2021/J_A+A_656_A110/tablea1.csv'

readmeToCSV(chornayReadMeFileName,chornayDataFileName,chornayCSVFileName)
