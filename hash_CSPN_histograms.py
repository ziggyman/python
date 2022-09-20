from pickle import STOP
import csvFree,csvData
import numpy as np
import matplotlib.pyplot as plt
from myUtils import angularDistancePyAsl,lonLatToRaDec,hmsToDeg,dmsToDeg,degToDMS,degToHMS,readmeToCSV

dataFileName = '/Users/azuri/daten/uni/HKU/HASH/allGalacticTLP_090822.csv'

data = csvFree.readCSVFile(dataFileName)

print('header = ',data.header)


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

labels = '%d T' % (nT), '%d L' % (nL), '%d P' % (nP)
sizes = [nT,nL,nP]
fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.savefig('/Users/azuri/daten/uni/HKU/HASH/allGalacticTLP_withHASHCSPN_090822_TLP.pdf', bbox_inches='tight')
plt.show()


labels = '%d E' % (nE), '%d R' % (nR), '%d B' % (nB), '%d A' % (nA), '%d S' % (nS), '%d I' % (nI), '%d not classified' % (noMainClass)
sizes = [nE,nR,nB,nA,nS,nI,noMainClass]
fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.savefig('/Users/azuri/daten/uni/HKU/HASH/allGalacticTLP_withHASHCSPN_090822_ERBIAS.pdf', bbox_inches='tight')
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

def makeAngDistHistograms(dataCSV, hashMainGPN, plotNameSuffix):

    dataAllGalacticPNeTLP = csvFree.readCSVFile(dataFileName)

    pngs = dataCSV.getData('PNG')
    idPNMainsdataCSV = getIDPNMainFromPNG(pngs, hashMainGPN)
    dataCSV.addColumn('idPNMain',idPNMainsdataCSV)

    dataAllGalacticPNeTLP = csvFree.readCSVFile(dataFileName)
    CS_Glons = []
    CS_Glats = []
    hashIDs = []
    mainClasses = []
    notFound = []
    nStellarPNe = 0
    for i in range(len(idPNMainsdataCSV)):
        idx = dataAllGalacticPNeTLP.find('idPNMain',idPNMainsdataCSV[i])[0]
        if idx >= 0:
            CS_Glons.append(dataAllGalacticPNeTLP.getData('CS_Glon', idx))
            CS_Glats.append(dataAllGalacticPNeTLP.getData('CS_Glat', idx))
            hashIDs.append(dataAllGalacticPNeTLP.getData('idPNMain', idx))
            if dataAllGalacticPNeTLP.getData('mainClass',idx) == 'S':
                nStellarPNe += 1
        else:
            print('ERROR: idPNMain ',idPNMainsdataCSV[i],' not found')
            notFound.append(idPNMainsdataCSV[i])
    print('notFound = ',len(notFound),': ',notFound)
    print('nStellarPNe = ',nStellarPNe)
    STOP

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
        idx = dataCSV.find('idPNMain',hashDRADECs[i][2])[0]
        distances.append(angularDistancePyAsl(hashDRADECs[i][0],
                                            hashDRADECs[i][1],
                                            float(dataCSV.getData('RAdeg',idx)),
                                            float(dataCSV.getData('DEdeg',idx))))
        distances[len(distances)-1] = [distances[len(distances)-1]*3600.,
                                    hashDRADECs[i][2],
                                    degToHMS(float(hashDRADECs[i][0])),
                                    degToDMS(float(hashDRADECs[i][1])),
                                    degToHMS(float(dataCSV.getData('RAdeg',idx))),
                                    degToDMS(float(dataCSV.getData('DEdeg',idx))),
                                    ]
    print('distances = ',distances)
    possibleWrongCSindataCSV = []
    for i in range(len(distances)):
        if distances[i][0] > 2.:
            print(distances[i])
            possibleWrongCSindataCSV.append(distances[i])
    print('possibleWrongCSindataCSV = ',len(possibleWrongCSindataCSV),': ',possibleWrongCSindataCSV)
    xmax = 178#381
    ymax = 50
    plt.hist([i[0] for i in distances],bins=int(xmax/2.5))
    plt.xlim([0,xmax])
    plt.ylim([0,ymax])
    plt.xlabel('Distance HASH to C+W 2021 ["]')
    plt.ylabel('Number of CSPN')
    plt.savefig('/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51//allGalacticTLP_distance_histogram_'+plotNameSuffix+'_y0-'+str(ymax)+'_x0-'+str(xmax)+'.pdf', bbox_inches='tight')
    plt.show()

    xmax=50
    plt.hist([i[0] for i in distances],bins=int(381/2.5))
    plt.xlim([0,xmax])
    plt.ylim([0,50])
    plt.xlabel('Distance HASH to C+W 2021 ["]')
    plt.ylabel('Number of CSPN')
    plt.savefig('/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51//allGalacticTLP_distance_histogram_'+plotNameSuffix+'_y0-'+str(ymax)+'_x0-'+str(xmax)+'.pdf', bbox_inches='tight')
    plt.show()

def getIdFromNames(csvHashNames, names):
    ids = []
    for name in names:
        idx = -1
        for i in range(csvHashNames.size()):
            if csvHashNames.getData('Name',i).strip().replace(' ','').lower() == name.strip().replace(' ','').lower():
                idx = i
        if idx < 0:
            print('could not find name <'+name+'>')
            ids.append('-1')
        else:
            ids.append(csvHashNames.getData('idPNMain',idx))
    return ids

def findDistanceToQuentin(dataCSV,csvQuentin, hashMainGPN, setName):
    pngs = dataCSV.getData('PNG')
    idPNMainsdataCSV = getIDPNMainFromPNG(pngs, hashMainGPN)
    dataCSV.addColumn('idPNMain',idPNMainsdataCSV)

    quentinIds = getIdFromNames(csvHashNames, csvQuentin.getData('ID1'))
    print('quentinIds = ',quentinIds)

    nNotFoundInGaiaSet = 0
    dists = []
    for i in range(len(quentinIds)):
        if int(quentinIds[i]) >= 0:
            idx = dataCSV.find('idPNMain',str(quentinIds[i]))[0]
            if idx < 0:
                print('did not find Quentins idPNMain ',quentinIds[i])
                nNotFoundInGaiaSet += 1
            else:
                print('found QuentinIds[',i,'] = ',quentinIds[i],' in ',setName,' at position ',idx)
                print('csvQuentin.getData(',i,') = ',csvQuentin.getData(i))
                qdRA = hmsToDeg(csvQuentin.getData('RA',i).replace('\r',''))
                print('qdRA = ',qdRA)
                qdDEC = dmsToDeg(csvQuentin.getData('DEC\r',i).replace('\r',''))
                print('qdDEC = ',qdDEC)
                gaiadRA = float(dataCSV.getData('RAdeg',idx))
                print('gaiadRA = ',gaiadRA)
                gaiadDEC = float(dataCSV.getData('DEdeg',idx))
                print('gaiadDEC = ',gaiadDEC)
                dists.append([quentinIds[i],
                              angularDistancePyAsl(qdRA,
                                                  qdDEC,
                                                  gaiadRA,
                                                  gaiadDEC) * 3600.,
                              qdRA,
                              qdDEC,
                              gaiadRA,
                              gaiadDEC])
                print('dists[',len(dists)-1,'] = ',dists[len(dists)-1])
    print('did not find ',nNotFoundInGaiaSet,' in set ',setName)
    #Stop
    nGTOne = 0
    nGTTwo = 0
    for dist in dists:
        print(dist)
        if dist[1] > 1.:
            nGTOne += 1
        if dist[1] > 2.:
            nGTTwo += 1
    print('nGTOne = ',nGTOne)
    print('nGTTwo = ',nGTTwo)
    #Stop

    dist = [d[1] for d in dists]
    xmax=int(np.max(dist))
    plt.hist(dist, bins=int(np.max(dist)))
    plt.xlim([0,xmax])
#    plt.ylim([0,ymax])
    plt.xlabel('Distance HASH to '+setName+' ["]')
    plt.ylabel('Number of CSPN')
    plt.savefig('/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51//allGalacticTLP_distance_histogram_'+setName+'_x0-'+str(xmax)+'.pdf', bbox_inches='tight')
    plt.show()

    xmax=int(np.max(dist))
    ymax=50
    plt.hist(dist, bins=int(np.max(dist)))
    plt.xlim([0,xmax])
    plt.ylim([0,ymax])
    plt.xlabel('Distance HASH to '+setName+' ["]')
    plt.ylabel('Number of CSPN')
    plt.savefig('/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51//allGalacticTLP_distance_histogram_'+setName+'_y0-'+str(ymax)+'_x0-'+str(xmax)+'.pdf', bbox_inches='tight')
    plt.show()

    return dists

if __name__ == '__main__':
    chornayReadMeFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A110_Chorney2021/J_A+A_656_A110/ReadMe'
    chornayDataFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A110_Chorney2021/J_A+A_656_A110/tablea1.dat'
    chornayCSVFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A110_Chorney2021/J_A+A_656_A110/tablea1.csv'

    #readmeToCSV(chornayReadMeFileName,chornayDataFileName,chornayCSVFileName)
    chornay = csvFree.readCSVFile(chornayCSVFileName)

    gonzalesFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51/tablea1.csv'
    hashMainGPNFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_MainGPN_200622.csv'
    hashMainGPN = csvFree.readCSVFile(hashMainGPNFileName)

    gonzales = csvFree.readCSVFile(gonzalesFileName)
#    makeAngDistHistograms(gonzales, hashMainGPN,'gonzales')
#    makeAngDistHistograms(chornay, hashMainGPN,'chornay')

    hashNamesFileName = '/Users/azuri/daten/uni/HKU/HASH/hash_tbCNames_200622.csv'
    quentinsCSPNFileName = '/Users/azuri/daten/uni/HKU/conferences/Tubingen08-2022/CSPN-added-130822.csv'
    csvHashNames = csvFree.readCSVFile(hashNamesFileName)
    csvQuentin = csvFree.readCSVFile(quentinsCSPNFileName)
    print('csvQuentin.header = ',csvQuentin.header)
    findDistanceToQuentin(gonzales,csvQuentin,hashMainGPN,'G-S')
    findDistanceToQuentin(chornay,csvQuentin,hashMainGPN,'CW')
