import numpy as np

import csvFree,csvData

def DfromLogRAndTheta(logr, theta):
    return 206.265 * (10**logr) / theta

def rFromDAndTheta(D, theta):
    return D * theta / 206.265

inFileName3 = '/Users/azuri/daten/uni/HKU/PNe/Frew - The H-alpha surface brightness - radius relation- a robust statistical distance indicator for planetary nebulae_Supplementary_Data/Table_A3_new.csv'
inFileName4 = '/Users/azuri/daten/uni/HKU/PNe/Frew - The H-alpha surface brightness - radius relation- a robust statistical distance indicator for planetary nebulae_Supplementary_Data/Table_A4_new.csv'

r = None
for i in range(2):
    if i == 0:
        inFileName = inFileName4
    else:
        inFileName = inFileName3
    data = csvFree.readCSVFile(inFileName,'&',False)
#    print('header = ',len(data.header),': ',data.header)

    if i == 0:
        strvecD = data.getData('D_mean [kpc]')
        strvecD = [d[1:d.find(' ')] for d in strvecD]
#        print('strvecD = ',strvecD)
        distances = np.array(csvFree.convertStringVectorToDoubleVector(strvecD))
        strvecLogR = data.getData('log r[pc]')
    else:
        strvecD = data.getData('D [pc]')
        strvecD = [d.replace(' ','').strip('\t').rstrip('\t') for d in strvecD]
        strvecD = [d[1:d.find('$')] for d in strvecD]
#        print('strvecD = ',strvecD)
        distances = np.array(csvFree.convertStringVectorToDoubleVector(strvecD)) / 1000.
        strvecLogR = data.getData('log r [pc] ')

    strvecLogR = [r.replace('$','').replace(' ','').strip('\t').rstrip('\t') for r in strvecLogR]
#    print('strvecLogR = ',strvecLogR)
    logr = np.array(csvFree.convertStringVectorToDoubleVector(strvecLogR))

    strvecA = data.getData('a [arcsec]')
    strvecA = [a.replace('$','').replace(' ','').strip('\t').rstrip('\t') for a in strvecA]
#    print('strvecA = ',strvecA)

    strvecB = data.getData('b [arcsec]')
    strvecB = [a.replace('$','').replace(' ','').strip('\t').rstrip('\t') for a in strvecB]
#    print('strvecB = ',strvecB)

    theta = (np.array(csvFree.convertStringVectorToDoubleVector(strvecA)) + np.array(csvFree.convertStringVectorToDoubleVector(strvecB)))/4.

    Dab = DfromLogRAndTheta(logr, theta)

#    for j in range(data.size()):
#        print('distances[',j,'] = ',distances[j],', Dab[',j,'] = ',Dab[j])

    if i == 0:
        r = rFromDAndTheta(Dab,theta)
    else:
#        print('r = ',r)
        r = np.concatenate([r,rFromDAndTheta(Dab,theta)])

    names = data.getData('Name')
    names = [n.replace('~',' ').strip(' ').rstrip(' ').strip('\t').rstrip('\t') for n in names]
    for iName in range(len(names)):
        if names[iName] == 'TK 1':
            print('i = ',i,': name = TK 1')
            print(data.header)
            print(data.getData(iName))
    if i == 0:
        allNames = names
    else:
        for name in names:
            allNames.append(name)

    if i == 1:
        print('r = ',r.shape)
        sortedIndices = np.argsort(r)
        for j in range(r.shape[0]):
            print('name = ',allNames[sortedIndices[j]],': r = ',r[sortedIndices[j]])
        print('max(r) = ',np.max(r))

