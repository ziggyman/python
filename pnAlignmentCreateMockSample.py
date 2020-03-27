import numpy as np
import os
from random import randint
import matplotlib.pyplot as plt
import scipy.stats as stats

import csvData,csvFree
from pnAlignment import findPNeWithPAinHASH

path = '/Users/azuri/daten/uni/HKU/PN alignment'

# @brief: change the GPAs for a certain percentage of a certain morphological class
#@param inFileName: string: name of input csvFile
#@param morphClassHeader: string: header string for the main morphological class
#@param morphClass: [string]: main morphological class for which to change the GPAs
#@param GPAHeader: string: header string for the GPA to change
#@param angleMean: [float]: mean GPA to adopt (degrees), same length as morphClass
#@param angleSDev: [float]: standard deviation around <angleMean> (degrees), same length as morphClass
#@param percentage: int or float: percentage of PNe with <morphClass> for which to change the GPA
#@param outFileName: string: output file name
def changePAs(inFileNameData,morphClassHeader,morphClass,GPAHeader,angleMean,angleSDev,percentage,outFileName):
    csv = csvFree.readCSVFile(inFileNameData)
    morphClasses = csv.getData(morphClassHeader)

    for iM in range(len(morphClass)):
        nPNeWithMorphClass = 0
        pNeWithMorphClass = []
        for i in range(csv.size()):
            if morphClasses[i] == morphClass[iM]:
                pNeWithMorphClass.append(i)
                nPNeWithMorphClass += 1

        randomIDs = []
        nPNeToChangeGPA = int(nPNeWithMorphClass * percentage / 100)
        while len(randomIDs) < nPNeToChangeGPA:
            rand = randint(0,nPNeToChangeGPA-1)
            if rand not in randomIDs:
                randomIDs.append(rand)
        print('nPNeToChangeGPA = ',nPNeToChangeGPA)
        print('randomIDs = ',len(randomIDs),': ',randomIDs)

        # for the circular problem rescale sigma and use random range [-1,1]
        a, b = -1., 1.
        mu, sigma = 0., angleSDev[iM] / 90.
        dist = stats.truncnorm((a - mu) / sigma, (b - mu) / sigma, loc=mu, scale=sigma)
        randAngles = angleMean[iM] + (dist.rvs(nPNeToChangeGPA)*90.)
        for iR in range(len(randAngles)):
            if randAngles[iR] < 0.:
                randAngles[iR] += 180.
            if randAngles[iR] >= 180.:
                randAngles[iR] -= 180.
        plt.hist(randAngles)
        plt.show()
        randAngles = np.array(randAngles)
        print('mean(randAngles) = ',randAngles.mean(),', sDev(randAngles) = ',randAngles.std())

        for i in range(len(randomIDs)):
            print('i = ',i,': HASH ID = ',csv.getData("idPNMain",pNeWithMorphClass[randomIDs[i]]),', mainClass=',csv.getData(morphClassHeader, pNeWithMorphClass[randomIDs[i]]),', GPA = ',csv.getData(GPAHeader,pNeWithMorphClass[randomIDs[i]]))
            csv.setData(GPAHeader,pNeWithMorphClass[randomIDs[i]],str(randAngles[i]))
            print('i = ',i,': HASH ID = ',csv.getData("idPNMain",pNeWithMorphClass[randomIDs[i]]),', mainClass=',csv.getData(morphClassHeader, pNeWithMorphClass[randomIDs[i]]),', GPA = ',csv.getData(GPAHeader,pNeWithMorphClass[randomIDs[i]]))
    csvFree.writeCSVFile(csv,outFileName)


if __name__ == '__main__':
    paFileName = os.path.join(path, 'PN-alignments.csv')
    hashFileName = os.path.join(path, 'HASH_bipolar+elliptical_true_PNe.csv')

    newHashFileName = hashFileName[:-4]+'_withPA.csv'
    findPNeWithPAinHASH(paFileName,hashFileName,newHashFileName)

    hashFileName = newHashFileName
    hashOutFileName = hashFileName[:hashFileName.rfind('/')]+'/mock'+hashFileName[hashFileName.rfind('/'):-4]
    print('hashOutFileName = <'+hashOutFileName+'>')

    mainClasses = ['B','E']
    mean = [15.,135]
    sdev = [10.,20.]
    for i in range(len(mainClasses)):
        hashOutFileName = hashOutFileName+'_'+mainClasses[i]+'mean%d_sdev%d' % (mean[i],sdev[i])
        print('hashOutFileName = <'+hashOutFileName+'>')

    hashOutFileName += '.csv'
    print('hashOutFileName = <'+hashOutFileName+'>')

    changePAs(hashFileName,'mainClass',['B','E'],'GPA',mean,sdev,100,hashOutFileName)
    print('finished writing <'+hashOutFileName+'>')
