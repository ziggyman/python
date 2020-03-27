import numpy as np
import os
from random import randint
import matplotlib.pyplot as plt
import scipy.stats as stats

import csvData,csvFree

path = '/Users/azuri/daten/uni/HKU/PN alignment'

def findPNeWithPAinHASH(inFileNamePAs, inFileNameHASH, outFileName):
    csvPAs = csvFree.readCSVFile(inFileNamePAs)
    print('csvPAs.header = ',csvPAs.header)

    csvHASH = csvFree.readCSVFile(inFileNameHASH)
    print('csvHASH.header = ',csvHASH.header)

    idsWithPA = csvPAs.getData('HASH ID')
    allIDs = csvHASH.getData('idPNMain')

    csvOut = csvData.CSVData()
    csvOut.header = csvHASH.header

    for iPA in range(len(idsWithPA)):
        found = False
        for iHASH in range(len(allIDs)):
            if idsWithPA[iPA] == allIDs[iHASH]:
                csvOut.append(csvHASH.getData(iHASH))
                found = True
        if not found:
            print('ERROR: ID ',idsWithPA[iPA],' not found in HASH data')
            STOP
    csvFree.writeCSVFile(csvOut,outFileName)

# @brief: change the GPAs for a certain percentage of a certain morphological class
#@param inFileName: string: name of input csvFile
#@param morphClassHeader: string: header string for the main morphological class
#@param morphClass: [string]: main morphological class for which to change the GPAs
#@param GPAHeader: string: header string for the GPA to change
#@param angleMean: [float]: mean GPA to adopt (degrees), same length as morphClass
#@param angleSDev: [float]: standard deviation around <angleMean> (degrees), same length as morphClass
#@param percentage: int or float: percentage of PNe with <morphClass> for which to change the GPA
#@param outFileName: string: output file name
def changePAs(inFileName,morphClassHeader,morphClass,GPAHeader,angleMean,angleSDev,percentage,outFileName):
    csv = csvFree.readCSVFile(inFileName)
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
        mu, sigma = 0., angleSDev[iM] / 180.
        dist = stats.truncnorm((a - mu) / sigma, (b - mu) / sigma, loc=mu, scale=sigma)
        randAngles = angleMean[iM] + (dist.rvs(nPNeToChangeGPA)*180.)
        for iR in range(len(randAngles)):
            if randAngles[iR] < 0.:
                randAngles[iR] += 360.
        plt.hist(randAngles)
        plt.show()
        randAngles = np.array(randAngles)
        print('mean(randAngles) = ',randAngles.mean(),', sDev(randAngles) = ',randAngles.std())

        for i in range(len(randomIDs)):
            csv.setData(GPAHeader,pNeWithMorphClass[randomIDs[i]],str(randAngles[i]))
    csvFree.writeCSVFile(csv,outFileName)


if __name__ == '__main__':
    paFileName = os.path.join(path, 'PN-alignments.csv')
    hashFileName = os.path.join(path, 'HASH_bipolar+elliptical_true_PNe.csv')

    findPNeWithPAinHASH(paFileName,hashFileName,hashFileName[:-4]+'_withPA.csv')

    hashFileName = hashFileName[:-4]+'_withPA.csv'
    hashOutFileName = hashFileName[:-4]+'_mock_Bmean%d_sdev%d_Emean%d_sdev%d.csv'
    mean = [45.,135]
    sdev = [10.,20.]
    changePAs(hashFileName,'mainClass',['B','E'],'GPA',mean,sdev,100,hashOutFileName%(int(mean[0]),int(sdev[0]),int(mean[1]),int(sdev[1])))
