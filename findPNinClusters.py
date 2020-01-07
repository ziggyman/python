import numpy as np
import os
import csvFree
import csvData

liuClusters = False
cantatClusters = False
lineOfSight = False

path = '/Users/azuri/daten/uni/HKU/PNe/star_clusters/gaia_ship/data/fof'
allFiles = os.path.join(path,'allClustersFiles.txt')
allPNeFileName = '/Users/azuri/daten/uni/HKU/PNe/star_clusters/all_true_PN.csv'
meanFilesName = os.path.join(path,'meanFiles.txt')
meanFileName = os.path.join(path,'meansSigmaRaDec.csv')
pnWithGaiaDistances = '/Users/azuri/daten/uni/HKU/PNe/inGaiaDr2/data.ssv'
cantatDataFile = '/Users/azuri/daten/uni/HKU/PNe/star_clusters/cantat-gaudin.ssv'

pnWithGaiaDistancesData = csvFree.readCSVFile(pnWithGaiaDistances,';',False)

#print('pnWithGaiaDistancesData.header = ',pnWithGaiaDistancesData.header)
#print('pnWithGaiaDistancesData.data = ',pnWithGaiaDistancesData.data)
sources = pnWithGaiaDistancesData.getData('Source')
#print('sources = ',sources)

nFound = 0

def writeMeans():
    meanFiles = []
    with open(allFiles,'r') as f:
        allFilesList = f.readlines()
    with open(meanFileName,'w') as mf:
        mf.write('scName,mean(ra),sigma(ra),mean(dec),sigma(dec)\n')
        for f in allFilesList:
            data = np.load(os.path.join(path,f.strip('\n')))
            ras = data['ra']
            decs = data['dec']
            scName = f[f.find('_')+1:f.find('.')]
            fNameOut = os.path.join(path,scName+'.csv')
            print('fNameOut = <'+fNameOut+'>')
            meanFiles.append(fNameOut)
            meanRA = np.mean(ras)
            sigmaRA = np.std(ras)
            meanDec = np.mean(decs)
            sigmaDec = np.std(decs)
            mf.write('%s,%.3f,%.3f,%.3f,%.3f\n' % (scName,meanRA,sigmaRA,meanDec,sigmaDec))
            with open(fNameOut, 'w') as fw:
                fw.write('mean(ra),sigma(ra),mean(dec),sigma(dec)\n')
                fw.write('%.3f,%.3f,%.3f,%.3f\n' % (meanRA,sigmaRA,meanDec,sigmaDec))
    with open(meanFilesName,'w') as fo:
        for f in meanFiles:
            fo.write(f+'\n')

def liuClustersMean():
    nFound = 0
    means = csvFree.readCSVFile(meanFileName,',',False)
    print('means.header = ',means.header)
    allPNe = csvFree.readCSVFile(allPNeFileName,',',False)
    print('allPNe.header = ',allPNe.header)
    for iPN in np.arange(0,allPNe.size(),1):
        for iSC in np.arange(0,means.size(),1):
            if ((float(allPNe.getData('DRAJ2000',iPN)) > float(means.getData('mean(ra)',iSC)) - (3.*float(means.getData('sigma(ra)',iSC))))
                and (float(allPNe.getData('DRAJ2000',iPN)) < float(means.getData('mean(ra)',iSC)) + (3.*float(means.getData('sigma(ra)',iSC))))
                and (float(allPNe.getData('DDECJ2000',iPN)) > float(means.getData('mean(dec)',iSC)) - (3.*float(means.getData('sigma(dec)',iSC))))
                and (float(allPNe.getData('DDECJ2000',iPN)) < float(means.getData('mean(dec)',iSC)) + (3.*float(means.getData('sigma(dec)',iSC))))):
                print(allPNe.getData('idPNMain',iPN),' found in ',means.getData('scName',iSC))
                nFound += 1
    return nFound

def liuClusters():
    with open(allFiles,'r') as f:
        allFilesList = f.readlines()
    #print('allFilesList = ',allFilesList)

    for source in sources:
        print('checking source id ',source)
        for f in allFilesList:
            data = np.load(os.path.join(path,f.strip('\n')))
    #        print('data[source_id] = ',data['source_id'])
            source_ids = data['source_id']
            for source_id in source_ids:
                print('source_id = ',source_id)
                if source == source_id:
                    print(source,' found')
                    nFound += 1

if cantatClusters:
    cantatData = csvFree.readCSVFile(cantatDataFile,';',False)
    source_ids = cantatData.getData('Source')
    for source in sources:
        for source_id in source_ids:
            if source_id == source:
                print('source_id = ',source_id,' found')
                nFound += 1

#writeMeans()
nFound = liuClustersMean()
print('nFound = ',nFound)
