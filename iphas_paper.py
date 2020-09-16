import numpy as np
from astroquery.vo_conesearch.conesearch import conesearch

import csvData
import csvFree
from myUtils import hmsToDeg, dmsToDeg

def coneSearch(csv):
    searches = []
    notFound = []
    for i in np.arange(0,csv.size(),1):
        try:
            search = conesearch(center=(hmsToDeg(csv.getData(' RA ',i).strip()), dmsToDeg(csv.getData(' DEC ',i).strip())),
                                radius=0.001,#in degrees
                                verb=3,
                                catalog_db="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&")
            print('dir(search) = ',dir(search))
            print('search.nrows() = ',search.nrows)
            print('search.array() = ',search.array)
            print('search.fields() = ',search.fields)
            print('search.params() = ',search.params)
            print('search.to_table() = ',search.to_table())
#            search.to_table().write('/Volumes/work/azuri/spectra/IPHAS_GTC/iphas-data.fits', format='fits')
            searches.append([csv.getData('Name ',i),search])
        except Exception as e:
            print('Exception caught: ',e)
            print(csv.getData('Name ',i)," not found")
            notFound.append(csv.getData('Name ',i))
            pass
    print('searches = ',len(searches),': ',searches)
    print('not found: ',len(notFound),': ',notFound)

def findHASHid(csvTablePaper, csvTableTargets):
    ids = []
    for iObs in range(csvTablePaper.size()):
        name = csvTablePaper.getData("Name ",iObs).strip(' ')
        if name == 'IPHASX J014238+600947':
            name = 'We 2-5'
        elif name == 'IPHASX J031058.8+624755':
            name = 'Sh 2-200'
        elif name == 'Pa 21':
            name = 'DSH J192315+270734'
        elif name == 'Pa 27':
            name = 'DSH J204858+321815'
        elif name == 'Pa 22':
            name = 'DSH J195813+395440'
        elif name == 'Pa 29':
            name = 'DSH J205943+345423'
        elif name == 'Kn 62':
            name = 'DSH J062355+381515'
        elif name == 'Pa 15':
            name = 'DSH J202907+231108'
        elif name == 'IPHASX J191104.8+060845':
            name = 'IRAS 19086+0603'
#        elif name == ''
#        print('Name =<'+name+'>')
        found = False
        for iTarget in range(csvTableTargets.size()):
            targetName = csvTableTargets.getData('Name',iTarget)
#            print('checking targetName <'+targetName+'>')
            if name == targetName:
#                print('found it!')
                ids.append([name,csvTableTargets.getData('idPNMain', iTarget)])
                found = True
        if not found:
            print('PROBLEM: did not find object named <'+name+'>')
            print('object data: ',csvTablePaper.getData(iObs))
    if len(ids) == csvTablePaper.size():
        print('FOUND ALL TARGETS!')
    else:
        print('HMMM, NOT ALL TARGETS FOUND... :(')
    return ids

def getImages(ids):
    print('getImages: ids = ',len(ids),': ',ids)
    hashPNMain = '/Users/azuri/daten/uni/HKU/HASH/PNMain_full_Aug-07-2020.csv'

    outFile = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/getIPHASImages.bash'
    with open(outFile,'w') as f:
        f.write('mkdir iphas-gtc-images\n')
        for survey in [['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_iphas.csv','/data/kegs/pngPNImages/%s/IPHAS/%s_%s_iphas3colour*.png',['idPNMain','idPNMain','run_id']],
                        ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_iquotHaSr.csv','/data/kegs/pngPNImages/%s/IPHAS/%s_%s_iquot*.png',['idPNMain','idPNMain','run_id']],
                        ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_2mass.csv','/data/kegs/pngPNImages/%s/TWOMASS/%s_*.png',['idPNMain','idPNMain']],
                        ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_shs.csv','/data/kegs/pngPNImages/%s/SHS/%s_thre*.png',['idPNMain','idPNMain']],
                        ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_quotHaSr.csv','/data/kegs/pngPNImages/%s/SHS/%s_quot*.png',['idPNMain','idPNMain']],
                        ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_vphas.csv','/data/kegs/pngPNImages/%s/VPHASplus/%s_vp*.png',['idPNMain','idPNMain']],
                        ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_vquotHaSr.csv','/data/fermenter/PNImages/%s/VPHASplus/quot_vphas*',['idPNMain']],
                        ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_wise.csv','/data/kegs/pngPNImages/%s/WISE/%s_wise321rgb*.png',['idPNMain','idPNMain']]]:
            csv = csvFree.readCSVFile(survey[0])
            for id in ids:
                hashID = id[1]
                lineIDs = csv.find('idPNMain',hashID)
                print('getImages: hashID = ',hashID,': lineIDs = ',lineIDs)
                if lineIDs[0] == -1:
                    print('getImages: PROBLEM: did not find hashID ',hashID,' in ',survey[0])
                    #STOP
                else:
                    for lineID in lineIDs:
                        if csv.getData('inuse',lineID) == '1':
                            print('survey[1] = ',survey[1],' survey[2] = ',survey[2])
                            tmp = [csv.getData(num,lineID) for num in survey[2]]
                            print('tmp = ',tmp)
                            print('dir(tmp) = ',dir(tmp))
                            print('(num for num in survey[2]) = ',(num for num in tmp))
                            print('survey[1] = ',survey[1],' % (',[csv.getData(num,lineID) for num in survey[2]],')')
                            f.write('cp '+survey[1] % [csv.getData(num,lineID) for num in survey[2]]+' iphas-gtc-images/\n')

if __name__ == '__main__':
    csvPaper = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table_paper.csv','&',False)
#    print(csvPaper.header)
#    print(csvPaper.data)

    csvTargets = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/Charts_IPHAS_Amateurs/Targets_GTC.cvs')
#    print(csvTargets.header)
#    print(csvTargets.data)
    #    coneSearch(csvPaper)

    ids = findHASHid(csvPaper, csvTargets)
    ids.append(['Ou 1','8458'])
    ids.append(['IPHASX J055242.8+262116','9824'])
    ids.append(['J190333','8506'])
    ids.append(['IPHASX J191707.3+020010','2502'])
    print(ids)
    getImages(ids)



