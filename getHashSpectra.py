import csv
import os

pnMainFile = '/Users/azuri/daten/uni/HKU/HASH/PNMain_full_Aug-07-2020.csv'
pnMainFitsFile = '/Users/azuri/daten/uni/HKU/HASH/PNSpectra_Sources_FitsFiles_full_Aug-07-2020.csv'
pnMainSpectraInfoFile = '/Users/azuri/daten/uni/HKU/HASH/PNSpectra_Sources_spectraInfo_full_Aug-07-2020.csv'

pnMain = csv.DictReader(open(pnMainFile))

outfile = '/Users/azuri/daten/uni/HKU/HASH/copyFits.sh'

with open(outfile,'w') as f:
    for rowPnMain in pnMain:
        if rowPnMain['PNstat'] == 'T':
            idPNMain = rowPnMain['idPNMain']
            print('idPNMain = ',idPNMain,' is a True PN')
            pnMainFits = csv.DictReader(open(pnMainFitsFile))
            for rowPnMainFits in pnMainFits:
                if rowPnMainFits['idPNMain'] == idPNMain:
                    setName = rowPnMainFits['setname']
                    fileName = rowPnMainFits['fileName']
                    print('idPNMain = ',idPNMain,': spectrum found in setName = <'+setName+'>, fileName = <'+fileName+'>')
                    pnMainSpectraInfo = csv.DictReader(open(pnMainSpectraInfoFile))
                    for rowPnMainSpectraInfo in pnMainSpectraInfo:
                        if rowPnMainSpectraInfo['name'] == setName:
                            path = rowPnMainSpectraInfo['path']
                            print('idPNMain = ',idPNMain,': path = <'+path+'>')
                            f.write('cp -p '+os.path.join('/data/mashtun/',path,fileName)+' simran/\n')
    f.write('tar -zcvf simran.tgz simran')