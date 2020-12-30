import os
import csv#Free,csvData
from myUtils import getHeader

pnMainFile = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_241220.csv'
fitsFile = '/Users/azuri/daten/uni/HKU/HASH/hash_fitsfiles_241220.csv'
namesFile = '/Users/azuri/daten/uni/HKU/HASH/hash_tbCNames_241220.csv'

inputSpectraDir = '/Users/azuri/spectra/saao/saao_mar2014/hash/'
searchFile = '/Users/azuri/daten/uni/HKU/HASH/hash_search.csv'
searchFileOut = '/Users/azuri/daten/uni/HKU/HASH/hash_found.csv'

outFile = '/Users/azuri/daten/uni/HKU/HASH/add_spectra.sql'
catalogName = 'SAAO_Mar2014'
catalogFile = '/Users/azuri/daten/uni/HKU/HASH/'+catalogName+'cat.sql'
hashpnFile = '/Users/azuri/daten/uni/HKU/HASH/'+catalogName+'_hashpn.sql'

#pnMain = csvFree.readCSVFile(pnMainFile)
#fits = csvFree.readCSVFile(fitsFile)
fitsStartId = 10221
reference = catalogName
instrument = 'SpUpNIC SIT1 1CD'
telescope = 'SAAO 1.9m'
sites = ['SA']#['inner_CR','outer_CR','CN','CR']#

def createHASHSearchInput():
    (_, _, filenames) = next(os.walk(inputSpectraDir))
    with open(searchFile,'w') as f:
        for fName in filenames:
            print('fName = <'+fName+'>')
            header = getHeader(os.path.join(inputSpectraDir,fName),0)
            f.write('%s,%s,%s\n' % (fName,header['RA'].strip(),header['DEC'].strip()))

def getIDs():
    reader = csv.DictReader(open(searchFileOut))
    outStr = ''
    for row in reader:
        outStr += row['pndb']+','
    print(outStr)

def addSpectraInfo():
    with open(outFile,'w') as f:
        f.write('USE PNSpectra_Sources;\n')
        with open(catalogFile,'w') as fb:
            fb.write('CREATE TABLE IF NOT EXISTS MainPNData.'+catalogName+' (\n')
            fb.write('id'+catalogName+' INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
            fb.write('idPNMain INT,\n')
            fb.write('mapFlag VARCHAR(1) NOT NULL\n')
            fb.write(');\n')

            fb.write("USE `MainPNData`;\n")

            fitsReader = csv.DictReader(open(fitsFile))
            ids = []
            nFound = 1
            for rowFits in fitsReader:
        #        print("int(rowFits['idFitsFiles']) = ",int(rowFits['idFitsFiles']))
                if (int(rowFits['idFitsFiles']) >= fitsStartId) and (rowFits['reference'] == reference):
                    name = rowFits['fileName']
                    for site in sites:
                        if name.find('_'+site) >= 0:
                            name = name[:name.find('_'+site)].replace('_','')
                            if name[:3] != 'PNG':
                                name = name.replace('-','')
                    if name[:-1] == 'NGC246':
                        name = 'NGC246'
                    elif name == 'RWT152above':
                        name = 'RWT152'
                    elif name == 'RWT152below':
                        name = 'RWT152'
                    elif name == 'SNRJOE':
                        name = 'SNRJoe'
                    elif name == 'SNRJOE1':
                        name = 'SNRJoe'
                    elif name == 'SNRJoe2':
                        name = 'SNRJoe'
                    elif name == 'SNRJoe3':
                        name = 'SNRJoe'
                    elif name == 'SNR15761':
                        name = 'SNR1576'
                    elif name == 'SNR15762':
                        name = 'SNR1576'
                    elif name == 'SNR15764':
                        name = 'SNR1576'
                    elif name == 'ha1761':
                        name = 'ha176'
                    elif name == 'SP17322807C':
                        name = 'SP17322807'
                    print(rowFits['idFitsFiles'],': name = <'+name+'>, bytes = ',bytes(name,'utf-8'))
                    found = False
                    if name[:3] != 'PNG':
                        namesReader = csv.DictReader(open(namesFile))
                        for rowNames in namesReader:
                            thisName = rowNames['Name'].replace('_','').replace('-','').replace(' ','')
            #                if rowNames['idPNMain'] == '2899':
            #                    print('thisName = <'+thisName+'>, bytes = ',bytes(thisName,'utf-8'))
                            if thisName == name:
                                found = True
                                idPNMain = rowNames['idPNMain']
                                print('found ',name,' in idPNMain = ',idPNMain,', rowNames[Name] = ',rowNames['Name'])
                    else:
                        pnMainReader = csv.DictReader(open(pnMainFile))
                        tempName = name[3:]
                        for rowMain in pnMainReader:
                            if rowMain['PNG'] == tempName:
                                found = True
                                idPNMain = rowMain['idPNMain']
                    if not found:
                        print("ERROR: name "+name+' not found')
                        STOP
                    ids.append([rowFits['idFitsFiles'],idPNMain])

                    pnMainReader = csv.DictReader(open(pnMainFile))
                    for rowMain in pnMainReader:
                        if rowMain['idPNMain'] == idPNMain:
                            f.write("UPDATE `FitsFiles` SET `reference` = '%s', `idPNMain` = %d, `object` = '%s', `instrument` = '%s', `telescope` = '%s', `RAJ2000` = '%s', `DECJ2000` = '%s', `DRAJ2000` = %.5f, `DDECJ2000` = %.5f, `convToText` = 'y' WHERE `idFitsFiles` = %d;\n"
                                     % (reference, int(idPNMain), rowNames['Name'], instrument, telescope, rowMain['RAJ2000'], rowMain['DECJ2000'], float(rowMain['DRAJ2000']), float(rowMain['DDECJ2000']), int(rowFits['idFitsFiles'])))
                            fb.write("INSERT INTO `"+catalogName+"`(`idPNMain`,`mapflag`)")
                            fb.write("VALUES (%d,'%s');\n" % (int(idPNMain),
                                                              'y'))
                    nFound+=1
#                    if not found:
#                        print('Problem: could not find name <'+name+'>')
#                        STOP
    return ids


def justMakeCatalog():
    with open(catalogFile,'w') as fb:
        fb.write('CREATE TABLE IF NOT EXISTS MainPNData.'+catalogName+' (\n')
        fb.write('id'+catalogName+' INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
        fb.write('idPNMain INT,\n')
        fb.write('mapFlag VARCHAR(1) NOT NULL\n')
        fb.write(');\n')

        fb.write("USE `MainPNData`;\n")

        fitsReader = csv.DictReader(open(fitsFile))
        ids = []
        nFound = 0
        for rowFits in fitsReader:
    #        print("int(rowFits['idFitsFiles']) = ",int(rowFits['idFitsFiles']))
            if (int(rowFits['idFitsFiles']) >= fitsStartId) and (rowFits['reference'] == reference):
                idPNMain = rowFits['idPNMain']
                if idPNMain not in ids:
                    fb.write("INSERT INTO `"+catalogName+"`(`idPNMain`,`mapflag`)")
                    fb.write("VALUES (%d,'%s');\n" % (int(idPNMain),
                                                      'y'))
                else:
                    print('id ',idPNMain,' already in catalog')
                ids.append(idPNMain)
                nFound += 1
    with open(hashpnFile,'w') as hf:
        hf.write('hashpn spetch all ')
        hf.write(str(ids[0]))
        for id in ids[1:]:
            hf.write(','+str(id))
        hf.write(' -w')

if __name__ == "__main__":
#    createHASHSearchInput()
#    getIDs()
#    ids = addSpectraInfo()
    justMakeCatalog()
    print("don't forget to add entry to MainPNData.DataInfo")
