import os
import csv#Free,csvData
import shutil
from myUtils import getHeader, find_nth, setHeaderKeyWord

pnMainFile = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_190221.csv'
fitsFile = '/Users/azuri/daten/uni/HKU/HASH/hash_fitsfiles_190221.csv'
namesFile = '/Users/azuri/daten/uni/HKU/HASH/hash_tbCNames_190221.csv'

inputSpectraDir = '/Users/azuri/spectra/MSO/wifes_jul2011/'#'/Users/azuri/daten/uni/HKU/HASH/FrenchAmateurs/new/Jan2021a/hash/'#
fibreFile = None#'556Rcomb.fibres.lis'
doneFile = '/Users/azuri/spectra/mash_spectra.txt'#None
searchFile = '/Users/azuri/daten/uni/HKU/HASH/hash_search.csv'
searchFileOut = '/Users/azuri/daten/uni/HKU/HASH/hash_found.csv'

outFile = '/Users/azuri/daten/uni/HKU/HASH/add_spectra.sql'
catalogName = 'WiFeS_Jul2011'#'FrenchAmateurs_Jan2021'#'MASH_REJSPEC_Jan2021'
setName = 'WiFeS_Jul2011'
catalogFile = '/Users/azuri/daten/uni/HKU/HASH/'+catalogName+'cat.sql'
hashpnFile = '/Users/azuri/daten/uni/HKU/HASH/'+catalogName+'_hashpn.txt'

#pnMain = csvFree.readCSVFile(pnMainFile)
#fits = csvFree.readCSVFile(fitsFile)
fitsStartId = 10625
reference = catalogName#'FrenchAmateurs'#
instrument = 'WiFeS'#'DBS'#
telescope = 'MSSO'#'AAO 2.3m'#
sites = ['MS']#['CN','PO','KO','CO']#['inner_CR','outer_CR','CN','CR']#

def getHashIDsFromNames():
    (_, _, filenames) = next(os.walk(inputSpectraDir))
    ids = []
    for fName in filenames:
        if fName[fName.rfind('.')+1:] == 'fits':
            name = fName[:fName.rfind('_')]
            if name[:2] == 'PK':
                name = name[:2]+'00'+name[2:]
            id = getHashIDfromName(name)
            if id == -1:
                print('ERROR: could not find name <'+name+'>')
            else:
                ids.append(id)
    strOut = ids[0]
    for id in ids[1:]:
        strOut += ','+id
    print(str)
    return ids

def createHASHSearchInput(fibreFile=None):
    if doneFile is not None:
        with open(doneFile,'r') as f:
            doneLines = f.readlines()
        doneLines = [doneLine.strip() for doneLine in doneLines]
        doneLines = [fixFileName(doneLine) for doneLine in doneLines]
    if fibreFile is None:
        (_, _, filenames) = next(os.walk(inputSpectraDir))
        with open(searchFile,'w') as f:
            for fName in filenames:
                if fName[fName.rfind('.')+1:] == 'fits':
                    print('fName = <'+fName+'>')
                    if (doneFile is not None) and (fName in doneLines):
                        print('found fName ',fName,' in doneFile')
                    else:
                        header = getHeader(os.path.join(inputSpectraDir,fName),0)
                        try:
                            print('fName = <'+fName+'>: RA = ',header['RA'],', DEC = ',header['DEC'])
                            f.write('%s,%s,%s\n' % (fName,header['RA'].strip(),header['DEC'].strip()))
                        except:
                            print('fName = <'+fName+'>: RA = ',header['MEANRA'],', DEC = ',header['MEANDEC'])
                            f.write('%s,%s,%s\n' % (fName,header['MEANRA'].strip(),header['MEANDEC'].strip()))
                        shutil.copyfile(os.path.join(inputSpectraDir,fName),os.path.join(os.path.join(inputSpectraDir,'hash'),fName))
    else:
        dictOut = []
        lines = None
        with open(fibreFile,'r') as f:
            lines = f.readlines()
        lines = [line.strip() for line in lines]
        with open(searchFile,'w') as f:
            for line in lines:
                if line[0] != '#':
                    name = line[3:32].replace(' ','')
                    if (name != 'Parked') and (not ('Sky' in name)):
                        if (doneFile is None):
                            f.write(name+','+line[40:52].strip(' ').replace(' ',':')+','+line[54:67].strip(' ').replace(' ',':')+'\n')
                        else:
                            """TODO: check for the same date as well"""
                            found = False
                            for doneLine in doneLines:
                                for site in sites:
                                    nameSite = name+'_'+site
                                    if nameSite in doneLine:
                                        found = True
                                        print('found name ',nameSite,' in files already done = ',)
                            if not found:
                                f.write(name+','+line[40:52].strip(' ').replace(' ',':')+','+line[54:67].strip(' ').replace(' ',':')+'\n')


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
                print("int(rowFits['idFitsFiles']) = ",int(rowFits['idFitsFiles']))
                if (int(rowFits['idFitsFiles']) >= fitsStartId) and (rowFits['setname'] == setName):
                    name = rowFits['fileName']
                    print('name = ',name)
                    if instrument != '':
                        instr = instrument
                    else:
                        header = getHeader(os.path.join(inputSpectraDir,name),0)
                        instr = header['BSS_INST'].strip()[find_nth(header['BSS_INST'],' ',2)+1:]
                        print('instr = ',instr)
                    if telescope != '':
                        tel = telescope
                    else:
                        header = getHeader(os.path.join(inputSpectraDir,name),0)
                        tel = header['BSS_INST'].strip()[:find_nth(header['BSS_INST'].strip(),' ',2)]
                        print('tel = ',tel)
                    for site in sites:
                        if name.find('_'+site) >= 0:
                            name = name[:name.find('_'+site)].replace('_','')
                            if name[:3] != 'PNG':
                                name = name.replace('-','')
                    if name[:-1] == 'NGC246':
                        name = 'NGC246'
                    elif name == 'PK29800B':
                        name = 'PK298001'
                    elif name == 'V432CAR':
                        name = 'V*V432Car'
                    elif name == 'V838MON':
                        name = 'V838Mon'
                    elif name == 'PHR072513473':
                        name = 'PHR07251347'
                    elif name == 'PHR143861403':
                        name = 'PHR14386140'
                    elif name == 'PK29800':
                        name = 'PK298001'
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
                    elif name == 'CRBB1NEB':
                        name = 'CRBB1'
                    elif name == 'SP17463413':
                        name = 'SPJ17463413'
                    elif name[:2] == 'PK':
                        name = name[1:]
                    elif name[:3] == 'HII':
                        name = 'PHR'+name[3:]
                    elif name[-1] == 'B':
                        name = name[:-1]
                    elif name[-1] == 'R':
                        name = name[:-1]
                    print(rowFits['idFitsFiles'],': name = <'+name+'>, bytes = ',bytes(name,'utf-8'))
                    found = False
                    object = ''
                    if name[:3] == 'PNG':
                        pnMainReader = csv.DictReader(open(pnMainFile))
                        tempName = name[3:]
                        for rowMain in pnMainReader:
                            if rowMain['PNG'] == tempName:
                                found = True
                                idPNMain = rowMain['idPNMain']
                                object = 'PNG '+tempName
                    elif name[:3] == 'PPA':
                        namesReader = csv.DictReader(open(namesFile))
                        for rowNames in namesReader:
                            thisName = rowNames['Name'].replace('_','').replace('-','').replace(' ','')
            #                if rowNames['idPNMain'] == '2899':
            #                    print('thisName = <'+thisName+'>, bytes = ',bytes(thisName,'utf-8'))
                            if thisName == name:
                                found = True
                                idPNMain = rowNames['idPNMain']
                                print('found ',name,' in idPNMain = ',idPNMain,', rowNames[Name] = ',rowNames['Name'])
                                object = rowNames['Name']
                        name = 'PHR'+name[3:]
                        namesReader = csv.DictReader(open(namesFile))
                        for rowNames in namesReader:
                            thisName = rowNames['Name'].replace('_','').replace('-','').replace(' ','')
            #                if rowNames['idPNMain'] == '2899':
            #                    print('thisName = <'+thisName+'>, bytes = ',bytes(thisName,'utf-8'))
                            if thisName == name:
                                found = True
                                idPNMain = rowNames['idPNMain']
                                print('found ',name,' in idPNMain = ',idPNMain,', rowNames[Name] = ',rowNames['Name'])
                                object = rowNames['Name']
                    else:
                        namesReader = csv.DictReader(open(namesFile))
                        for rowNames in namesReader:
                            thisName = rowNames['Name'].replace('_','').replace('-','').replace(' ','')
            #                if rowNames['idPNMain'] == '2899':
            #                    print('thisName = <'+thisName+'>, bytes = ',bytes(thisName,'utf-8'))
                            if thisName == name:
                                found = True
                                idPNMain = rowNames['idPNMain']
                                print('found ',name,' in idPNMain = ',idPNMain,', rowNames[Name] = ',rowNames['Name'])
                                object = rowNames['Name']
                    if not found:
                        print("ERROR: name "+name+' not found')
                        STOP
                    ids.append([rowFits['idFitsFiles'],idPNMain])

                    pnMainReader = csv.DictReader(open(pnMainFile))
                    for rowMain in pnMainReader:
                        if rowMain['idPNMain'] == idPNMain:
#                            print('rowNames["Name"] = <'+rowNames['Name']+'>')
                            f.write("UPDATE `FitsFiles` SET `reference` = '%s', `idPNMain` = %d, `object` = '%s', `instrument` = '%s', `telescope` = '%s', `RAJ2000` = '%s', `DECJ2000` = '%s', `DRAJ2000` = %.5f, `DDECJ2000` = %.5f, `convToText` = 'y' WHERE `idFitsFiles` = %d;\n"
                                     % (reference, int(idPNMain), object, instr, tel, rowMain['RAJ2000'], rowMain['DECJ2000'], float(rowMain['DRAJ2000']), float(rowMain['DDECJ2000']), int(rowFits['idFitsFiles'])))
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
            if (int(rowFits['idFitsFiles']) >= fitsStartId) and (rowFits['setname'] == setName):
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

def getHashIDfromName(name):
    nam = name.replace(' ','').replace('_','').replace('-','')
    tbCNames = csv.DictReader(open(namesFile))
    for line in tbCNames:
        n = line['Name'].replace(' ','').replace('_','').replace('-','')
        if n == nam:
            print('getHashIDfromName: name = <'+name+'> found with idPNMain = ',line['idPNMain'])
            return line['idPNMain']
    print('getHashIDfromName: name = <'+name+'> NOT FOUND in tbCNames')
    return -1

def getRaDec(idPNMain):
    pnMain = csv.DictReader(open(pnMainFile))
    for row in pnMain:
        if row['idPNMain'] == idPNMain:
            return (row['RAJ2000'],row['DECJ2000'])
    return ('0','0')

def fixRAandDEC():
    (_, _, filenames) = next(os.walk(inputSpectraDir))
    for fName in filenames:
        if fName[fName.rfind('.')+1:] == 'fits':
            fullName = os.path.join(inputSpectraDir,fName)
            header = getHeader(fullName,0)
            print('fName = <'+fName+'>: RA = ',header['RA'],', DEC = ',header['DEC'])
            idPNMain = -1
            if isinstance(header['RA'],float) or isinstance(header['DEC'],float):
                name = fName[:fName.find('_')].replace('rej','')
                idPNMain = getHashIDfromName(name)
                if idPNMain == -1:
                    idPNMain = getHashIDfromName(name[:3]+'J'+name[3:])
                    if idPNMain == -1:
                        print('moving to not_found/')
                        shutil.move(fullName,os.path.join(inputSpectraDir,'not_found'))
                if idPNMain != -1:
                    print('changing RA and DEC')
                    ra, dec = getRaDec(idPNMain)
                    print('fName = <'+fName+'>: name = <'+name+'>: idPNMain = ',idPNMain,': ra = ',ra,', dec = ',dec)
                    setHeaderKeyWord(fullName,'RA',ra,0)
                    setHeaderKeyWord(fullName,'DEC',dec,0)
            else:
                print('fName = <'+fName+'>: RA and DEC are okay')
            print(' ')

def fixFileName(fName):
    replace = [['_NF',''],
               ['NF',''],
               ['HD_','_'],
               ['HD',''],
               ['_SA_','_SA'],
               ['_paperclipSNR',''],
               ['-probably',''],
               ['N_','_'],
               ['_S0','_SA0'],
               ['_S1','_SA1'],
               ['_S2','_SA2'],
               ['_S3','_SA3'],
               ['SF_','_'],
              ]
    fNameNew = fName
    for rep in replace:
        fNameNew = fNameNew.replace(rep[0],rep[1])
    print('fName = ',fName,', fNameNew = ',fNameNew)
    return fNameNew

def fixFileNames():
    (_, _, filenames) = next(os.walk(inputSpectraDir))
    for fName in filenames:
        if fName[fName.rfind('.'):] == '.fits':
            fNameNew = fixFileName(fName)
            if fName != fNameNew:
                os.rename(os.path.join(inputSpectraDir,fName),os.path.join(inputSpectraDir,fNameNew))

if __name__ == "__main__":
#    getHashIDsFromNames()
#    fixRAandDEC()
#    fixFileNames()
#    createHASHSearchInput()
#    createHASHSearchInput(os.path.join(inputSpectraDir,fibreFile))
#    getIDs()
#    ids = getHashIDsFromNames()
#    print(ids)
    ids = addSpectraInfo()
#    justMakeCatalog()
#    print("don't forget to add entry to MainPNData.DataInfo")
