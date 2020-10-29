import csv#Free,csvData

pnMainFile = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_231020.csv'
fitsFile = '/Users/azuri/daten/uni/HKU/HASH/hash_fitsfiles_231020.csv'
namesFile = '/Users/azuri/daten/uni/HKU/HASH/hash_tbCNames_231020.csv'

outFile = '/Users/azuri/daten/uni/HKU/HASH/add_spectra.sql'
catalogName = 'FrenchAmateurs'
catalogFile = '/Users/azuri/daten/uni/HKU/HASH/'+catalogName+'cat.sql'

#pnMain = csvFree.readCSVFile(pnMainFile)
#fits = csvFree.readCSVFile(fitsFile)
fitsStartId = 10141
reference = catalogName
instrument = 'Alpy 600 - 23 micron slit'
telescope = 'Newton Skywatcher 200 mm F/5'
sites = ['inner_CR','outer_CR','CN','CR']#['SA']

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
                            name = name[:name.find('_'+site)].replace('_','').replace('-','')
                    print(rowFits['idFitsFiles'],': name = <'+name+'>, bytes = ',bytes(name,'utf-8'))
                    found = False
                    namesReader = csv.DictReader(open(namesFile))
                    for rowNames in namesReader:
                        thisName = rowNames['Name'].replace('_','').replace('-','').replace(' ','')
        #                if rowNames['idPNMain'] == '2899':
        #                    print('thisName = <'+thisName+'>, bytes = ',bytes(thisName,'utf-8'))
                        if thisName == name:
                            found = True
                            idPNMain = rowNames['idPNMain']
                            print('found ',name,' in idPNMain = ',idPNMain)
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
                    if not found:
                        print('Problem: could not find name <'+name+'>')
                        STOP
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

if __name__ == "__main__":
    ids = addSpectraInfo()
    justMakeCatalog()
    print("don't forget to add entry to MainPNData.DataInfo")