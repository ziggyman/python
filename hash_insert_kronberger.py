import csv

inputFile = '/Users/azuri/daten/uni/HKU/HASH/DSH_PNe_20160303  (updated version_  Aug. 7, 2020).csv'

idPNMain_start = 32709
idAngDiam_start = 39407
idtbCNames_start = 121596

def getPNGName(lon,lat):
    png = '%010.6f' % lon
    png = png[:png.find('.')+2]
    #png = png.zfill(3)
    if lat > 0:
        png = png+'+'
    png = png + '%08.6g' % lat
    png = png[:png.rfind('.')+2]
    print('lon = ',lon,', lat = ',lat,', png = <'+png+'>')
    return png

with open('/Users/azuri/daten/uni/HKU/HASH/ingestDSH.sql','w') as ingest:
    with open('/Users/azuri/daten/uni/HKU/HASH/ingestDSHNames.sql','w') as ingestNames:
        with open('/Users/azuri/daten/uni/HKU/HASH/createDSHcatalogue.sql','w') as cat:
            with open('/Users/azuri/daten/uni/HKU/HASH/ingestDSH.hash','w') as hash:
                hash.write('hashpn fetch ')
                cat.write('CREATE TABLE IF NOT EXISTS MainPNData.DSHPNe_Aug2020 (\n')
                cat.write('idDSHPNe_Aug2020 INT AUTO_INCREMENT PRIMARY KEY UNIQUE,\n')
                cat.write('idPNMain INT,\n')
                cat.write('PNG VARCHAR(25),\n')
                cat.write('Name VARCHAR(25),\n')
                cat.write('PNstat VARCHAR(20),\n')
                cat.write('RAJ2000 VARCHAR(15),\n')
                cat.write('DECJ2000 VARCHAR(15),\n')
                cat.write('DRAJ2000 FLOAT,\n')
                cat.write('DDECJ2000 FLOAT,\n')
                cat.write('Glon FLOAT,\n')
                cat.write('Glat FLOAT,\n')
                cat.write('Catalogue VARCHAR(25),\n')
                cat.write('MajDiam VARCHAR(25),\n')
                cat.write('MinDiam VARCHAR(25),\n')
                cat.write('Comments VARCHAR(255),\n')
                cat.write('mapFlag VARCHAR(1) NOT NULL\n')
                cat.write(');\n')

                cat.write("USE `MainPNData`;\n")
                ingest.write("USE `MainGPN`;\n")
                ingestNames.write("USE `MainGPN`;\n")

                dsh = csv.DictReader(open(inputFile))
                print('dir(dsh) = ',dir(dsh))
                print('dsh.fieldnames() = ',dsh.fieldnames)
                j=0
                for rowDSH in dsh:
                    j+=1
                    idPNMain = idPNMain_start+j
                    hash.write(str(idPNMain)+',')
                    cat.write("INSERT INTO `DSHPNe_Aug2020`(`idDSHPNe_Aug2020`,`idPNMain`,`PNG`,`Name`,`PNstat`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDECJ2000`,`Glon`,`Glat`,`Catalogue`,`MajDiam`,`MinDiam`,`Comments`,`mapflag`) ")
                    print("j = ",j,": Comments1 = <",type(rowDSH['Comments1']),">, Comments2 = <",type(rowDSH['Comments2']),">, Comments3 = <",type(rowDSH['Comments3']),">, Comments4 = <",type(rowDSH['Comments4']),">")
                    print("j = ",j,": Comments1 = <",rowDSH['Comments1'],">, Comments2 = <",rowDSH['Comments2'],">, Comments3 = <",rowDSH['Comments3'],">, Comments4 = <",rowDSH['Comments4'],">")
                    comments = str(rowDSH['Comments1'])
                    if str(rowDSH['Comments2']) != '':
                        comments += '; '+str(rowDSH['Comments2'])
                    if str(rowDSH['Comments3']) != '':
                        comments += '; '+str(rowDSH['Comments3'])
                    if str(rowDSH['Comments4']) != '':
                        comments += '; '+str(rowDSH['Comments4'])
                    cat.write("VALUES (%d,%d,'%s','%s','%s','%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s');\n" % (j,
                                                                                                                           idPNMain,
                                                                                                                           rowDSH['PNG'],
                                                                                                                           rowDSH['Name'],
                                                                                                                           rowDSH['PNstat'],
                                                                                                                           rowDSH['RAJ2000'],
                                                                                                                           rowDSH['DECJ2000'],
                                                                                                                           float(rowDSH['DRAJ2000']),
                                                                                                                           float(rowDSH['DDECJ2000']),
                                                                                                                           float(rowDSH['Glon']),
                                                                                                                           float(rowDSH['Glat']),
                                                                                                                           rowDSH['Catalogue'],
                                                                                                                           rowDSH['MajDiam'],
                                                                                                                           rowDSH['MinDiam'],
                                                                                                                           comments,
                                                                                                                           'y'))
                    ingest.write("INSERT INTO `PNMain`(`idPNMain`,`PNG`,`refPNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDECJ2000`,`Glon`,`Glat`,`refCoord`,`Catalogue`,`refCatalogue`,`userRecord`,`domain`,`refDomain`,`PNstat`,`refPNstat`,`refSimbadID`,`show`) ")
                    ingest.write("VALUES (%d,'%s','%s', '%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n" % (idPNMain,
                                                                                                                                          getPNGName(float(rowDSH['Glon']),
                                                                                                                                                     float(rowDSH['Glat'])),
                                                                                                                                          'sys',
                                                                                                                                          rowDSH['RAJ2000'],
                                                                                                                                          rowDSH['DECJ2000'],
                                                                                                                                          float(rowDSH['DRAJ2000']),
                                                                                                                                          float(rowDSH['DDECJ2000']),
                                                                                                                                          float(rowDSH['Glon']),
                                                                                                                                          float(rowDSH['Glat']),
                                                                                                                                          'DSHPNe_Aug2020',
                                                                                                                                          'DSHPNe_Aug2020',
                                                                                                                                          'ziggy',
                                                                                                                                          'ziggy',
                                                                                                                                          'Galaxy',
                                                                                                                                          'ziggy',
                                                                                                                                          rowDSH['PNstat'],
                                                                                                                                          'ziggy',
                                                                                                                                          'sys',
                                                                                                                                          'y'))
                    if rowDSH['MajDiam'] != '':

                        if rowDSH['MinDiam'] != '':
                            ingest.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`MinDiam`,`reference`,`refTable`,`refRecord`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`tempflag`) ")
                            ingest.write("VALUES (%d,%d,%d,'%s','%s',%d,%d,'%s','%s',%d,'%s');\n" % (idAngDiam_start+j,
                                                                                                     int(rowDSH['MajDiam']),
                                                                                                     int(rowDSH['MinDiam']),
                                                                                                     'DSHPNe_Aug2020',
                                                                                                     'DSHPNe_Aug2020',
                                                                                                     j,
                                                                                                     1,
                                                                                                     'ziggy',
                                                                                                     'ziggy',
                                                                                                     idPNMain,
                                                                                                     'n'))
                        else:
                            ingest.write("INSERT INTO `tbAngDiam`(`idtbAngDiam`,`MajDiam`,`reference`,`refTable`,`refRecord`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`tempflag`) ")
                            ingest.write("VALUES (%d,%d,'%s','%s',%d,%d,'%s','%s',%d,'%s');\n" % (idAngDiam_start+j,
                                                                                                    int(rowDSH['MajDiam']),
                                                                                                    'DSHPNe_Aug2020',
                                                                                                    'DSHPNe_Aug2020',
                                                                                                    j,
                                                                                                    1,
                                                                                                    'ziggy',
                                                                                                    'ziggy',
                                                                                                    idPNMain,
                                                                                                    'n'))
                        ingest.write("INSERT INTO `PNMain_tbAngDiam`(`idPNMain`,`idtbAngDiam`) ")
                        ingest.write("VALUES (%d,%d);\n" % (idPNMain,
                                                            idAngDiam_start+j))

                    ingestNames.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`refTable`,`refRecord`,`InUse`,`refInUse`,`userRecord`,`idPNMain`,`comments`) ")
                    ingestNames.write("VALUES (%d,'%s','%s',%d,%d,'%s','%s',%d,'%s');\n" % (idtbCNames_start+j,
                                                                                       rowDSH['Name'],
                                                                                       'DSHPNe_Aug2020',
                                                                                       j,
                                                                                       1,
                                                                                       'ziggy',
                                                                                       'ziggy',
                                                                                       idPNMain,
                                                                                       'From DSH Kronberger list'))

                    ingestNames.write("INSERT INTO `PNMain_tbCNames`(`idPNMain`,`idtbCNames`) Values (%d,%d);\n" % (idPNMain, idtbCNames_start+j))

