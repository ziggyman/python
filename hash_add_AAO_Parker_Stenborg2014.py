import csv
import math
import os
from shutil import copyfile

from drUtils import getHeader#Value
from myUtils import degToHMS,degToDMS,getPNGName,raDecToLonLat,hmsToDeg,dmsToDeg,setHeaderKeyWord

path = '/Users/azuri/spectra/AAO/2014'
csvFileName = os.path.join(path,'AAOmega-Stenborg-Parker-service-run2014.csv')

hashFileIn = os.path.join(path,'hash_check.csv')
hashFileOut = os.path.join(path,'hash_found.csv')
sqlFileOut = os.path.join(path,'hash_add_SP.sql')

pnMainFile = '/Users/azuri/daten/uni/HKU/HASH/hash_PNMain_180221.csv'
tbCNamesFile = '/Users/azuri/daten/uni/HKU/HASH/hash_tbCNames_180221.csv'

rows = csv.DictReader(open(pnMainFile))
for row in rows:
    idPNMainStart = int(row['idPNMain'])+1
rows = csv.DictReader(open(tbCNamesFile))
for row in rows:
    idtbCNamesStart = int(row['idtbCNames'])+1

(_, _, filenames) = next(os.walk(path))
names = []
ids = ''
with open(hashFileIn,'w') as f:
    for filename in filenames:
#        print('filename[-5:] = ',filename[-5:])
        if (filename[-5:] == '.fits') and (filename[:3] == 'Fib'):
            print('filename = <'+filename+'>')
            fibreID = filename[3:]
            if fibreID[fibreID.rfind('.')-1] == 'B':
                fibreID = fibreID[:fibreID.find('B')]
            elif fibreID[fibreID.rfind('.')-1] == 'C':
                fibreID = fibreID[:fibreID.find('R')]
            rows = csv.DictReader(open(csvFileName))
            print('fibreID = ',fibreID)
            for row in rows:
                #print('row[#Fibre] = ',row['#Fibre'])
                if row['#Fibre'] == fibreID:
                    ra = math.degrees(float(row['RA']))
                    dec = math.degrees(float(row['DEC']))
                    raHMS = degToHMS(ra)
                    decDMS = degToDMS(dec)
                    print('ra = ',ra,', raHMS = ',raHMS)
                    print('dec = ',dec,', decDMS = ',decDMS)
                    if 'B' in filename:
                        arm = 'B'
                    else:
                        arm = 'R'
                    newFileName = os.path.join(path,row['Name']+arm+'_2D180614.fits')
                    copyfile(os.path.join(path,filename),newFileName)
                    setHeaderKeyWord(newFileName, 'RA', value=raHMS, hdu=0)
                    setHeaderKeyWord(newFileName, 'DEC', value=decDMS, hdu=0)
                    setHeaderKeyWord(newFileName, 'DATE-OBS', value='18-06-2014', hdu=0)
                    print(newFileName+': set RA to ',getHeader(newFileName)['RA'],' and DEC to ',getHeader(newFileName)['DEC'])
#                    print(newFileName+': set RA to ',getHeaderValue(newFileName, 'RA', hduNum=0),' and DEC to ',getHeaderValue(newFileName, 'DEC', hduNum=0))
                    if row['Name'] not in names:
                        f.write(row['Name']+','+raHMS+','+decDMS+'\n')
                        names.append(row['Name'])

with open(sqlFileOut,'w') as f:
    rows = csv.DictReader(open(hashFileOut))
    for row in rows:
        if row['pndb'] == '':
            print(row['id']+": new PN")
            rowsRaDec = csv.DictReader(open(os.path.join(path,'raAndDecs.csv')))
            ra = ''
            dec = ''
            for rowRaDec in rowsRaDec:
                if row['id'] == rowRaDec['id']:
                    ra = rowRaDec['RA']
                    dec = rowRaDec['DEC']
            lon,lat = raDecToLonLat(hmsToDeg(ra),dmsToDeg(dec))
            f.write("USE `MainGPN`;\n")
            f.write("INSERT INTO `PNMain`(`idPNMain`,`PNG`,`refPNG`,`RAJ2000`,`DECJ2000`,`DRAJ2000`,`DDECJ2000`,`Glon`,`Glat`,`refCoord`,`Catalogue`,`refCatalogue`,`userRecord`,`domain`,`refDomain`,`PNstat`,`refPNstat`,`refSimbadID`,`show`)")
            f.write("VALUES (%d,'%s','%s','%s','%s',%.5f,%.5f,%.5f,%.5f,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s');\n"
                    % (idPNMainStart,
                    getPNGName(lon,lat),
                    'sys',
                    ra,
                    dec,
                    hmsToDeg(ra),
                    dmsToDeg(dec),
                    lon,
                    lat,
                    'qparker',
                    'tstenborg',
                    'ziggy',
                    'ziggy',
                    'Galaxy',
                    'ziggy',
                    'c',
                    'ziggy',
                    'sys',
                    'y'
                    ))
            f.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`)")
            f.write("VALUES (%d,'%s','%s',%d,'%s','%s',%d);\n"
                    % (idtbCNamesStart,
                    row['id'][:2]+' '+row['id'][2:],
                    'qparker',
                    1,
                    'qparker',
                    'ziggy',
                    idPNMainStart
                    ))
            f.write("INSERT INTO `PNMain_tbCNames`(`idPNMain`,`idtbCNames`)")
            f.write("VALUES (%d,%d);\n"
                    % (idPNMainStart,
                    idtbCNamesStart
                    ))
            ids += str(idPNMainStart)+','
            idPNMainStart += 1
            idtbCNamesStart += 1
        else:
            print(row['id']+": "+row['pndb'])
            nameFound = False
            tbCNames = csv.DictReader(open(tbCNamesFile))
            for tbCName in tbCNames:
                if tbCName['Name'] == row['id'][:2]+' '+row['id'][2:]:
                    nameFound = True
            if not nameFound:
                if row['id'][:2] != 'PN':
                    f.write("INSERT INTO `tbCNames`(`idtbCNames`,`Name`,`reference`,`InUse`,`refInUse`,`userRecord`,`idPNMain`)")
                    f.write("VALUES (%d,'%s','%s',%d,'%s','%s',%d);\n"
                            % (idtbCNamesStart,
                            row['id'][:2]+' '+row['id'][2:],
                            'qparker',
                            0,
                            'qparker',
                            'ziggy',
                            int(row['pndb']),
                            ))
                    f.write("INSERT INTO `PNMain_tbCNames`(`idPNMain`,`idtbCNames`)")
                    f.write("VALUES (%d,%d);\n"
                            % (int(row['pndb']),
                            idtbCNamesStart
                            ))
                    idtbCNamesStart += 1
            ids += row['pndb']+','

print('ids = ',ids)
