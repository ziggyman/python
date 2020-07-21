import csv
from myUtils import hmsToDeg, dmsToDeg, raDecToLonLat, lonLatToRaDec, angularDistancePyAsl, angularDistance

fNameHashPNMain = '/Users/azuri/daten/uni/HKU/HASH/PNMain_full_July-10-2020.csv'
#fNameHashPNMain = '/Users/azuri/daten/uni/HKU/HASH/PNMain_full_July-03-2020.csv'
fNameHashCSPN = '/Users/azuri/daten/uni/HKU/interns_projects/simba/hash-CSCoords_July-10-2020.csv'
#fNameHashCSPN = '/Users/azuri/daten/uni/HKU/interns_projects/simba/hash-CSCoords.csv'
fNameHashTbCSPN = '/Users/azuri/daten/uni/HKU/interns_projects/simba/hash-tbCSCoords_July-10-2020.csv'
#fNameHashTbCSPN = '/Users/azuri/daten/uni/HKU/interns_projects/simba/hash-tbCSCoords.csv'

sqlFileOut = '/Users/azuri/daten/uni/HKU/interns_projects/simba/fixCSPN.sql'
sqlFileOutZiggy = '/Users/azuri/daten/uni/HKU/interns_projects/simba/fixCSPN_ziggy.sql'

pnMain = csv.DictReader(open(fNameHashPNMain))
cspn = csv.DictReader(open(fNameHashCSPN))
tbCSPN = csv.DictReader(open(fNameHashTbCSPN))

#for row in pnMain:
#    print(row)
#for row in cspn:
#    print(row)
#for row in tbCSPN:
#    print(row)

with open(sqlFileOutZiggy,'w') as f:
    f.write('CREATE TABLE IF NOT EXISTS MainPNData.CSPN_ziggy (\n')
    f.write('idCSPN_ziggy INT UNIQUE AUTO_INCREMENT PRIMARY KEY,\n')
    f.write('idPNMain INT UNIQUE,\n')
    f.write('mapFlag  VARCHAR(1) NOT NULL\n')
    f.write(');\n')
    f.write('USE MainPNData;\n')

def getRaDecForPN(idPNMain):
    pnMainTemp = csv.DictReader(open(fNameHashPNMain))
    for row in pnMainTemp:
        if row['idPNMain'] == idPNMain:
            return [hmsToDeg(row['RAJ2000']), dmsToDeg(row['DECJ2000'])]
    return None

def getidtbCSCoords(idPNMain):
    tbCSPNTemp = csv.DictReader(open(fNameHashTbCSPN))
    for row in tbCSPNTemp:
        if row['idPNMain'] == idPNMain:
            return row['idtbCSCoords']
    return None

def getidPNMain(idtbCSCoords):
#    print('searching for idtbCSCoords = <'+idtbCSCoords+'>')
    tbCSPNTemp = csv.DictReader(open(fNameHashTbCSPN))
    for row in tbCSPNTemp:
#        print("row['idtbCSCoords'] = <"+row['idtbCSCoords']+'> != <'+idtbCSCoords+'>')
        if row['idtbCSCoords'] == idtbCSCoords:
            return row['idPNMain']
    return None

def checkCSPNCoordsAndAddMissingAngles():
#    pnMain = csv.DictReader(open(fNameHashPNMain))
    cspn = csv.DictReader(open(fNameHashCSPN))
#    tbCSPN = csv.DictReader(open(fNameHashTbCSPN))

    iZiggy = 1
    nProblems = 0
    nWrongCoords = 0
    nWrongidPNMain = 0
    with open(sqlFileOut,'w') as f:
        f.write("USE `MainGPN`;\n")
        for line in cspn:
            ra = line['CS_RAJ2000']
            dec = line['CS_DECJ2000']
            dra = line['CS_DRAJ2000']
            ddec = line['CS_DDECJ2000']
            l = line['CS_Glon']
            b = line['CS_Glat']
            idPNMain = line['idPNMain']
            idtbCSCoords = line['idtbCSCoords']
            print('idPNMain = ',idPNMain,', idtbCSCoords = ',idtbCSCoords)
        #    print('idPNMain = ',idPNMain,': line = ',line)
            print('idPNMain = ',idPNMain,': ra = ',ra,', dec = ',dec)
            print('idPNMain = ',idPNMain,': dra = ',dra,', ddec = ',ddec)
            print('idPNMain = ',idPNMain,': l = ',l,', b = ',b)

            isProblem = False
            wrongCoords = False
            if (ra != 'NULL') and (line['userRecord'] != 'chandra0101'):
                draCalc = hmsToDeg(ra)
                ddecCalc = dmsToDeg(dec)
                lCalc, bCalc = raDecToLonLat(draCalc,ddecCalc)
                if dec == 'NULL':
                    print('Problem with idPNMain = ',idPNMain,': ra is not empty but dec is => aborting')
                    isProblem = True
                    STOP
                if dra == 'NULL':
                    isProblem = True
                else:
                    if abs((angularDistancePyAsl(float(dra), float(ddec), draCalc, ddecCalc)*3600.) - angularDistance(float(dra), float(ddec), draCalc, ddecCalc)) > 0.1:
                        print('separation = ',abs((angularDistancePyAsl(float(dra), float(ddec), draCalc, ddecCalc)*3600.) - angularDistance(float(dra), float(ddec), draCalc, ddecCalc)))
                        print('problem: angular distances differ by more than 0.1 arcsec')
                        wrongCoords = True
                    if angularDistance(float(dra), float(ddec), draCalc, ddecCalc) > 0.1:
                        print('angular distance = ',angularDistance(float(dra), float(ddec), draCalc, ddecCalc))
                        print('problem: angular distances between stated and calculated dRA and dDEC differ by more than 0.1 arcsec')
                        wrongCoords = True
                if ddec == 'NULL':
                    isProblem = True
                if l == 'NULL':
                    isProblem = True
                else:
                    raFromlb, decFromlb = lonLatToRaDec(float(l), float(b))
                    raFromlbCalc, decFromlbCalc = lonLatToRaDec(lCalc, bCalc)
                    if angularDistance(raFromlb, decFromlb, raFromlbCalc, decFromlbCalc) > 0.1:
                        print('lCalc = ',lCalc,', bCalc = ',bCalc)
                        print('angular distance = ',angularDistance(raFromlb, decFromlb, raFromlbCalc, decFromlbCalc))
                        print('problem: angular distances between stated and calculated l and b differ by more than 0.1 arcsec')
                        wrongCoords = True
                if b == 'NULL':
                    isProblem = True
            else:
                #remove from both tables
    #            if getidtbCSCoords(idPNMain) != line['idtbCSCoords']:
    #                print("PROBLEM: getidtbCSCoords(idPNMain)=",getidtbCSCoords(idPNMain)," != line['idtbCSCoords']=",line['idtbCSCoords'])
                f.write('DELETE FROM `tbCSCoords` WHERE `idtbCSCoords` = %d;\n' % int(idtbCSCoords))
                print('deleted record from tbCSCoords')
                if line['InUse'] == '1':
                    f.write('DELETE FROM `PNMain_tbCSCoords` WHERE `idtbCSCoords` = %d;\n' % int(idtbCSCoords))
                    print('deleted record from PNMain_tbCSCoords')
                cspnTemp = csv.DictReader(open(fNameHashCSPN))
                setInUse = 0
                for lineTemp in cspnTemp:
                    if lineTemp['idPNMain'] == idPNMain:
                        if lineTemp['userRecord'] != 'chandra0101':
                            if lineTemp['InUse'] == '1':
                                print('Found another entry for idPNMain=',idPNMain,' which was set InUse')
                                setInUse = 1
                if setInUse == 0:
                    cspnTemp = csv.DictReader(open(fNameHashCSPN))
                    for lineTemp in cspnTemp:
                        if lineTemp['idPNMain'] == idPNMain:
                            if (lineTemp['userRecord'] != 'chandra0101') and (setInUse == 0):
                                print('lineTemp[userRecord] = ',lineTemp['userRecord'])
                                setInUse = 1
                                f.write('UPDATE `tbCSCoords` SET `InUse` = 1 WHERE `idtbCSCoords` = %d;\n' % (int(lineTemp['idtbCSCoords'])))
                                f.write("INSERT INTO `PNMain_tbCSCoords`(`idPNMain`,`idtbCSCoords`)")
                                f.write(" VALUES (%d,%d);\n" % (int(idPNMain),
                                                                int(lineTemp['idtbCSCoords'])))
                                print('updated tbCSCoords and added entry to PNMain_tbCSCoords')
                                #STOP
                            elif (lineTemp['userRecord'] != 'chandra0101') and (setInUse == 1):
                                print('Problem: found another entry with InUse == 1')
                                STOP


    #        idtbCSCoords = getidtbCSCoords(idPNMain)
    #        print('idtbCSCoords = ',idtbCSCoords)
            addedTotbCSCoords = False
            if isProblem:
                if getidtbCSCoords(idPNMain) != idtbCSCoords:
                    print("PROBLEM: getidtbCSCoords(idPNMain)=",getidtbCSCoords(idPNMain)," != line['idtbCSCoords']=",line['idtbCSCoords'])
    #                idtbCSCoords = getidtbCSCoords(idPNMain)
                if getidtbCSCoords(idPNMain) is None:
                    if line['InUse'] == '1':
                        f.write("INSERT INTO `PNMain_tbCSCoords`(`idPNMain`,`idtbCSCoords`)")
                        f.write(" VALUES (%d,%d);\n" % (int(idPNMain),
                                                        int(idtbCSCoords)))
                        print('added entry to PNMain_tbCSCoords')
                        addedTotbCSCoords = True
    #                STOP

                nProblems += 1
                f.write('UPDATE `tbCSCoords` SET `CS_DRAJ2000` = %.5f, `CS_DDECJ2000` = %.5f, `CS_Glon` = %.5f, `CS_Glat` = %.5f WHERE `idtbCSCoords` = %d;\n' % (draCalc, ddecCalc, lCalc, bCalc, int(idtbCSCoords)))
                print('updated tbCSCoords')

            if wrongCoords:
                if line['InUse'] == '1':
                    if getidtbCSCoords(idPNMain) != idtbCSCoords:
                        print("PROBLEM: getidtbCSCoords(idPNMain)=",getidtbCSCoords(idPNMain)," != line['idtbCSCoords']=",line['idtbCSCoords'])
                        STOP
    #                idtbCSCoords = getidtbCSCoords(idPNMain)
                nWrongCoords += 1

                #check RA and DEC with the PN
                draPN, ddecPN = getRaDecForPN(idPNMain)
                if angularDistance(draPN, ddecPN, draCalc, ddecCalc) > 20.:
                    print('draPN = ',draPN,', ddecPN = ',ddecPN)
                    print('angular distance = ',angularDistance(draPN, ddecPN, draCalc, ddecCalc))
                    print('problem: angular distances between stated and calculated l and b differ by more than 20. arcsec')
                    #STOP
                f.write('UPDATE `tbCSCoords` SET `CS_DRAJ2000` = %.5f, `CS_DDECJ2000` = %.5f, `CS_Glon` = %.5f, `CS_Glat` = %.5f WHERE `idtbCSCoords` = %d;\n' % (draCalc, ddecCalc, lCalc, bCalc, int(idtbCSCoords)))
            if line['InUse'] == '1':
                if getidPNMain(line['idtbCSCoords']) is None:
                    print('idtbCSCoords = '+idtbCSCoords+' not found in PNMain_tbCSCoords')
                    if not addedTotbCSCoords:
                        f.write("INSERT INTO `PNMain_tbCSCoords`(`idPNMain`,`idtbCSCoords`)")
                        f.write(" VALUES (%d,%d);\n" % (int(idPNMain),
                                                        int(line['idtbCSCoords'])))
                        print('added entry to PNMain_tbCSCoords')
                        addedTotbCSCoords = True
                        STOP
                elif getidPNMain(line['idtbCSCoords']) != idPNMain:
                    print('PROBLEM: wrong idPNMain(=',idPNMain,') in PNMain_tbCSCoords(getidPNMain(',line['idtbCSCoords'],') = ',getidPNMain(line['idtbCSCoords']),')')
                    nWrongidPNMain += 1
                    f.write('DELETE FROM `PNMain_tbCSCoords` WHERE `idtbCSCoords` = %d;\n' % int(line['idtbCSCoords']))
                    print('deleted entry from PNMain_tbCSCoords')
                    f.write("INSERT INTO `PNMain_tbCSCoords`(`idPNMain`,`idtbCSCoords`)")
                    f.write(" VALUES (%d,%d);\n" % (int(idPNMain),
                                                    int(idtbCSCoords)))
                    print('added entry to PNMain_tbCSCoords')
    #                STOP
            if line['userRecord'] == 'ziggy':
                with open(sqlFileOutZiggy,'a') as fZiggy:
                    fZiggy.write("INSERT INTO `CSPN_ziggy`(`idCSPN_ziggy`,`idPNMain`,`mapFlag`)")
                    fZiggy.write(" VALUES (%d,%d,'%s');\n" % (iZiggy,
                                                            int(idPNMain),
                                                            'y'))
                iZiggy+=1


    print('found ',nProblems,' problematic CSPN and ',nWrongCoords,' CSPN with problematic coordinates, ',nWrongidPNMain,' wrong idPNMain in PNMain_tcCSCoords')

def fixAccuracy():
    with open(sqlFileOut[:sqlFileOut.rfind('.')]+'_accuracy.sql','w') as f:
        f.write("USE `MainGPN`;\n")
        cspn = csv.DictReader(open(fNameHashCSPN))
        for row in cspn:
            problem = False
            ra = row['CS_RAJ2000']
            dec = row['CS_DECJ2000']
            if '\t' in ra:
                print('found a tab in ra')
                ra = ra.strip()
                problem = True
            if '\t' in dec:
                print('found a tab in dec')
                problem = True
                dec = dec.strip()

            rah, ram, ras = [i for i in ra.split(':')]
            accuracyRA = len(ras[ras.rfind('.')+1:])
            decd, decm, decs = [i for i in dec.split(':')]
            accuracyDEC = len(decs[decs.rfind('.')+1:])
            print('ra = '+ra+': accuracyRA = ',accuracyRA)
            print('dec = '+dec+': accuracyDEC = ',accuracyDEC)
            if accuracyRA > 2:
                problem = True
                ra = rah+':'+ram+':'+'%.2f' % float(ras)
                print('new ra = ',ra)
            if accuracyDEC > 1:
                problem = True
                dec = decd+':'+decm+':'+'%.1f' % float(decs)
                print('new dec = ',dec)
            if problem:
                draCalc = hmsToDeg(ra)
                ddecCalc = dmsToDeg(dec)
                lCalc, bCalc = raDecToLonLat(draCalc,ddecCalc)
                f.write("UPDATE `tbCSCoords` SET `CS_RAJ2000` = '%s', `CS_DECJ2000` = '%s', `CS_DRAJ2000` = %.5f, `CS_DDECJ2000` = %.5f, `CS_Glon` = %.5f, `CS_Glat` = %.5f WHERE `idtbCSCoords` = %d;\n" % (ra, dec, draCalc, ddecCalc, lCalc, bCalc, int(row['idtbCSCoords'])))

def addSimbasComments():
    simbaFName = '/Users/azuri/daten/uni/HKU/interns_projects/simba/All_list_v2.csv'
    simba = csv.DictReader(open(simbaFName))

    for row in simba:
        if row['Comments detailed'] != '':
            idPNMain = row['Hash id']
            cspn = csv.DictReader(open(fNameHashCSPN))
            for line in cspn:
                if (line['idPNMain'] == idPNMain) and (row[''])




if __name__ == '__main__':
    checkCSPNCoordsAndAddMissingAngles()
    fixAccuracy()
