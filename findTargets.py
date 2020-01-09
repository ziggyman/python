from astroplan import Observer
#from astroplan import download_IERS_A
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon
import numpy as np

from astropy.utils.iers import conf
conf.auto_max_age = None

midnightTime = Time('2020-8-1 0:00:00')

# string = xx:yy:zz.zzz
def hmsToDeg(string):
    h, m, s = [float(i) for i in string.split(':')]
    return (s / 240.) + (m / 4.) + (h * 15.)

# string = xx:yy:zz.zzz
def dmsToDeg(string):
    d, m, s = [float(i) for i in string.split(':')]
    if d < 0.:
        d = 0. - d
        return 0. - (s / 3600. + m / 60. + d)
    return s / 3600. + m / 60. + d

def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

def readCSV(fName):
    lines = readFileToArr(fName)
    print('read ',len(lines),' lines from ',fName)
    lines = [line.split(',') for line in lines]
    for iLine in np.arange(0,len(lines),1):
        lines[iLine] = [x.strip(' ') for x in lines[iLine]]
    catalogueKeys = []
    for key in lines[0]:
        catalogueKeys.append(key)
#    print('catalogueKeys = ',catalogueKeys)
    catLines = []
    for iLine in np.arange(1,len(lines),1):
        catLine = {catalogueKeys[0]:lines[iLine][0]}
        for iKey in np.arange(1,len(catalogueKeys),1):
            catLine.update({catalogueKeys[iKey]:lines[iLine][iKey]})
#        print('catLine = ',catLine)
        catLines.append(catLine)
#    print('lines = ',lines)
#    print('lines[0] = ',lines[0])
#    print('lines[1] = ',lines[1])
#    print('catLines = ',catLines)
    return catLines

def writeCSVHeader(dataLine, fileOut):
    keys = list(dataLine.keys())
    line = keys[0]
    for key in keys[1:]:
        line += ', '+key
    line += '\n'
    fileOut.write(line)

def writeCSVLine(dataLine, keys, fileOut):
    line = dataLine[keys[0]]
    for key in keys[1:]:
        line += ', '+dataLine[key]
    line += '\n'
    fileOut.write(line)

def writeCSV(data, fName, sortKey = None):
    with open(fName,'w') as f:
        if len(data) > 0:
            dataSorted = data
            keys = list(data[0].keys())
            if sortKey is not None:
                dataSorted = sorted(data, key=lambda k: k[sortKey])
            writeCSVHeader(data[0],f)
            for dataLine in dataSorted:
                writeCSVLine(dataLine, keys, f)
        else:
            f.write('sorry no targets found')
    print('wrote ',fName,' with ',len(data),' data lines')

def removeAlreadyObserved(data, fName='/Users/azuri/daten/uni/HKU/observing/already_observed.list'):
    alreadyObserved = readFileToArr(fName)
    dataOut = []
    nRemoved=0
    for line in data:
        targetName = line['Name'].strip('"')
        found = False
        for name in alreadyObserved:
            if name == targetName:
                nRemoved += 1
                found = True
#                print('found name <'+name+'> in targetName = <'+targetName+'>')
        if not found:
            dataOut.append(line)
    print('removed ',nRemoved,' targets which have already been observed but the spectra have not yet been ingested into the HASH database')
    return dataOut

def moveTargetsStartingWithToNewList(fNameIn, namePrefix, fNameGoodOut, fNameRejectedOut, append=True):
    all = readCSV(fNameIn)
    if len(all) > 0:
        keys = list(all[0].keys())
        nRemoved=0
        with open(fNameGoodOut, 'w') as fGood:
            writeCSVHeader(all[0], fGood)
            appen = 'w'
            if append:
                appen = 'a'
            with open(fNameRejectedOut, appen) as fReject:
                if appen == 'w':
                    # write header
                    writeCSVHeader(all[0], fReject)
                for target in all:
                    if target['Name'].strip('"')[:len(namePrefix)] == namePrefix:
                        nRemoved += 1
                        writeCSVLine(target, keys, fReject)
                    else:
                        writeCSVLine(target, keys, fGood)
        print('removed ',nRemoved,' targets which start with <'+namePrefix+'>')

def writeSAAOTargetList(fNameIn, fNameOut):
    all = readCSV(fNameIn)
    print('writeSAAOTargetList: read ',len(all),' lines')
    with open(fNameOut,'w') as f:
        for lineIn in all:
            lineOut = lineIn['Name'].replace(' ','_').replace('"','')+' '+lineIn['RAJ2000']+' '+lineIn['DECJ2000']+' J2000\n'
            f.write(lineOut)

def findTargetsVisibleAt(targets, timeUTC, location, utcoffset, minimumAltitude):
    targetsOut = []
    for line in targets:
        ra = float(line['DRAJ2000'])
        dec = float(line['DDECJ2000'])
        raHMS = line['RAJ2000']
        decDMS = line['DECJ2000']
#        print(line['Name']+': ra = ',raHMS,', dec = ',decDMS)
        targetCoord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        altaz = targetCoord.transform_to(AltAz(obstime=timeUTC,location=location))
        altitude = float('{0.alt:.2}'.format(altaz).split(' ')[0])
        localTime = timeUTC+utcoffset
#        print('altitude of '+line['Name']+' at '+localTime.strftime("%H:%M")+' is '+str(altitude))
        if altitude > minimumAltitude:
            line.update({'Altitude': '%d' % int(altitude)})
            targetsOut.append(line)
#            print('object is visible')
    return targetsOut

def takeOnlyWhatIsInBothInputLists(fileA, fileB, fileNameOut):

    targetsA = readCSV(fileA)
    targetsB = readCSV(fileB)
    print('len('+fileA+') = ',len(targetsA),', len('+fileB+') = ',len(targetsB))

    goodTargets = []

    for lineA in targetsA:
        for lineB in targetsB:
            if lineA['idPNMain'] == lineB['idPNMain']:
                goodTargets.append(lineA)

    writeCSV(goodTargets, fileNameOut, 'DRAJ2000')
    print('found ',len(goodTargets),' in both files')
    return goodTargets

def getDistanceToMoon(location, starCoord, time):
#    print('dir(ICRS) = ',dir(ICRS))
    moon = get_moon(time, location=location)
#    print('moon = ',moon)
#    moonICRS = moon.transform_to(ICRS)
#    print('starCoord = ',starCoord)
    frameNight = AltAz(obstime=time,#midnight+delta_midnight,
                       location=location)
#    print('type(frameNight) = ',type(frameNight))
#    print('frameNight = ',frameNight)
    targetaltazsNight = starCoord.transform_to(frameNight)
    moonaltazsNight = moon.transform_to(frameNight)
#    print('targetaltazsNight = ',targetaltazsNight)
#    print('moonaltazsNight = ',moonaltazsNight)
    separation = targetaltazsNight.separation(moonaltazsNight)
#    print('separation = ',separation)
#    print('dir(separation) = ',dir(separation))
#    print('separation = ',separation.degree,' degrees')
#    print('dir(separation.degree) = ',dir(separation.degree))
#    STOP
    return separation.degree

def main():
    for site in ['SAAO']:#,'SSO']:
        location = None
        for noDiameterPN in [False]:#, True]:
            for priority in [True, False]:#

                if priority:
                    allPossibleTargets = '/Users/azuri/daten/uni/HKU/observing/all_targets_truePN_MPA.csv'
                else:
                    allPossibleTargets = '/Users/azuri/daten/uni/HKU/observing/all_targets_with_catalogue_without_DFrew_noTrue.csv'

                ssoObs = Observer.at_site("Siding Spring Observatory")#, timezone='Eastern Standard Time')
                sso = EarthLocation(lat=-31.2749*u.deg, lon=149.0685*u.deg, height=1165*u.m)
                sso_utcoffset = 10*u.hour # Australian Eastern Standard Time
                ssoOutFileName = allPossibleTargets[0:allPossibleTargets.rfind('/')]+'/targets_SSO'
                ssoOutFileName += '_'+str(midnightTime.datetime.date())
                if priority:
                    ssoOutFileName += '_priority'
                ssoOutFileName += '_new.csv'
                ssoMinimumAltitude = 30.
                ssoMinimumMajorDiameter = 0.
                ssoMaximumMajorDiameter = 50.

                saaoObs = Observer.at_site("Southern African Large Telescope")#, timezone='Eastern Daylight Time')
                saao = EarthLocation(lat=-32.3783*u.deg, lon=20.8105*u.deg, height=1750*u.m)
                saao_utcoffset = 2*u.hour  # Eastern Daylight Time
                saaoOutFileName = allPossibleTargets[0:allPossibleTargets.rfind('/')]+'/targets_SAAO'
                saaoOutFileName += '_'+str(midnightTime.datetime.date())
                if priority:
                    saaoOutFileName += '_priority'
                saaoOutFileName += '.csv'
                saaoMinimumAltitude = 30.
                saaoMinimumMajorDiameter = 0.
                saaoMaximumMajorDiameter = 30.

                if site == 'SSO':
                    obs = ssoObs
                    location = sso
                    utcoffset = sso_utcoffset
                    minimumAltitude = ssoMinimumAltitude
                    minimumMajorDiameter = ssoMinimumMajorDiameter
                    maximumMajorDiameter = ssoMaximumMajorDiameter
                    outFileName = ssoOutFileName
                elif site == 'SAAO':
                    obs = saaoObs
                    location = saao
                    utcoffset = saao_utcoffset
                    minimumAltitude = saaoMinimumAltitude
                    minimumMajorDiameter = saaoMinimumMajorDiameter
                    maximumMajorDiameter = saaoMaximumMajorDiameter
                    outFileName = saaoOutFileName
                else:
                    print('could not identify site "'+site+'"')
                    STOP

                if noDiameterPN:
                    outFileName = outFileName[0:outFileName.rfind('.')]+'_noDiamGiven.csv'

                midnight = midnightTime - utcoffset

                #print('astronomical twilight as Observatory: %s - %s' % (obs.twilight_evening_astronomical(midnight), obs.twilight_morning_astronomical(midnight)))
                observationStartTime = obs.twilight_evening_astronomical(midnight)#Time('2019-9-5 19:10:00') - utcoffset
                observationEndTime = obs.twilight_morning_astronomical(midnight)#Time('2019-9-6 4:55:00') - utcoffset
                observationStartTime.format = 'iso'
                observationEndTime.format = 'iso'
                print('observationStartTime in local time = ',observationStartTime+utcoffset)
                print('observationEndTime in local time = ',observationEndTime+utcoffset)
                print('observationStartTime in UT = ',observationStartTime)
                print('observationEndTime in UT = ',observationEndTime)

                #observationStartTime = Time('2019-9-5 19:10:00') - utcoffset
                #observationEndTime = Time('2019-9-6 4:55:00') - utcoffset

                catLines = readCSV(allPossibleTargets)

                time = observationStartTime
                times = [time]
                fullHours = []
                while time < observationEndTime:
                    time += 1*u.minute
#                    print('time = ',time)
                    times.append(time)
#                    print('time = ',time)
#                    print('type(time) = ',type(time))
#                    print('dir(time) = ',dir(time))
#                    print('type(time.to_datetime()) = ',type(time.to_datetime()))
#                    print('dir(time.to_datetime()) = ',dir(time.to_datetime()))
                    if time.to_datetime().strftime("%M") == '00':
#                        print('full hour found at ',time)
                        fullHours.append(time)

                goodTargets = []
                for line in catLines:
                    diamStr = line['MajDiam']
                    print('diamStr = ',diamStr)
                    if (not noDiameterPN) and (diamStr != ''):
                        diam = float(diamStr)
                        if (diam >= minimumMajorDiameter) and (diam <= maximumMajorDiameter):
                            ra = float(line['DRAJ2000'])
                            dec = float(line['DDECJ2000'])
#                            print('ra = ',ra,', dec = ',dec)
                            targetCoord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
                            moonDistance = getDistanceToMoon(location, targetCoord, midnight)
                            print('moonDistance = ',moonDistance)
                            nGoodHours = 0
                            altaz = targetCoord.transform_to(AltAz(obstime=np.array(times),location=location))
                            altitude = [float('{0.alt:.2}'.format(alt).split(' ')[0]) for alt in altaz]
                            #print(line['Name']+': RA = '+line['RAJ2000']+', DEC = '+line['DECJ2000']+': altitude = ',altitude)
                            maxAltitude = np.max(altitude)
                            if line['Name'] == '':
                                line.update({'Name': 'PNG'+line['PNG']})
                            line.update({'maxAlt':'%.1f'%maxAltitude})
                            line.update({'moonDist':'%.0f'%moonDistance})
                #            print('altitude = ',altitude)
                            print(line['Name']+': RA = '+line['RAJ2000']+', DEC = '+line['DECJ2000']+': max(altitude) = ',maxAltitude)
                            whereGTminAlt = np.where(np.array(altitude) > minimumAltitude)[0]
                #            print('whereGTminAlt = ',whereGTminAlt)
                            nGoodHours = len(whereGTminAlt) / 60.
                            print('object can be observed for ',nGoodHours,' hours')
                            if nGoodHours > 1.5:
                                goodTargets.append(line)
#                                print('target is possible to observe')
                    elif noDiameterPN and (diamStr == ''):
                        print('no diameter given')
                        if line['MajDiam'] != '':
                            STOP
                        ra = float(line['DRAJ2000'])
                        dec = float(line['DDECJ2000'])
#                        print('ra = ',ra,', dec = ',dec)
                        targetCoord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
                        moonDistance = getDistanceToMoon(location, targetCoord, midnight)
                        nGoodHours = 0
                        altaz = targetCoord.transform_to(AltAz(obstime=np.array(times),location=location))
                        altitude = [float('{0.alt:.2}'.format(alt).split(' ')[0]) for alt in altaz]
                        if line['Name'] == '':
                            line.update({'Name': 'PNG'+line['PNG']})
                        maxAltitude = np.max(altitude)
                        line.update({'maxAlt':'%.1f'%maxAltitude})
                        line.update({'moonDist':'%.0f'%moonDistance})
                        print(line['Name']+': RA = '+line['RAJ2000']+', DEC = '+line['DECJ2000']+': altitude = ',maxAltitude)
            #            print('altitude = ',altitude)
                        print('max(altitude) = ',maxAltitude)
                        whereGTminAlt = np.where(np.array(altitude) > minimumAltitude)[0]
            #            print('whereGTminAlt = ',whereGTminAlt)
                        nGoodHours = len(whereGTminAlt) / 60.
                        print('object can be observed for ',nGoodHours,' hours')
                        if nGoodHours > 1.5:
                            goodTargets.append(line)
#                            print('target is possible to observe')

                #            if float(line['DRAJ2000']) < 199:
                #                STOP
                goodLength = len(goodTargets)
                print('found ',len(goodTargets),' good targets')

                goodTargets = removeAlreadyObserved(goodTargets)
                print('removed ',goodLength - len(goodTargets),' targets already observed')
                writeCSV(goodTargets,outFileName,'DRAJ2000')

                goodFileName = outFileName[0:outFileName.rfind('.')]+'_good.csv'

                moveTargetsStartingWithToNewList(outFileName,
                                                 'DeGaPe',
                                                 goodFileName,
                                                 outFileName[0:outFileName.rfind('.')]+'_DeGaPe.csv',
                                                 append=False)

                moveTargetsStartingWithToNewList(goodFileName,
                                                 'MGE',
                                                 goodFileName,
                                                 outFileName[0:outFileName.rfind('.')]+'_MGE.csv',
                                                 append=False)

                if (not noDiameterPN) and (obs == ssoObs):
                    fileA = goodFileName
                    fileB = outFileName[:outFileName.rfind('_')]+'_good.csv'

                    goodFileName = fileB[:fileB.rfind('.')]+'_new.csv'
                    print('fileA = <'+fileA+'>, fileB = <'+fileB+'>, goodFileName = <'+goodFileName+'>')
                    takeOnlyWhatIsInBothInputLists(fileA, fileB, goodFileName)


                writeSAAOTargetList(goodFileName,goodFileName[0:goodFileName.rfind('.')]+'.dat')

                print('fullHours = ',fullHours)
                for fName in [goodFileName, outFileName[0:outFileName.rfind('.')]+'_DeGaPe.csv', outFileName[0:outFileName.rfind('.')]+'_MGE.csv']:
                    targets = readCSV(fName)
                    for time in fullHours:
                        visibleAt = findTargetsVisibleAt(targets, time, location, utcoffset, minimumAltitude)
                        localTime = time + utcoffset
                        writeCSV(visibleAt, fName[:fName.rfind('.')]+'_visible_at_'+localTime.strftime("%H-%M")+'.csv', 'DRAJ2000')

if __name__ == "__main__":
    main()
