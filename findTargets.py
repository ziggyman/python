from astroplan import Observer
from astroplan import download_IERS_A
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
from time import clock

#sso = Observer.at_site("Siding Spring Observatory")#, timezone='Eastern Standard Time')
#saao = Observer.at_site("Southern African Large Telescope")#, timezone='Eastern Daylight Time')

sso = EarthLocation(lat=-31.2749*u.deg, lon=149.0685*u.deg, height=1165*u.m)
sso_utcoffset = 10*u.hour # Australian Eastern Standard Time
saao = EarthLocation(lat=-32.3783*u.deg, lon=20.8105*u.deg, height=1750*u.m)
saao_utcoffset = 2*u.hour  # Eastern Daylight Time

location = sso
utcoffset = sso_utcoffset

observationStartTime = Time('2019-9-5 19:10:00') - utcoffset
observationEndTime = Time('2019-9-6 4:55:00') - utcoffset

minimumAltitude = 30.

minimumMajorDiameter = 10.
maximumMajorDiameter = 100.

allPossibleTargets = '/Users/azuri/daten/uni/HKU/observing/all_targets.csv'
outFileName = allPossibleTargets[0:allPossibleTargets.rfind('/')]+'/targets_SSO.csv'
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
    lines = [line.split(',') for line in lines]
    catalogueKeys = []
    for key in lines[0]:
        catalogueKeys.append(key)
    print('catalogueKeys = ',catalogueKeys)
    catLines = []
    for iLine in np.arange(1,len(lines),1):
        catLine = {catalogueKeys[0]:lines[iLine][0]}
        for iKey in np.arange(1,len(catalogueKeys),1):
            catLine.update({catalogueKeys[iKey]:lines[iLine][iKey]})
        print('catLine = ',catLine)
        catLines.append(catLine)
    print('lines = ',lines)
    print('lines[0] = ',lines[0])
    print('lines[1] = ',lines[1])
    print('catLines = ',catLines)
    return catLines

def writeCSV(data, fName, sortKey = None):
    dataSorted = data
    if sortKey is not None:
        dataSorted = sorted(data, key=lambda k: k[sortKey])
    print('dir(data[0]) = ',dir(data[0]))
    print('data[0].keys() = ',data[0].keys())
    keys = list(data[0].keys())
    with open(fName,'w') as f:
        line = keys[0]
        for key in keys[1:]:
            line += ','+key
        line += '\n'
        f.write(line)

        for dataLine in dataSorted:
            line = dataLine[keys[0]]
            for key in keys[1:]:
                line += ','+dataLine[key]
            line += '\n'
            f.write(line)

catLines = readCSV(allPossibleTargets)

time = observationStartTime
times = [time]
while time < observationEndTime:
    time += 1*u.hour
    print('time = ',time)
    times.append(time)


goodTargets = []
for line in catLines:
    diamStr = line['MajDiam']
    if diamStr != '':
        diam = float(diamStr)
        if (diam >= minimumMajorDiameter) and (diam <= maximumMajorDiameter):
            ra = float(line['DRAJ2000'])
            dec = float(line['DDECJ2000'])
            print('ra = ',ra,', dec = ',dec)
            targetCoord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
            nGoodHours = 0
#            for time in times:
            if True:
                print(clock())
                altaz = targetCoord.transform_to(AltAz(obstime=np.array(times),location=location))
#                altaz = targetCoord.transform_to(AltAz(obstime=time,location=location))
                print(clock())
                altitude = [float('{0.alt:.2}'.format(alt).split(' ')[0]) for alt in altaz]
#                print("target's Altitude = {0.alt:.2}".format(altaz))
                print('altitude = ',altitude)
                whereGTminAlt = np.where(np.array(altitude) > minimumAltitude)[0]
                print('whereGTminAlt = ',whereGTminAlt)
                nGoodHours = len(whereGTminAlt)
                print('object can be observed for ',nGoodHours,' hours')
#                    nGoodHours += 1
            if nGoodHours > 0:
                goodTargets.append(line)
                print('target is possible to observe')
print('found ',len(goodTargets),' good targets')

writeCSV(goodTargets,outFileName,'DRAJ2000')