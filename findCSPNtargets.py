from astroplan import Observer
#from astroplan import download_IERS_A
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon
import numpy as np

from astropy.utils.iers import conf
conf.auto_max_age = None

from findTargets import findTargetsVisibleAt,getDistanceToMoon,readCSV,writeCSV
from plot_obs_planning import plot_target

allPossibleTargets = '/Users/azuri/daten/uni/HKU/HASH/CSPN_GTC_targets/allTruePNe.csv'
targetListOut = allPossibleTargets[:-4]+'_good.csv'

date = '2020-12-01'
midnightTime = Time(date+' 0:00:00')

def findCSPNe(csvTargets):
    csvOut = []
    for line in csvTargets:
        if (line['CS_DRAJ2000'] != '') and (line['CS_DDECJ2000'] != ''):
            csvOut.append(line)
    print('found ',len(csvOut),' CSPN')
    return csvOut

def findGoodTargets(targetsIn,location,raString='CS_DRAJ2000',decString='CS_DDECJ2000'):
    goodTargets = []
    for line in targetsIn:
        ra = float(line[raString])
        dec = float(line[decString])
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
        print(line['Name']+': RA = '+line[raString]+', DEC = '+line[decString]+': max(altitude) = ',maxAltitude)
        whereGTminAlt = np.where(np.array(altitude) > minimumAltitude)[0]
#            print('whereGTminAlt = ',whereGTminAlt)
        nGoodHours = len(whereGTminAlt) / 60.
        print('object can be observed for ',nGoodHours,' hours')
        if nGoodHours > 1.5:
            goodTargets.append(line)
    return goodTargets

def writeTargetList(fNameIn, fNameOut, raString='CS_RAJ2000', decString='CS_DECJ2000'):
    all = readCSV(fNameIn)
    print('writeTargetList: read ',len(all),' lines')
    with open(fNameOut,'w') as f:
        for lineIn in all:
            lineOut = lineIn['Name'].replace(' ','_').replace('"','')+' '+lineIn[raString]+' '+lineIn[decString]+' J2000\n'
            f.write(lineOut)

def createTargetLists():
    allSiteNames = EarthLocation.get_site_names()
    observatoryName = "Roque de los Muchachos"

    gtcObs = Observer.at_site(observatoryName)#, timezone='Eastern Standard Time')
    print('dir(gtcObs) = ',dir(gtcObs))
    print('gtcObs.location = ',gtcObs.location)
#    gtc = EarthLocation(lat=dmsToDeg()*u.deg, lon=149.0685*u.deg, height=1165*u.m)
    gtc_utcoffset = 1*u.hour # Australian Eastern Standard Time

    minimumAltitude = 30.

    midnight = midnightTime - gtc_utcoffset

    #print('astronomical twilight as Observatory: %s - %s' % (obs.twilight_evening_astronomical(midnight), obs.twilight_morning_astronomical(midnight)))
    observationStartTime = gtcObs.twilight_evening_astronomical(midnight)#Time('2019-9-5 19:10:00') - utcoffset
    observationEndTime = gtcObs.twilight_morning_astronomical(midnight)#Time('2019-9-6 4:55:00') - utcoffset
    observationStartTime.format = 'iso'
    observationEndTime.format = 'iso'

    catLines = readCSV(allPossibleTargets)
    targetsWithCS = findCSPNe(catLines)

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

    goodTargets = findGoodTargets(targetsWithCS,gtcObs.location,raString='CS_DRAJ2000',decString='CS_DDECJ2000')
    writeCSV(goodTargets,targetListOut,'DRAJ2000')

    allGoodTargets = []
    allGoodTargetIDs = []
    for time in fullHours:
        visibleAt = findTargetsVisibleAt(goodTargets, time, gtcObs.location, gtc_utcoffset, minimumAltitude)
        localTime = time + gtc_utcoffset
        writeCSV(visibleAt, targetListOut[:-4]+'_visible_at_'+date+'_'+localTime.strftime("%H-%M")+'.csv', 'DRAJ2000')

#        print('visibleAt[0] = ',visibleAt[0])
        ra = float(visibleAt[0]['CS_DRAJ2000'])
        dec = float(visibleAt[0]['CS_DDECJ2000'])
        targetCoord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
#        print('ra, dec = ',ra,dec)
#        plot_target(targetCoord, gtcObs.location, observatoryName, gtc_utcoffset, date, False)
        print('ra, dec = ',ra,dec)

        for line in visibleAt:
            if line['idPNMain'] not in allGoodTargetIDs:
                allGoodTargetIDs.append(line['idPNMain'])
                allGoodTargets.append(line)
    writeCSV(allGoodTargets,targetListOut[:-4]+'_visible_at_'+date+'.csv')

if __name__ == '__main__':
    #createTargetLists()
    visibilityLists = '/Users/azuri/daten/uni/HKU/HASH/CSPN_GTC_targets/allTruePNe_withCS_visible.list'
    with open(visibilityLists,'r') as vf:
        lists = [f.rstrip() for f in vf.readlines()]
    allCSPN = []
    allCSPNIDs = []
    for list in lists:
        csv = readCSV(list)
        for line in csv:
            if line['idPNMain'] not in allCSPNIDs:
                allCSPNIDs.append(line['idPNMain'])
                allCSPN.append(line)
    writeCSV(allCSPN,visibilityLists[:-4]+'csv')
