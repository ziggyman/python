import os
from distutils.dir_util import copy_tree
import numpy as np
from PIL import Image
import shutil

from astroplan import Observer, FixedTarget
from astropy.time import Time,TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon
import astropy.units as u
from astroplan import AltitudeConstraint, AirmassConstraint,AtNightConstraint
from astroplan import observability_table

from plot_obs_planning import plot_target
import csvFree,csvData
from myUtils import angularDistancePyAsl,hmsToDeg,dmsToDeg


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


def getDarkTimesUT(midnight):
    time = midnight - 8*u.hour
    times = [time]
    while time < midnight + 8*u.hour:
        time += 1*u.minute
        times.append(time)
    utTimes = Time([timeX.value for timeX in times])
    frames = AltAz(obstime=utTimes, location=observatoryLocation)
    sunaltazs = get_sun(utTimes).transform_to(frames)
    utTimes = utTimes[np.where(sunaltazs.alt < -0*u.deg)]
    return utTimes

allSiteNames = EarthLocation.get_site_names()
for siteName in allSiteNames:
    print(siteName)
    obsTest = Observer.at_site(siteName)
    print('obsTest = ',obsTest)
    if 'Alto' in siteName:
        STOP
#allSiteNames = Observer.get_site_names()
#for siteName in allSiteNames:
#    print(siteName)
#    if 'Alto' in siteName:
#        STOP
#STOP
latDeg = dmsToDeg('37:13:46')#37.22361
lonDeg = 360.-dmsToDeg('2:32:30')#2.54611

elev = 2168*u.m

print('latDeg = ',latDeg,', lonDeg = ',lonDeg)
observatoryLocation = EarthLocation(lat=latDeg*u.deg, lon=lonDeg*u.deg, height=elev)
utcoffset = 0*u.hour

observatoryName = "Calar Alto"
observer = Observer(longitude=lonDeg, latitude=latDeg, elevation=elev, name=observatoryName)


pa30_dra = 13.29667
pa30_ddec = 67.50067
targetCoord = SkyCoord(ra=pa30_dra*u.deg, dec=pa30_ddec*u.deg, frame='icrs')
target = FixedTarget(coord=targetCoord, name='Pa 30')

time_range = Time(["2023-07-01 00:00", "2023-12-31 00:00"])
oneDay = TimeDelta(24.*60.*60.,format='sec')

date = time_range[0]
dates = [date]
print(time_range[0] < time_range[1])
while dates[len(dates)-1] < time_range[1]:
    dates.append(dates[len(dates)-1]+oneDay)
print('dates = ',dates)
print('dates[0] = ',dates[0])
print('type(dates[0]) = ',type(dates[0]))
print('dir(dates[0]) = ',dir(dates[0]))
datesStr = dates[0].strftime("%Y-%m-%d")
print('datesStr = <'+datesStr+'>')
#STOP
#newDate = Time(date)+oneDay
#print('date = ',date,', newDate = ',newDate)
#STOP
#constraints = [AtNightConstraint.twilight_civil()]
#table = observability_table(constraints, observer, [target], time_range=time_range)
#print('table = ',table)
minAltitude = 30.0
goodMinutesPerDate = []
maxGoodMinutes = 0
if False:
    bestDate = None
    for date in dates:
        midnight = date - utcoffset
        utTimes = getDarkTimesUT(midnight)
        print('utTimes = ',utTimes)
        frames = AltAz(obstime=utTimes, location=observatoryLocation)
        targetaltazs = targetCoord.transform_to(frames)
        altitudes = targetaltazs.alt.to_value()
        print('altitudes = ',altitudes)
        nGoodMinutes = len(np.where(altitudes > 30.)[0])
        print('nGoodMinutes = ',nGoodMinutes)
        goodMinutesPerDate.append([date,nGoodMinutes])
        moonDistance = getDistanceToMoon(observatoryLocation, targetCoord, utTimes)
        print('monnDistance = ',moonDistance)
        if nGoodMinutes > maxGoodMinutes:
            maxGoodMinutes = nGoodMinutes
            bestDate = date

    print('bestDate = ',bestDate,': ',maxGoodMinutes,' good minutes')
    bestDateStr = bestDate.strftime("%Y-%m-%d")
    print('bestDateStr = <'+bestDateStr+'>')
visName = '/Users/azuri/daten/uni/HKU/observing_proposals/Pa30_PMAS/visibility_best_dates_start_2023-10-13.png'
maxAltitudeTime = plot_target(targetCoord,
                                observatoryLocation,
                                observatoryName,
                                utcoffset,
                                '2023-10-13',
                                False,
                                outFileName=visName)
maxAltitudeHour = maxAltitudeTime.hour
print('maxAltitudeHour = ',maxAltitudeHour)


if False:
#    target =

    observationStartTime = observer.twilight_evening_astronomical(midnight)#Time('2019-9-5 19:10:00') - utcoffset
    observationEndTime = observer.twilight_morning_astronomical(midnight)#Time('2019-9-6 4:55:00') - utcoffset
    observationStartTime.format = 'iso'
    observationEndTime.format = 'iso'
    print('observationStartTime in local time = ',observationStartTime+utcoffset)
    print('observationEndTime in local time = ',observationEndTime+utcoffset)
    print('observationStartTime in UT = ',observationStartTime)
    print('observationEndTime in UT = ',observationEndTime)


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
