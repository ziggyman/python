from astroplan import Observer, FixedTarget
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u
from astroplan import AltitudeConstraint, AirmassConstraint,AtNightConstraint
from astroplan import observability_table

from myUtils import hmsToDeg,dmsToDeg
import csvFree,csvData

siteNames = EarthLocation.get_site_names()
for siteName in siteNames:
    print(siteName)

lon = 1.0161*u.deg
lat = 40.0420*u.deg
elev = 1957*u.m
#location = EarthLocation(lat=, lon=lon, height=elev)
observer = Observer(longitude=lon, latitude=lat, elevation=elev, name="Javalambre")
time_range = Time(["2023-04-01 06:00", "2024-04-01 06:00"])
airmassConstraints = [1.5,1.4,1.3,1.2]

targetListPNeFName = '/Users/azuri/daten/uni/HKU/J-ALFIN/targetsPNe.csv'
targetListPNeHaloesFName = '/Users/azuri/daten/uni/HKU/J-ALFIN/targetsPNeHaloes.csv'

hashAngDiamsFName = '/Users/azuri/daten/uni/HKU/J-ALFIN/hash_tbAngDiam.csv'
hashCNamesFName = '/Users/azuri/daten/uni/HKU/J-ALFIN/hash_tbCNames.csv'
hashPNMainFName = '/Users/azuri/daten/uni/HKU/J-ALFIN/hash_PNMain.csv'

hashAngDiams = csvFree.readCSVFile(hashAngDiamsFName)
hashCNames = csvFree.readCSVFile(hashCNamesFName)
hashPNMain = csvFree.readCSVFile(hashPNMainFName)

for targetListFName in [targetListPNeFName,targetListPNeHaloesFName]:
    targetList = csvFree.readCSVFile(targetListFName)
    targetList.addColumn('AngDiam [arcsec]')

    targets = []
    for i in range(targetList.size()):
        name = targetList.getData('Name',i)
        ra = hmsToDeg(targetList.getData('RA',i).replace(' ',':'))
        dec = dmsToDeg(targetList.getData('DEC',i).replace(' ',':'))
        targets.append(FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name))

        targetList.setData('RA',i,targetList.getData('RA',i).replace(' ',':'))
        targetList.setData('DEC',i,targetList.getData('DEC',i).replace(' ',':'))
        idxHashCNames = hashCNames.find('Name',name)
        if idxHashCNames[0] < 0:
            idxHashCNames = hashCNames.find('Name',name.replace(' ',''))
        if idxHashCNames[0] < 0:
            print('did not find name <'+name+'> in hashCNames')
            STOP
        idPNMain = hashCNames.getData('idPNMain',idxHashCNames[0])
        idxHashAngDiams = hashAngDiams.find('idPNMain',idPNMain)
        for idx in idxHashAngDiams:
            if hashAngDiams.getData('InUse',idx) == '1':
                targetList.setData('AngDiam [arcsec]',i,hashAngDiams.getData('MajDiam',idx))

    for airmass in airmassConstraints:
        fNameOut = targetListFName[:-4]+'_airmass_lt_%.1f.txt' % (airmass)
        constraints = [AirmassConstraint(airmass), AtNightConstraint.twilight_civil()]

        table = observability_table(constraints, observer, targets, time_range=time_range)
        print('targetListFName = ',targetListFName,': airmass = ',airmass)
        table.remove_column('always observable')
        table.add_column(targetList.getData('RA'),name='RA')
        table.add_column(targetList.getData('DEC'),name='DEC')
        table.add_column(targetList.getData('priority'),name='priority')
        table.add_column(targetList.getData('AngDiam [arcsec]'),name='AngDiam [arcsec]')
        print('table = ',table)
        table.write(fNameOut, format='ascii', overwrite=True)
