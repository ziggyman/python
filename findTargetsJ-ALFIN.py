from astroplan import Observer
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u

from myUtils import hmsToDeg,dmsToDeg
import csvFree,csvData

location = EarthLocation(lat=40.0420*u.deg, lon=1.0161*u.deg, height=1957*u.m)

targetListPNeFName = '/Users/azuri/daten/uni/HKU/J-ALPHIN/targetsPNe.csv'
targetListPNeHaloesFName = '/Users/azuri/daten/uni/HKU/J-ALPHIN/targetsPNeHaloes.csv'

targetsPNe = csvFree.readCSVFile(targetListPNeFName)
targetsPNHaloes = csvFree.readCSVFile(targetListPNeHaloesFName)

for targetList in [targetsPNe,targetsPNHaloes]:
    for i in range(targetList.size()):
        ra = hmsToDeg(targetList.getData('RA',i).replace(' ',':'))
        dec = dmsToDeg(targetList.getData('DEC',i).replace(' ',':'))
