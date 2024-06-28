from plot_obs_planning import plot_target
from myUtils import hmsToDeg,dmsToDeg
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun

ra = '12:55:34'
dec = '25:53:31'
observatoryName = "SAAO"
observatoryLocation = EarthLocation(lat=-32.3783*u.deg, lon=20.8105*u.deg, height=1750*u.m)
utcoffset = 2*u.hour
date = '2024-06-08'
midnight = Time(date+' 00:00:00') - utcoffset
minAltitude = 30.0

targetCoord = SkyCoord(ra=hmsToDeg(ra)*u.deg,
                        dec=dmsToDeg(dec)*u.deg,
                        frame='icrs')

plot_target(targetCoord,
            observatoryLocation,
            observatoryName,
            utcoffset,
            date)
