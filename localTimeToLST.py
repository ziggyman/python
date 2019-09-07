import ephem
from astroplan import Observer
#from astroplan import download_IERS_A
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

gatech = ephem.Observer()
#ssoObs = Observer.at_site("Siding Spring Observatory")#, timezone='Eastern Standard Time')
sso_utcoffset = 2*u.hour  # Eastern Daylight Time
gatech.lat, gatech.lon = '-31.2749', '149.0685'

saaoObs = Observer.at_site("Southern African Large Telescope")#, timezone='Eastern Daylight Time')
lat = -32.3783
lon = 20.8105
saao = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=1750*u.m)
saao_utcoffset = 2*u.hour  # Eastern Daylight Time

#localTime = Time('2019-9-6 19:41:00') - saao_utcoffset
localTime = Time('2019-9-6 19:41:00') - sso_utcoffset

#gatech.lat, gatech.lon = str(lat), str(lon)#'-32.3783', '20.8105'

timeStr = localTime.strftime("%y/%m/%d %H:%M")
print('timeStr = ',timeStr)
gatech.date = timeStr#'2019/9/6 1:30'
print(gatech.sidereal_time())
