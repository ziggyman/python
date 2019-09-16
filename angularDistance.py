#import myUtils
from astropy.coordinates import SkyCoord
#import astropy.units as u
#18:03:13.90 -46:15:10.08
mm1 = SkyCoord(ra='18h03m13.90s', dec='-46d15m10.08s', frame='icrs')
mm2 = SkyCoord(ra='18h03m13.992s', dec='-46d15m39.379s', frame='icrs')
sep = mm1.separation(mm2).arcsecond
print('separation = ',sep)