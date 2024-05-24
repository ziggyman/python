import astropy.units as u
import astropy.time
import astropy.coordinates

from myUtils import raDecToLonLat,lonLatToRaDec,hmsToDeg,dmsToDeg

#ecliptic_angle_b = 60.2 #degrees
#ecliptic_angle_l = 120. #degrees
ecliptic_angle_ra = '18:00:00'
ecliptic_angle_dec = '66:34:00'
moon_inclination = 5.14 # degrees
moon_axial_tilt_to_orbit = 6.68 # degrees
moon_obliquity_to_ecliptic = 1.54 # degrees

galactic_center_RA = '17:45:40.04'
galactic_center_DEC = '-29:00:28.1'

galactic_center_l,galactic_center_b = raDecToLonLat(hmsToDeg(galactic_center_RA),dmsToDeg(galactic_center_DEC))
print('galactic_center_l = ',galactic_center_l,', galactic_center_b = ',galactic_center_b)

earth_axis_ra = 0.
earth_axis_dec = 90.

l,b = raDecToLonLat(earth_axis_ra,earth_axis_dec)
print('l = ',l,', b = ',b)

earth_axis_ra = 180.
earth_axis_dec = 90.

ecliptic_angle_l,ecliptic_angle_b = raDecToLonLat(hmsToDeg(ecliptic_angle_ra),dmsToDeg(ecliptic_angle_dec))
print('ecliptic_angle_l = ',ecliptic_angle_l,', ecliptic_angle_b = ',ecliptic_angle_b)
STOP
ra,dec = lonLatToRaDec(ecliptic_angle_l,ecliptic_angle_b)
print('ra = ',ra,', dec = ',dec)

import numpy as np
import math
import astropy.units as u
from astropy.coordinates import SkyCoord,spherical_to_cartesian,cartesian_to_spherical

ra, dec = lonLatToRaDec(ecliptic_angle_l,ecliptic_angle_b)
c1 = SkyCoord(ra=ra*u.degree, dec = dec*u.degree, frame='fk5')
c2 = SkyCoord(ra=hmsToDeg(galactic_center_RA)*u.degree,dec=dmsToDeg(galactic_center_DEC)*u.degree)

print('separation from Skycoord = ',c1.position_angle(c2).to(u.deg))

x1,y1,z1 = spherical_to_cartesian(1.,math.radians(ecliptic_angle_b),math.radians(ecliptic_angle_l))
print('x1 = ',x1,', y1 = ',y1,', z1 = ',z1)

x2,y2,z2 = spherical_to_cartesian(1.,0.,math.radians(180.))
print('x2 = ',x2,', y2 = ',y2,', z2 = ',z2)

dotProduct = np.dot([x1,y1,z1],[x2,y2,z2])
print('dotProduct = ',dotProduct)
print('dir(dotProduct) = ',dir(dotProduct))
sq1 = np.sqrt(x1**2+y1**2+z1**2)
sq2 = np.sqrt(x2**2+y2**2+z2**2)
print('sq1 = ',sq1)
print('dir(sq1) = ',dir(sq1))
print('sq2 = ',sq2)
sep = math.degrees(np.arccos(sq1.value*sq2.value / dotProduct))
print('sep = ',sep)
