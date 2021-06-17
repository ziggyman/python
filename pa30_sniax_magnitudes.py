import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy import pyasl

Mv_max = -19. #mag
Mv_min = -13. #mag

Mv_dd = -11.3

Mv_nova = -8.8

Av = 2.4 #mag


delta_m_15_B_max = 2.4 # mag / 15 days (Stritzinger 2015)
delta_m_15_B_min = 1.25 # mag / 15 days (Stritzinger 2015)

delta_m_15_R_max = 1.0 # mag / 15 days (Jha 2017)
delta_m_15_R_min = 0.2 # mag / 15 days (Jha 2017)

radius = 100. #arcsec
eRadius = 10. #arcsec
expansionVelocity = 1100. #km/s
eExpansionVelocity = 110. #km/s
distance = 2300. # pc
eDistance = 140 # pc

ages = []
m_max = []
m_min = []
m_dd = []
m_nova = []
radii_pc = []

for rad,dist,expVelo in [[radius,distance,expansionVelocity],
                         [radius+eRadius,distance+eDistance,expansionVelocity+eExpansionVelocity],
                         [radius+eRadius,distance+eDistance,expansionVelocity-eExpansionVelocity],
                         [radius+eRadius,distance-eDistance,expansionVelocity+eExpansionVelocity],
                         [radius+eRadius,distance-eDistance,expansionVelocity-eExpansionVelocity],
                         [radius-eRadius,distance+eDistance,expansionVelocity+eExpansionVelocity],
                         [radius-eRadius,distance+eDistance,expansionVelocity-eExpansionVelocity],
                         [radius-eRadius,distance-eDistance,expansionVelocity+eExpansionVelocity],
                         [radius-eRadius,distance-eDistance,expansionVelocity-eExpansionVelocity],
                        ]:
    add = (2.5 * np.log10(np.square(dist) / 100.))
    m_max.append(Mv_max + add + Av)
    m_min.append(Mv_min + add + Av)
    m_dd.append(Mv_dd + add + Av)
    m_nova.append(Mv_nova + add + Av)

    radiusDeg = rad / 3600.
    radiusRad = np.radians(radiusDeg)
    radiusPc = dist * np.sin(radiusRad)
    radii_pc.append(radiusPc)
    radiusKm = radiusPc * 3.086e+13

    ageSec = radiusKm / expVelo #seconds
    ageMin = ageSec / 60.
    ageHours = ageMin / 60.
    ageDays = ageHours / 24.
    ageYrs = ageDays / 365.2425
    ages.append(ageYrs)
    print('dist = ',dist,': ageYrs = ',ageYrs,', m_max = ',m_max[len(m_max)-1],', m_min = ',m_min[len(m_min)-1],', m_dd = ',m_dd[len(m_dd)-1],', radius = ',radii_pc[len(radii_pc)-1],' pc')

ages = np.array(ages)
m_max = np.array(m_max)
m_min = np.array(m_min)
m_dd = np.array(m_dd)
radii_pc = np.array(radii_pc)

print('min(ages) = ',np.min(ages))
print('max(ages) = ',np.max(ages))
print('min(m_min) = ',np.min(m_min))
print('max(m_max) = ',np.max(m_max))
print('min(m_dd) = ',np.min(m_dd),', max(m_dd) = ',np.max(m_dd))
print('min(m_nova) = ',np.min(m_nova),', max(m_nova) = ',np.max(m_nova))

print('age = ',(np.max(ages)+np.min(ages))/2.,'+-',np.max(ages)-((np.max(ages)+np.min(ages))/2.))
print('radius = ',(np.max(radii_pc)+np.min(radii_pc))/2.,'+-',np.max(radii_pc)-((np.max(radii_pc)+np.min(radii_pc))/2.))
print('m_min = ',(np.max(m_min)+np.min(m_min))/2.,'+-',np.max(m_min)-((np.max(m_min)+np.min(m_min))/2.))
print('m_max = ',(np.max(m_max)+np.min(m_max))/2.,'+-',np.max(m_max)-((np.max(m_max)+np.min(m_max))/2.))
print('m_dd = ',(np.max(m_dd)+np.min(m_dd))/2.,'+-',np.max(m_dd)-((np.max(m_dd)+np.min(m_dd))/2.))
print('m_nova = ',(np.max(m_nova)+np.min(m_nova))/2.,'+-',np.max(m_nova)-((np.max(m_nova)+np.min(m_nova))/2.))


# double-degenerate model by Kashyap 2018:
kickVelocity = 90.# km/s
distanceTraveled_km = kickVelocity * 3600. * 24. * 365.2425 * 1000.#km
distanceTraveled_pc = distanceTraveled_km / 3.086e+13
print('distanceTraveled_pc = ',distanceTraveled_pc)

