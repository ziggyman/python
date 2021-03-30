import matplotlib.pyplot as plt
import numpy as np

import csvFree,csvData

distance = 2300. #pc
Mv_max = -19. #mag
Mv_min = -13. #mag

add = (2.5 * np.log10(np.square(distance) / 100.))
mv_max = Mv_max + add
mv_min = Mv_min + add
print('mv_max = ',mv_max,', mv_min = ',mv_min)

delta_m_15_B_max = 2.4 # mag / 15 days (Stritzinger 2015)
delta_m_15_B_min = 1.25 # mag / 15 days (Stritzinger 2015)

delta_m_15_R_max = 1.0 # mag / 15 days (Jha 2017)
delta_m_15_R_min = 0.2 # mag / 15 days (Jha 2017)

time = np.arange(0,300,0.1)

table3Name = '/Users/azuri/daten/uni/HKU/Pa30/Stritzinger2015/table3.dat'
table3 = csvFree.readCSVFile(table3Name,'\t',False)

print('table3.header = ',table3.header)
dates = 2450000.+np.array([float(i) for i in table3.getData('JD2450000+')])

Bmag = np.array(table3.getData('B (mag)'))
print('Bmag = ',Bmag)
print('Bmag[1] = ',Bmag[1])
idx = np.where(Bmag != '-')
print('idx = ',idx)
print('idx[0] = ',idx[0])
Bdates = dates[idx]
BmagPlot = Bmag[idx]
BmagPlot = [float(b[:b.find('(')]) for b in BmagPlot]
Bdates = Bdates-(Bdates[np.where(BmagPlot == np.min(BmagPlot))[0]])
plt.scatter(Bdates,BmagPlot,c='b',label='B')

Vmag = np.array(table3.getData('V (mag)'))
print('Vmag = ',Vmag)
idx = np.where(Vmag != '-')
print('idx = ',idx)
print('idx[0] = ',idx[0])
Vdates = dates[idx]
VmagPlot = Vmag[idx]
VmagPlot = [float(b[:b.find('(')]) for b in VmagPlot]
Vdates = Vdates-(Vdates[np.where(VmagPlot == np.min(VmagPlot))[0]])
plt.scatter(Vdates,VmagPlot,c='m',label='V')

rmag = np.array(table3.getData('r (mag)'))
print('rmag = ',rmag)
idx = np.where(rmag != '-')
print('idx = ',idx)
print('idx[0] = ',idx[0])
rdates = dates[idx]
rmagPlot = rmag[idx]
rmagPlot = [float(b[:b.find('(')]) for b in rmagPlot]
rdates = rdates-(rdates[np.where(rmagPlot == np.min(rmagPlot))[0]])
plt.scatter(rdates,rmagPlot,c='r',label='r')

Vfit = np.polyfit(Vdates[:-1],VmagPlot[:-1], 6)
Vx = np.arange(Vdates[0],Vdates[len(Vdates)-2],1)
Vy = np.poly1d(Vfit)(Vx)
plt.plot(Vx,Vy)

plt.ylim([np.min(rmagPlot),np.max(BmagPlot)])
plt.legend()
plt.gca().invert_yaxis()
plt.show()

radius = 100. #arcsec
expansionVelocity = 1100. #km/s
distance = 2.3 # +- 0.4 kpc
eDistance = 0.4

radiusDeg = radius / 3600.
radiusRad = np.radians(radiusDeg)

for dist in [distance,distance+eDistance,distance-eDistance]:
    radiusPc = dist * 1000. * np.sin(radiusRad)
    radiusKm = radiusPc * 3.086e+13
    print('dist = ',dist,': radiusPc = ',radiusPc,', radiusKm = ',radiusKm)
    ageSec = radiusKm / expansionVelocity #seconds
    print('dist = ',dist,': ageSec = ',ageSec)
    ageMin = ageSec / 60.
    ageHours = ageMin / 60.
    ageDays = ageHours / 24.
    ageYrs = ageDays / 365.2425
    print('dist = ',dist,': ageYrs = ',ageYrs)
    print('dist = ',dist,': age in Yrs = ',ageSec * 3.154e7)
