import numpy as np
import matplotlib.pyplot as plt

t = np.arange(1,1000,1)
r = t ** (2/5)

v = r/t

#plt.plot(t,r,label = 'r')
plt.plot(t,v,label = 'v')
plt.legend()
plt.xlabel('t')
plt.ylabel('v')
plt.show()


E = (10. ** 49) / (10.**51.) # erg Explosion energy in units of 10 ** 51 erg
n = 1. # cm ** -3
r = 2.2 # pc

k = 1.3806 * (10. ** (-23.)) # m^2 kg s^-2 K^-1
T = 20.E6 #K

kT = 25.7E-6 * 1.E-9

ageE = (r / (0.31 * (E/n)**.2)) ** (5./2.)
print('ageE = ',ageE)

ageT = (1.8E5 * r**2 / (kT)) ** 0.5
print('ageT = ',ageT)
