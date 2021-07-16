import numpy as np
import matplotlib.pyplot as plt

mSolar = 1.989e30 #kg
mStar = 1.5 * mSolar
oxygen16Mass = 2.66e-26# kg
atomicMass = oxygen16Mass / 15.99491461956
sulfurAtomMass = 32.066 * atomicMass# kg
ejectionSpeed = 1100000# m/s
nHours = 1000. * 365.25 * 24.# * 3600.# s
G = 6.6743e-11#m^3 * kg^-1 * s^-2

print('atomicMass = ',atomicMass)
print('nHours = ',nHours)

radius = np.zeros(int(nHours))
radius[0] = 10000000.
velocity = np.zeros(int(nHours))
velocity[0] = ejectionSpeed
time = np.arange(0,nHours,1)
for i in np.arange(1,nHours,1):
    radius[int(i)] = radius[int(i-1)] + (velocity[int(i-1)] * 3600.)
    print('radius[',i,'] = ',radius[int(i)])
    force = G * mStar * sulfurAtomMass / (radius[int(i-1)]**2)
    print('force = ',force)
    acceleration = force / sulfurAtomMass
    print('acceleration = ',acceleration)
    velocity[int(i)] = velocity[int(i-1)] - (acceleration * 3600.)
    print('velocity[',i,'] = ',velocity[int(i)])
    if i > 10:
        STOP

plt.plot(time,velocity)
plt.show()
plt.plot(time,radius)
plt.show()
