import numpy as np
import pyneb as pn
import matplotlib.pyplot as plt

S2 = pn.Atom('S',2)

ratios = np.arange(0.5,2.,0.05)

dens = S2.getTemDen(int_ratio=ratios,tem=10000.,wave1=6716,wave2=6731,maxIter=100)

print('dens = ',dens)
dens[np.argwhere(np.isnan(dens))] = 0.

plt.plot(ratios,dens)
plt.show()
