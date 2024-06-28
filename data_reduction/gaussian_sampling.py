import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import simps
from drUtils import gauss

x = np.arange(-5,5.01,0.01)
y = gauss(x,1.,0.,1.5)
plt.plot(x,y)
print('x = ',len(x),': ',x)
print('y = ',len(y),': ',y)

x_bins = np.arange(-4.5,5.5,1)
y_bins = np.array([])
integ = np.array([])
for i in range(len(x_bins)):
    idx = np.where(x >= x_bins[i]-0.5)[0]
    print('idx = ',idx)
    idxb = np.where(x[idx] < x_bins[i]+0.5)[0]
    print('idxb = ',idxb)
    print('idx[idxb] = ',idx[idxb])
    print('y[idx[idxb]] = ',y[idx[idxb]])
    y_bins = np.append(y_bins,simps(y[idx[idxb]]) * 0.01)
    integ = np.append(integ,np.sum(y[idx[idxb]] * 0.01))
plt.scatter(x_bins,y_bins)
plt.scatter(x_bins,integ)

xPlus = np.arange(-0.5,0.51,0.01)
for xp in xPlus:
    x_binsp = x_bins + xp
    y_binsp = np.array([])
    for i in range(len(x_bins)):
        idx = np.where(x >= x_binsp[i]-0.5)[0]
        print('idx = ',idx)
        idxb = np.where(x[idx] < x_binsp[i]+0.5)[0]
        print('idxb = ',idxb)
        print('idx[idxb] = ',idx[idxb])
        print('y[idx[idxb]] = ',y[idx[idxb]])
        y_binsp = np.append(y_binsp,simps(y[idx[idxb]]) * 0.01)
        integ = np.append(integ,np.sum(y[idx[idxb]] * 0.01))
    plt.scatter(x_binsp,y_binsp)
plt.show()
