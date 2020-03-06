import numpy as np
import matplotlib.pyplot as plt

from drUtils import gauss

amplitudes = np.arange(100.,1000.,1.)
sigmas = [1.,3.,10.]
thresholds = [10,50,100]
x = np.arange(-100.,100.,1.)

for threshold in thresholds:
    for sigma in sigmas:
        cog = []
        for amplitude in amplitudes:
            g = gauss(x,amplitude,0.,sigma,0.)
            indices = np.where(g > threshold)
            print('indices = ',indices)
            cog.append(indices[0].shape[0])
        plt.plot(amplitudes,cog)
        plt.show()