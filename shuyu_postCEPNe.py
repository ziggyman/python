import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.stats import circstats
import pycircstat

import csvFree,csvData
from pnOrientationUtils import calcMoments,calcMean,rose_plot,vectorDiagram,plotHammerProjection,linearOrderDiagram,plotEPA
from myUtils import hmsToDeg,dmsToDeg,findClosestObjectTo

dataFileName = '/Users/azuri/daten/uni/HKU/publications/shuyu/close_binary_GPAs.txt'
with open(dataFileName,'r') as f:
    lines = f.readlines()
GPAs = [float(line.strip()) for line in lines]
print('GPAs = ',len(GPAs),': ',GPAs)

gpaRad = np.radians(2.*np.array(GPAs))
pRad = circstats.rayleightest(gpaRad)
print('pRad = ',pRad)

moments = calcMoments(GPAs)
if moments[0][0] < 0.:
    moments[0][0] += 180.
print('moments = ',moments)

rosePlotMean = moments[0][0]
rosePlotSigma = moments[1][0]

fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
rose_plot(ax, GPAs, bins=16, density=True, offset=0, lab_unit="degrees",
            start_zero=True, fill=True, color='blue', max_count=None,
            max_size=None, smooth=False, mean=rosePlotMean, sigma=rosePlotSigma)
fig.tight_layout()
fig.savefig('/Users/azuri/daten/uni/HKU/publications/shuyu/post-CEbinaryPNe_rosePlot.png', bbox_inches='tight')
plt.show()
print('rosePlotMean = ',rosePlotMean,', rosePlotSigma = ',rosePlotSigma)

vectorDiagram(2.*np.array(GPAs), fNameOut='/Users/azuri/daten/uni/HKU/publications/shuyu/post-CEbinaryPNe_vectorDiagram.png')
linearOrderDiagram(2.*np.array(GPAs), fNameOut='/Users/azuri/daten/uni/HKU/publications/shuyu/post-CEbinaryPNe_linearOrderDiagram.png')
