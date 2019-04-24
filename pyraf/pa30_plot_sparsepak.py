import astropy.io.fits as pyfits
import os
import numpy as np
import matplotlib.pyplot as plt
from myUtils import getWavelength

fname = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_all_fibres_sum.fits'

hdulist = pyfits.open(fname)
header = hdulist[0].header
data = hdulist[0].data
wavelength = getWavelength(header, 1)

print('data.shape = ',data.shape)

for i in range(data.shape[0]):
    plt.plot(wavelength, data[i,:])
plt.show()
