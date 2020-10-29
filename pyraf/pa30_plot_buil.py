from astropy.io import fits
import astropy.io.fits as pyfits
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

from myUtils import getWavelength
#execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

#ebv = 1.206
ebv = 0.79 # 3D map Bayestar15
#ebv = 0.99 # 3D map Bayestar17
R_V = 3.1  # (Fitzpatrick and Massa 2007)
#R_V = 3.1
A_V = 2.4 #Claire
ebv = A_V / R_V

somme_file = "/Users/azuri/daten/uni/HKU/Pa30/_pa30_somme6_scaled.fits"
somme_hdulist = pyfits.open(somme_file)
somme_header = somme_hdulist[0].header
somme_wavelength = getWavelength(somme_header,1)
somme_spectrum = fits.getdata(somme_file)
print('somme_wavelength = ',somme_wavelength)

somme_spectrum_smoothed = ndimage.median_filter(somme_spectrum, size=21)
somme_spectrum_smoothed_dereddened = pyasl.unred(somme_wavelength, somme_spectrum_smoothed, ebv=ebv, R_V = R_V)
plt.plot(somme_wavelength,np.log10(somme_spectrum_smoothed),'m-', label='Buil scaled and smoothed')
plt.plot(somme_wavelength,np.log10(somme_spectrum_smoothed_dereddened),'m-')#, label='Somme scaled, smoothed, and dereddened')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/pa30_buil_uvuex.pdf'
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()
