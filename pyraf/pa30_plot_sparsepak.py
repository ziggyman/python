import astropy.io.fits as pyfits
import os
import numpy as np
import matplotlib.pyplot as plt

fname = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_sum.fits'

hdulist = pyfits.open(fname)
header = hdulist[0].header
data = hdulist[0].data

print('data.shape = ',data.shape)


