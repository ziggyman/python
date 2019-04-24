import astropy.io.fits as pyfits
import numpy as np

inFile = '/Users/azuri/daten/uni/HKU/Pa30/KittPeak/Pa30_OIII_astro_sum.fits'
outFile = '/Users/azuri/daten/uni/HKU/Pa30/KittPeak/Pa30_OIII_astro_sum_binned.fits'

hdulist = pyfits.open(inFile)
header = hdulist[0].header
inIm = hdulist[0].data

outIm = np.zeros([int(inIm.shape[0]/4), int(inIm.shape[1]/4)])

for i in range(outIm.shape[0]):
    for j in range(outIm.shape[1]):
        outIm[i,j] = np.sum(inIm[i*4:(i+1)*4, j*4:(j+1)*4])

pyfits.writeto(outFile, outIm, header)
