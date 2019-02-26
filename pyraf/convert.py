#!/usr/bin/python
import os
import sys
from astropy.io import fits as pyfits
import numpy as np
from pyraf import iraf
#
#infile = "/Users/Lea/Desktop/CS_centre_sum_spec.fits"

infile = "/Volumes/obiwan/azuri/spectra/Kamila/CSTotal.fits"
hdulist = pyfits.open( infile )

# print column information
hdulist[1].columns

# get to the data part (in extension 1)
scidata = hdulist[1].data

print('type(scidata[0]) = ',type(scidata[0]))
print('dir(scidata[0] = ',dir(scidata[0]))
print('scidata[0].array = ',scidata[0].array)
print('type(scidata[0].array) = ',type(scidata[0].array))
print('dir(scidata[0].array) = ',dir(scidata[0].array))

arr = scidata[0].array
print('arr.columns = ',arr.columns)
print('arr.names = ',arr.names)
print('arr.data = ',arr.data)
#STOP
#myList = arr.tolist()
#print('myList = ',myList)

wave = None
flux = None
if len(arr.names) > 2:
    wave = scidata[0][arr.names[0]]#'wave']
    flux = scidata[0][arr.names[1]]#'flux']
    try:
        print('len(wave) = ',len(wave))
        print('len(flux) = ',len(flux))
    except:
        print('something went wrong, trying something different')
        wave = arr['wave']
        flux = arr['flux']
        print('wave = ',wave)
        print('flux = ',flux)
        print('len(wave) = ',len(wave))
        print('len(flux) = ',len(flux))
else:
    wave = arr[arr.names[0]]#myList[:][0]#scidata[0]['wave']
    flux = arr[arr.names[1]]#myList[:][1]#scidata[0]['flux']
#arr2 = scidata[0]['ERR_REDUCED']
# etc.
# where flux will contain the data corresponding to the column named: hdulist[1].columns[1]
# where arr2 will contain the data corresponding to the column named: hdulist[1].columns[2]
# etc.
print('wave = ',wave)
print('flux = ',flux)
# To plot using maptplotlib:

#import matplotlib.pyplot as plt

#plt.plot(wave, flux)

#plt.show()

outFile = infile[0:infile.rfind('.')]
if os.path.exists(outFile+'.text'):
    os.remove(outFile+'.text')
f = open(outFile+'.text','w')
for line in range(len(wave)):
    f.write(str(wave[line])+' '+str(flux[line])+'\n')
f.close()

if os.path.exists(outFile+'_new.fits'):
    os.remove(outFile+'_new.fits')
iraf.rspectext(outFile+'.text',outFile+'_new.fits',dtype="interp")

