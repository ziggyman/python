#! /usr/bin/env python
import sys
import numpy as np
import astropy.io.fits as pyfits

def getHeader(fName, hduNum=0):
    hdulist = pyfits.open(fName)
    header = hdulist[hduNum].header
    hdulist.close()
    return header

def getImageData(fname,hduNum=0):
    hdulist = pyfits.open(fname)
    scidata = hdulist[hduNum].data
    hdulist.close()
    return np.array(scidata)

if __name__ == '__main__':
    try:
        headerNum = int(sys.argv[2])
    except:
        headerNum = 0
    data = getImageData(sys.argv[1],headerNum)
    header = getHeader(sys.argv[1],headerNum)
    print('fits file name,nPix,min,max,mean,stddev')
    outStr = str(sys.argv[1])+','
    naxis = int(header['NAXIS'])
    nPix = int(header['NAXIS1'])
    for i in range(naxis-1):
        nPix = nPix * int(header['NAXIS%d' % (i+2)])
    outStr += '%d,' % (nPix)
    outStr += '%f,%f,%f,%f' % (np.min(data),np.max(data),np.mean(data),np.std(data))
    print(outStr)
