#! /usr/bin/env python
import sys
import astropy.io.fits as pyfits

def getHeader(fName, hduNum=0):
    hdulist = pyfits.open(fName)
    header = hdulist[hduNum].header
    hdulist.close()
    return header

if __name__ == '__main__':
    try:
        headerNum = int(sys.argv[2])
    except:
        headerNum = 0
    header = getHeader(sys.argv[1],headerNum)
    #print(header.keys)
    outStr = str(sys.argv[1])+'['
    naxis = int(header['NAXIS'])
    outStr += '%d' % (header['NAXIS1'])
    for i in range(naxis-1):
        outStr += ',%d' % (header['NAXIS%d' % (i+2)])
    outStr += ']: ' + header['IMAGETYP'] + ': ' + header['OBJECT'] + ' %ds' % (header['EXPTIME'])
    print(outStr)
