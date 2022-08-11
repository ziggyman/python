import astropy.io.fits as pyfits
import os
import matplotlib.pyplot as plt
from drUtils import getWavelengthArr,getImageData,getHeader,getHeaderValue

path = '/Users/azuri/daten/uni/HKU/interns_projects/julian/'
fileA = os.path.join(path,'6dF_Aug2004_PNG000.2-04.6_6D100804_id20.fits')
fileBBlue = os.path.join(path,'AAOmega_bulge_Sa3-117_BLUEx_300508_id20.fits')
fileBRed = os.path.join(path,'AAOmega_bulge_Sa3-117_RED_2dF300508_id20.fits')

wLenA = getWavelengthArr(fileA,hduNum=0)
fluxA = getImageData(fileA,hduNum=0)

plt.plot(wLenA,fluxA)
plt.show()
