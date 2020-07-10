from astropy.nddata import CCDData
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import os

from drUtils import getWavelengthArr

path = '/Volumes/work/azuri/spectra/saao/saao-nov20160-reduced/02nov/final'
objectSpectraIn = ['NGC246norm_SA021116.fits','NGC246norm_SA021116B.fits']

for iSpec in range(len(objectSpectraIn)):
    img = CCDData.read(os.path.join(path, objectSpectraIn[iSpec]), unit=u.adu)
    wLen = getWavelengthArr(os.path.join(path, objectSpectraIn[iSpec]),0) * u.angstrom
    print('img.data = ',img.data[0][0])

    plt.plot(wLen, img.data[0][0], label=objectSpectraIn[iSpec][:objectSpectraIn[iSpec].rfind('.')])
plt.legend()
plt.xlabel('wavelength [$\AA$]')
plt.ylabel('flux')
plt.show()
