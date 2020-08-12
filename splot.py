from astropy.nddata import CCDData
from astropy import units as u
import matplotlib.pyplot as plt
import os

from drUtils import getWavelengthArr,readFileToArr

#path = '/Volumes/work/azuri/spectra/saao/saao-nov20160-reduced/02nov/final'
objectSpectraIn = readFileToArr('/Users/azuri/spectra/saao/saao_sep2019/20190906/SCIENCE_otzfifEc.list')#['NGC246norm_SA021116.fits','NGC246norm_SA021116B.fits']
print('objectSpectraIn = ',objectSpectraIn)
for iSpec in range(len(objectSpectraIn)):
    img = CCDData.read(objectSpectraIn[iSpec], unit=u.adu)
    print('img = ',img)
    print('img.data = ',img.data)
    wLen = getWavelengthArr(objectSpectraIn[iSpec],0) * u.angstrom
    plt.plot(wLen, img.data, label=objectSpectraIn[iSpec][objectSpectraIn[iSpec].rfind('/')+1:objectSpectraIn[iSpec].rfind('.')])
plt.legend()
plt.xlabel('wavelength [$\AA$]')
plt.ylabel('flux')
plt.show()
