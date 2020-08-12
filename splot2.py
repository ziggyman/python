from astropy.nddata import CCDData
from astropy import units as u
import matplotlib.pyplot as plt
import os

from drUtils import getWavelengthArr,readFileToArr

objectSpectraInA = readFileToArr('/Users/azuri/spectra/saao/saao_sep2019/20190906/SCIENCE_otzfifEc.list')
objectSpectraInB = readFileToArr('/Users/azuri/spectra/saao/saao_sep2019/20190906/SCIENCE_otzfifEc_oldExtraction.list')
for iSpec in range(len(objectSpectraInA)):
    imgA = CCDData.read(objectSpectraInA[iSpec], unit=u.adu)
    imgB = CCDData.read(objectSpectraInB[iSpec], unit=u.adu)
    wLenA = getWavelengthArr(objectSpectraInA[iSpec],0) * u.angstrom
    wLenB = getWavelengthArr(objectSpectraInB[iSpec],0) * u.angstrom
    plt.plot(wLenA, imgA.data, label=objectSpectraInA[iSpec][objectSpectraInA[iSpec].rfind('/')+1:objectSpectraInA[iSpec].rfind('.')])
    plt.plot(wLenB, imgB.data, label=objectSpectraInB[iSpec][objectSpectraInB[iSpec].rfind('/')+1:objectSpectraInB[iSpec].rfind('.')])
    plt.legend()
    plt.xlabel('wavelength [$\AA$]')
    plt.ylabel('flux')
    plt.show()
