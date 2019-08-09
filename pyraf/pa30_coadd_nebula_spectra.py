import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

from myUtils import getWavelength

calibratedSpectraFile = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/pa30_zdtsEcndr-skyMedian_cal.fits'
calibratedSpectraFile = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/pa30_zd_Ecmsndr-skyMedian_cal.fits'
starSpectra = [0,4,5,6,9,18,20,21,30,33,]

hdulist = pyfits.open(calibratedSpectraFile)
header = hdulist[0].header
wavelength = getWavelength(header,1)
spectrum = hdulist[0].data
print('spectrum = ',spectrum.shape,': ',spectrum)
print('wavelength = ',wavelength.shape,': ',wavelength)

coadd = np.zeros(wavelength.shape[0])
nSpec = 0
for iSpec in np.arange(0,spectrum.shape[0],1):
    specMean = np.mean(spectrum[iSpec,:])
    print('mean of spectrum ',iSpec,' = ',specMean)
    if abs(specMean) < 2.0e-16:
        print('using spectrum ',iSpec)
        coadd += spectrum[iSpec,:]
        nSpec += 1

print('mean of ',nSpec,' spectra in coadd = ',np.mean(coadd))
plt.plot(wavelength,coadd)
plt.show()

hdulist[0].data = coadd
print('hdulist[0].header = ',hdulist[0].header)
hdulist.writeto(calibratedSpectraFile[0:calibratedSpectraFile.rfind('.')]+'_coadd.fits')


gtcFile = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_av_x_wl_flt_cal_mSky_obs_not_smoothed_minComb.fits'
columns = np.arange(1790,1827,1)
hdulist = pyfits.open(gtcFile)
wavelength = getWavelength(header,1)
spectra = hdulist[0].data
print('spectra.shape=',spectra.shape)
spectrum = np.zeros(shape=(spectra.shape[0]))
for row in np.arange(spectra.shape[0]):
    spectrum[row] = np.sum(spectra[row,columns[0]:columns[len(columns)-1]])
hdulist[0].data = spectrum
hdulist.writeto(gtcFile[0:gtcFile.rfind('.')]+'_%d-%d.fits' % (columns[0], columns[len(columns)-1]))

