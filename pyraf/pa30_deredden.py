import PyAstronomy.pyasl
execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned_scaled_absoluteFlux.fits"
wiyn_file = "/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned_absoluteFlux.fits"#Pa30_WN151014_cal_sum_cleaned.fits"
leDu_file = "/Users/azuri/daten/uni/HKU/Pa30/LeDu/_pa30_20181009_941_PLeDu_scaled.fits"
somme_file = "/Users/azuri/daten/uni/HKU/Pa30/_pa30_somme6_scaled.fits"

gtc_hdulist = pyfits.open(gtc_file)
wiyn_hdulist = pyfits.open(wiyn_file)
leDu_hdulist = pyfits.open(leDu_file)
somme_hdulist = pyfits.open(somme_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum_absolute = fits.getdata(gtc_file)
print('gtc_wavelength = ',gtc_wavelength)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1)
wiyn_spectrum_absolute = fits.getdata(wiyn_file)
print('wiyn_wavelength = ',wiyn_wavelength)

leDu_header = leDu_hdulist[0].header
leDu_wavelength = getWavelength(leDu_header,1)
leDu_spectrum = fits.getdata(leDu_file)
leDu_spectrum_absolute = calibratedFluxToAbsoluteFlux(leDu_spectrum, 3370.0)
print('leDu_wavelength = ',leDu_wavelength)

somme_header = somme_hdulist[0].header
somme_wavelength = getWavelength(somme_header,1)
somme_spectrum = fits.getdata(somme_file)
somme_spectrum_absolute = calibratedFluxToAbsoluteFlux(somme_spectrum, 3370.0)
print('somme_wavelength = ',somme_wavelength)

leDu_spectrum_smoothed = boxCarMeanSmooth(leDu_spectrum, 0, 21)
somme_spectrum_smoothed = boxCarMeanSmooth(somme_spectrum, 0, 21)
leDu_spectrum_absolute_smoothed = boxCarMeanSmooth(leDu_spectrum_absolute, 0, 21)
somme_spectrum_absolute_smoothed = boxCarMeanSmooth(somme_spectrum_absolute, 0, 21)

ebv = 0.79
R_V = 3.0

gtc_spectrum_absolute_dereddened = pyasl.unred(gtc_wavelength, gtc_spectrum_absolute, ebv=ebv, R_V = R_V)
plt.plot(gtc_wavelength, np.log10(gtc_spectrum_absolute), 'r-', label='GTC absolute scaled')
plt.plot(gtc_wavelength, np.log10(gtc_spectrum_absolute_dereddened), 'b-', label='GTC absolute scaled dereddened')
#plt.xlabel('wavelength [\AA]')
#plt.ylabel('log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])')
#plt.legend()
#plt.show()

wiyn_spectrum_absolute_dereddened = pyasl.unred(wiyn_wavelength, wiyn_spectrum_absolute, ebv=ebv, R_V = R_V)
plt.plot(wiyn_wavelength, np.log10(wiyn_spectrum_absolute), 'g-', label='WIYN absolute')
plt.plot(wiyn_wavelength, np.log10(wiyn_spectrum_absolute_dereddened), 'c-', label='WIYN absolute dereddened')
plt.xlabel('wavelength [\AA]')
plt.ylabel('log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC_absolute_dereddened.pdf'
#plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()


