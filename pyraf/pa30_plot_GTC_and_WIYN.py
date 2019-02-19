execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

calibrated = True
useMean = True

gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_sum_n.fits"
wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_sum_n.fits"
wiyn_file = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15T08-22-49_botzfxsEcBld_rebinned_52_n.fits'
if calibrated:
    gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned.fits"
    wiyn_file = "/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned.fits"#Pa30_WN151014_cal_sum_cleaned.fits"
    if not useMean:
        wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMedian_52_cal_cleaned.fits"#Pa30_WN151014_cal_sum_cleaned.fits"


gtc_hdulist = pyfits.open(gtc_file)
wiyn_hdulist = pyfits.open(wiyn_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum = fits.getdata(gtc_file)
print('gtc_wavelength = ',gtc_wavelength)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1)
wiyn_spectrum = fits.getdata(wiyn_file)
print('wiyn_wavelength = ',wiyn_wavelength)

plt.plot(gtc_wavelength,gtc_spectrum,'g-', label='GTC')
plt.plot(wiyn_wavelength[5:len(wiyn_wavelength)-5],wiyn_spectrum[5:len(wiyn_wavelength)-5],'b-', label='WIYN')
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'counts [ADU]'
if calibrated:
    yLabel = 'calibrated flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plt.xlim(3600.,7900.)
plt.ylim(0.,4.9e-15)
plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC_'
if calibrated:
    plotname = plotname + 'cal_'
plotname = plotname + 'sum.pdf'
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()


gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_rebinned_4000px_cleaned.fits"#Pa30_GT080716_cal_sum_cleaned.fits"
wiyn_file = "/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned_scaled_with_constant.fits"#Pa30_WN151014_cal_sum_rebinned_scaled_down_cleaned.fits"

gtc_hdulist = pyfits.open(gtc_file)
wiyn_hdulist = pyfits.open(wiyn_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum = fits.getdata(gtc_file)
print('gtc_wavelength = ',gtc_wavelength)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1)
wiyn_spectrum = fits.getdata(wiyn_file)
print('wiyn_wavelength = ',wiyn_wavelength)

plt.plot(gtc_wavelength,gtc_spectrum,'g-', label='Pa30 GTC')
plt.plot(wiyn_wavelength,wiyn_spectrum,'b-', label='Pa30 WIYN scaled down')

plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'calibrated flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plt.xlim(4280.,7095.)
plt.ylim(0.0, 1.5e-15)
plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WNscaled_with_constant+GTC.pdf'
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()
