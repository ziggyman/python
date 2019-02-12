execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

calibrated = False

gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_sum_n.fits"
wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_sum_n.fits"
wiyn_file = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15T08-22-49_botzfxsEcBld_rebinned_52_n.fits'
if calibrated:
    gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum.fits"
    wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_cal_sum.fits"


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
plt.plot(wiyn_wavelength,wiyn_spectrum,'b-', label='WIYN')
plt.xlabel('wavelength [A]')
yLabel = 'counts [ADU]'
if calibrated:
    yLabel = 'calibrated flux [erg/cm2/s/A]'
plt.ylabel(yLabel)
plt.xlim(4280.,7095.)
plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN+GTC_'
if calibrated:
    plotname = plotname + 'cal_'
plotname = plotname + 'sum_n.png'
plt.savefig(plotname)
plt.show()
