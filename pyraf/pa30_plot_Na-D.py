execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...


calibrated = True

gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum.fits"
wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_cal_sum_rebinned_scaled_down.fits"
hiltner102_before_file = "/Volumes/obiwan/azuri/spectra/sparsepak/2014B-0176/Hiltner102_before_combined_cal.fits"
hiltner102_after_file = "/Volumes/obiwan/azuri/spectra/sparsepak/2014B-0176/Hiltner102_after_combined_cal.fits"
fibre78_file = "/Volumes/obiwan/azuri/spectra/sparsepak/2014B-0176/kwb_141015_080034_ori_otxzd_78_cal.fits"

gtc_hdulist = pyfits.open(gtc_file)
wiyn_hdulist = pyfits.open(wiyn_file)
hiltner102_before_hdulist = pyfits.open(hiltner102_before_file)
hiltner102_after_hdulist = pyfits.open(hiltner102_after_file)
fibre78_hdulist = pyfits.open(fibre78_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum = fits.getdata(gtc_file)
print('gtc_wavelength = ',gtc_wavelength)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1)
wiyn_spectrum = fits.getdata(wiyn_file)
print('wiyn_wavelength = ',wiyn_wavelength)

hiltner102_before_header = hiltner102_before_hdulist[0].header
hiltner102_before_wavelength = getWavelength(hiltner102_before_header,1)
hiltner102_before_spectrum = fits.getdata(hiltner102_before_file)

hiltner102_after_header = hiltner102_after_hdulist[0].header
hiltner102_after_wavelength = getWavelength(hiltner102_after_header,1)
hiltner102_after_spectrum = fits.getdata(hiltner102_after_file)

fibre78_header = fibre78_hdulist[0].header
fibre78_wavelength = getWavelength(fibre78_header,1)
fibre78_spectrum = fits.getdata(fibre78_file)


plt.plot(gtc_wavelength,gtc_spectrum,'g-', label='Pa30 GTC')
plt.plot(wiyn_wavelength,wiyn_spectrum,'b-', label='Pa30 WIYN scaled down')
plt.plot(hiltner102_before_wavelength,hiltner102_before_spectrum / 20.,'r-', label='Hiltner102 before / 20')
plt.plot(hiltner102_after_wavelength,hiltner102_after_spectrum,'c-', label='Hiltner102 after')
plt.plot(fibre78_wavelength,fibre78_spectrum,'m-', label='fibre 78')

plt.xlabel('wavelength [A]')
yLabel = 'calibrated flux [erg/cm2/s/A]'
plt.ylabel(yLabel)
plt.xlim(4280.,7095.)
plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN+GTC+Hiltner102+fibre78_cal.png'
plt.savefig(plotname)
plt.show()
