execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

calibrated = True
useMean = True

gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_sum_n.fits"
wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_sum_n.fits"
wiyn_file = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15T08-22-49_botzfxsEcBld_rebinned_52_n.fits'
if calibrated:
    gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned.fits"
    wiyn_file = "/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned.fits"#Pa30_WN151014_cal_sum_cleaned.fits"
    leDu_file = "/Users/azuri/daten/uni/HKU/Pa30/LeDu/_pa30_20181009_941_PLeDu_scaled.fits"
    somme_file = "/Users/azuri/daten/uni/HKU/Pa30/_pa30_somme6_scaled.fits"
    if not useMean:
        wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMedian_52_cal_cleaned.fits"#Pa30_WN151014_cal_sum_cleaned.fits"


gtc_hdulist = pyfits.open(gtc_file)
wiyn_hdulist = pyfits.open(wiyn_file)
leDu_hdulist = pyfits.open(leDu_file)
somme_hdulist = pyfits.open(somme_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum = fits.getdata(gtc_file)
print('gtc_wavelength = ',gtc_wavelength)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1)
wiyn_spectrum = fits.getdata(wiyn_file)
wiyn_spectrum_absolute = calibratedFluxToAbsoluteFlux(wiyn_spectrum, 3370.0)
print('wiyn_hdulist[0].data = ',wiyn_hdulist[0].data)
wiyn_hdulist[0].data = wiyn_spectrum_absolute
specOutName = wiyn_file[0:wiyn_file.rfind('.')]+'_absoluteFlux.fits'
wiyn_hdulist.writeto(specOutName, clobber=True)
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

xLim = [3675,7500]
xLimStr = [str(x) for x in xLim]
yLim = [0.,4.9e-15]
yLimStr = [str(y) for y in yLim]
print('yLimStr = ',yLimStr)
yLimStr[0] = yLimStr[0][0:yLimStr[0].find('.')]
yLimStr[1] = yLimStr[1].replace('.','_')


#plt.plot(leDu_wavelength,leDu_spectrum,'g-', label='PLeDu scaled and smoothed')
#plt.plot(somme_wavelength,somme_spectrum,'c-', label='Somme scaled and smoothed')
plt.plot(gtc_wavelength,gtc_spectrum,'r-', label='GTC')
plt.plot(wiyn_wavelength[5:len(wiyn_wavelength)-5],wiyn_spectrum[5:len(wiyn_wavelength)-5],'b-', label='WIYN')
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'counts [ADU]'
if calibrated:
    yLabel = 'calibrated flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])
plt.ylim(yLim[0],yLim[1])
#plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC_'
if calibrated:
    plotname = plotname + 'cal_'
plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()


#gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_rebinned_4000px_cleaned.fits"#Pa30_GT080716_cal_sum_cleaned.fits"
#wiyn_file = "/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned_scaled_with_constant.fits"#Pa30_WN151014_cal_sum_rebinned_scaled_down_cleaned.fits"
#gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_rebinned_4281-7084_scaled_order2.fits"#Pa30_GT080716_cal_sum_cleaned_scaled.fits"
gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned_scaled.fits"

gtc_hdulist = pyfits.open(gtc_file)
#wiyn_hdulist = pyfits.open(wiyn_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum = fits.getdata(gtc_file)
gtc_spectrum_absolute = calibratedFluxToAbsoluteFlux(gtc_spectrum, 3370.0)
gtc_hdulist[0].data = gtc_spectrum_absolute
specOutName = gtc_file[0:gtc_file.rfind('.')]+'_absoluteFlux.fits'
gtc_hdulist.writeto(specOutName, clobber=True)

print('gtc_wavelength = ',gtc_wavelength)




plt.plot(gtc_wavelength,gtc_spectrum_absolute,'r-', label='GTC scaled')
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'counts [ADU]'
if calibrated:
    yLabel = 'absolute flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_GTC_absoluteFlux_'
plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()



plt.plot(gtc_wavelength,np.log10(gtc_spectrum_absolute),'r-', label='GTC scaled')
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'counts [ADU]'
if calibrated:
    yLabel = 'log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_GTC_logAbsoluteFlux_'
plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()





plt.plot(leDu_wavelength,leDu_spectrum_absolute_smoothed,'g-', label='LeDu scaled and smoothed')
#plt.plot(somme_wavelength,somme_spectrum_absolute_smoothed,'c-', label='Somme scaled and smoothed')
plt.plot(wiyn_wavelength,wiyn_spectrum_absolute,'b-', label='WIYN')
plt.plot(gtc_wavelength,gtc_spectrum_absolute,'r-', label='GTC scaled')
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'absolute flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_GTC_WIYN_LeDu-smoothed_absoluteFlux_'
plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()



plt.plot(leDu_wavelength,np.log10(leDu_spectrum_absolute_smoothed),'g-', label='LeDu scaled and smoothed')
#plt.plot(somme_wavelength,np.log10(somme_spectrum_absolute_smoothed),'c-', label='Somme scaled and smoothed')
plt.plot(wiyn_wavelength,np.log10(wiyn_spectrum_absolute),'b-', label='WIYN')
plt.plot(gtc_wavelength,np.log10(gtc_spectrum_absolute),'r-', label='GTC scaled')
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'counts [ADU]'
if calibrated:
    yLabel = 'log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_GTC_WIYN_LeDu-smoothed_logAbsoluteFlux_'
plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()




plt.plot(leDu_wavelength,leDu_spectrum_absolute_smoothed,'g-', label='LeDu scaled and smoothed')
plt.plot(somme_wavelength,somme_spectrum_absolute_smoothed,'c-', label='Somme scaled and smoothed')
plt.plot(wiyn_wavelength,wiyn_spectrum_absolute,'b-', label='WIYN')
plt.plot(gtc_wavelength,gtc_spectrum_absolute,'r-', label='GTC scaled')
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'absolute flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_GTC_WIYN_LeDu-smoothed_Somme-smoothed_absoluteFlux_'
plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()



plt.plot(leDu_wavelength,np.log10(leDu_spectrum_absolute_smoothed),'g-', label='LeDu scaled and smoothed')
plt.plot(somme_wavelength,np.log10(somme_spectrum_absolute_smoothed),'c-', label='Somme scaled and smoothed')
plt.plot(wiyn_wavelength,np.log10(wiyn_spectrum_absolute),'b-', label='WIYN')
plt.plot(gtc_wavelength,np.log10(gtc_spectrum_absolute),'r-', label='GTC scaled')
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'counts [ADU]'
if calibrated:
    yLabel = 'log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_GTC_WIYN_LeDu-smoothed_Somme-smoothed_logAbsoluteFlux_'
plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()






#wiyn_header = wiyn_hdulist[0].header
#wiyn_wavelength = getWavelength(wiyn_header,1)
#wiyn_spectrum = fits.getdata(wiyn_file)
#print('wiyn_wavelength = ',wiyn_wavelength)

plt.plot(leDu_wavelength,leDu_spectrum,'g-', label='LeDu scaled and smoothed')
#plt.plot(somme_wavelength,somme_spectrum,'c-', label='Pa30 Somme scaled')
plt.plot(gtc_wavelength,gtc_spectrum,'r-', label='GTC scaled')
plt.plot(wiyn_wavelength[2:],wiyn_spectrum[2:],'b-', label='WIYN')

xLim = [3700,7500]
xLimStr = [str(x) for x in xLim]
yLim = [0,1.2e-14]
yLimStr = [str(x) for x in yLim]
yLimStr[0] = yLimStr[0][0:yLimStr[0].find('.')]
yLimStr[1] = yLimStr[1].replace('.','_')

plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'calibrated flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])#4280.,7095.)
#plt.xlim(3600.,4100.)#7900.)#4280.,7095.)
#plt.ylim(0.0, 1.5e-15)
#plt.ylim(yLim[0], yLim[1])#1.5e-15)
#plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+LeDu_scaled_with_constant_x_fit_order1_'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
#plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+LeDu_scaled_with_constant_x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

#plt.plot(leDu_wavelength,leDu_spectrum,'g-', label='LeDu scaled and smoothed')
#plt.plot(somme_wavelength,somme_spectrum,'c-', label='Pa30 Somme scaled')
plt.plot(gtc_wavelength,gtc_spectrum,'r-', label='GTC scaled')
plt.plot(wiyn_wavelength[2:],wiyn_spectrum[2:],'b-', label='WIYN')

plt.plot([4245.,4430.],[1.25e-15,1.25e-15],'k')
plt.plot([4245.,4245.],[1.22e-15,1.28e-15],'k')
plt.plot([4337.5,4337.5],[1.22e-15,1.28e-15],'k')
plt.plot([4430.,4430.],[1.22e-15,1.28e-15],'k')#12800 km/s 4.3%

plt.plot([4338,4522],[1.1e-15,1.1e-15],'m')
plt.plot([4338,4338],[1.07e-15,1.13e-15],'m')
plt.plot([4430.,4430.],[1.07e-15,1.13e-15],'m')
plt.plot([4522.,4522.],[1.07e-15,1.13e-15],'m')#12500 km/s 4.2%

plt.plot([4420.,4590.],[1.2e-15,1.2e-15],'k')
plt.plot([4420.,4420.],[1.17e-15,1.23e-15],'k')
plt.plot([4505.,4505.],[1.17e-15,1.23e-15],'k')
plt.plot([4590.,4590.],[1.17e-15,1.23e-15],'k')#11300 km/s 3.8%

plt.plot([4505.,4676.],[1.15e-15,1.15e-15],'m')
plt.plot([4505.,4505.],[1.12e-15,1.18e-15],'m')
plt.plot([4590.,4590.],[1.12e-15,1.18e-15],'m')#11200 km/s 3.7%
plt.plot([4676.,4676.],[1.12e-15,1.18e-15],'m')

plt.plot([4590.,4734.],[1.25e-15,1.25e-15],'k')
plt.plot([4590.,4590.],[1.22e-15,1.28e-15],'k')
plt.plot([4662.,4662.],[1.22e-15,1.28e-15],'k')#9300 km/s 3.1%
plt.plot([4734.,4734.],[1.22e-15,1.28e-15],'k')

plt.plot([5140.,5420.],[1.1e-15,1.1e-15],'k')
plt.plot([5140.,5140.],[1.07e-15,1.13e-15],'k')
plt.plot([5280.,5280.],[1.07e-15,1.13e-15],'k')#15900 km/s 5.3%
plt.plot([5420.,5420.],[1.07e-15,1.13e-15],'k')

plt.plot([5070.,5490.],[1.2e-15,1.2e-15],'k')
plt.plot([5070.,5070.],[1.17e-15,1.23e-15],'k')
plt.plot([5280.,5280.],[1.17e-15,1.23e-15],'k')#23800 km/s 8%
plt.plot([5490.,5490.],[1.17e-15,1.23e-15],'k')

xLim = [4000,5600]
xLimStr = [str(x) for x in xLim]
yLim = [1.0e-15,3.2e-15]
yLimStr = [str(x) for x in yLim]
yLimStr[0] = yLimStr[0][0:yLimStr[0].find('.')]
yLimStr[1] = yLimStr[1].replace('.','_')

plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'calibrated flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plt.xlim(xLim[0],xLim[1])#4280.,7095.)
#plt.xlim(3600.,4100.)#7900.)#4280.,7095.)
#plt.ylim(0.0, 1.5e-15)
plt.ylim(yLim[0], yLim[1])#1.5e-15)
#plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+LeDu_scaled_with_constant_x_fit_order1_'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
#plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+LeDu_scaled_with_constant_x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

ratio_file = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_WIYNdivbyGTC.fits'

ratio_hdulist = pyfits.open(ratio_file)
#wiyn_hdulist = pyfits.open(wiyn_file)

ratio_header = ratio_hdulist[0].header
ratio_wavelength = getWavelength(ratio_header,1)
ratio_spectrum = fits.getdata(ratio_file)

plt.plot(ratio_wavelength, ratio_spectrum)
plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'calibrated flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WIYNdivbyGTC.pdf'
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()





plt.plot(gtc_wavelength,gtc_spectrum,'r-', label='GTC scaled')
#plt.plot(wiyn_wavelength[2:],wiyn_spectrum[2:],'b-', label='WIYN')
plt.plot(leDu_wavelength,leDu_spectrum,'g-', label='LeDu scaled and smoothed')

xLim = [5800,6260]
xLimStr = [str(x) for x in xLim]
yLim = [0,1.2e-14]
yLimStr = [str(x) for x in yLim]
yLimStr[0] = yLimStr[0][0:yLimStr[0].find('.')]
yLimStr[1] = yLimStr[1].replace('.','_')

plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'calibrated flux [$\mathrm{erg/cm^2/s/\AA}$]'
plt.ylabel(yLabel)
#plt.xlim(xLim[0],xLim[1])#4280.,7095.)
#plt.xlim(3600.,4100.)#7900.)#4280.,7095.)
#plt.ylim(0.0, 1.5e-15)
#plt.ylim(yLim[0], yLim[1])#1.5e-15)
#plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_GTC+LeDu_scaled_smoothed_'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
#plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+LeDu_scaled_with_constant_x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()


plt.plot(somme_wavelength,np.log10(somme_spectrum_absolute_smoothed),'c-', label='Somme scaled and smoothed')
plt.plot(leDu_wavelength,np.log10(leDu_spectrum_absolute_smoothed),'g-', label='LeDu scaled and smoothed')
plt.plot(gtc_wavelength,np.log10(gtc_spectrum_absolute),'r-', label='GTC scaled')
plt.plot(wiyn_wavelength[2:],np.log10(wiyn_spectrum_absolute[2:]),'b-', label='WIYN')

xLim = [5800,6260]
xLimStr = [str(x) for x in xLim]
yLim = [0,1.2e-14]
yLimStr = [str(x) for x in yLim]
yLimStr[0] = yLimStr[0][0:yLimStr[0].find('.')]
yLimStr[1] = yLimStr[1].replace('.','_')

plt.xlabel('wavelength [$\mathrm{\AA}$]')
yLabel = 'log(absolute flux / [$\mathrm{erg/cm^2/s/\AA}$])'
plt.ylabel(yLabel)
#plt.xlim(xLim[0],xLim[1])#4280.,7095.)
#plt.xlim(3600.,4100.)#7900.)#4280.,7095.)
#plt.ylim(0.0, 1.5e-15)
#plt.ylim(yLim[0], yLim[1])#1.5e-15)
#plt.title('Pa30')
plt.legend()
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+Somme+LeDu_absolute_scaled_smoothed_'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
#plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+LeDu_scaled_with_constant_x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='pdf', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()
