from astropy.io import fits
import astropy.io.fits as pyfits
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
import numpy as np

from myUtils import getWavelength,calibratedFluxToAbsoluteFlux,boxCarMeanSmooth
#execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

#ebv = 1.206
ebv = 0.79 # 3D map Bayestar15
#ebv = 0.99 # 3D map Bayestar17
R_V = 3.1  # (Fitzpatrick and Massa 2007)
#R_V = 3.1
A_V = 2.4 #Claire
ebv = A_V / R_V
print('A_V = ',A_V,' => ebv = ',ebv)

if False:
    gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned_scaled_absoluteFlux.fits"
    wiyn_file = "/Volumes/work/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned_absoluteFlux.fits"#Pa30_WN151014_cal_sum_cleaned.fits"
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

    leDu_spectrum_absolute_smoothed_dereddened = pyasl.unred(leDu_wavelength, leDu_spectrum_absolute_smoothed, ebv=ebv, R_V = R_V)
    plt.plot(leDu_wavelength,np.log10(leDu_spectrum_absolute_smoothed),'-', color='#64a030', label='LeDu scaled and smoothed')
    plt.plot(leDu_wavelength,np.log10(leDu_spectrum_absolute_smoothed_dereddened),'-', color='#64a030')#, label='Le Du scaled, smoothed, and dereddened')

    somme_spectrum_absolute_smoothed_dereddened = pyasl.unred(somme_wavelength, somme_spectrum_absolute_smoothed, ebv=ebv, R_V = R_V)
    plt.plot(somme_wavelength,np.log10(somme_spectrum_absolute_smoothed),'m-', label='Somme scaled and smoothed')
    plt.plot(somme_wavelength,np.log10(somme_spectrum_absolute_smoothed_dereddened),'m-')#, label='Somme scaled, smoothed, and dereddened')

    gtc_spectrum_absolute_dereddened = pyasl.unred(gtc_wavelength, gtc_spectrum_absolute, ebv=ebv, R_V = R_V)
    plt.plot(gtc_wavelength, np.log10(gtc_spectrum_absolute), 'r-', label='GTC scaled')
    plt.plot(gtc_wavelength, np.log10(gtc_spectrum_absolute_dereddened), 'r-')#, label='GTC scaled')
    #plt.xlabel('wavelength [\AA]')
    #plt.ylabel('log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])')
    #plt.legend()
    #plt.show()

    wiyn_spectrum_absolute_dereddened = pyasl.unred(wiyn_wavelength, wiyn_spectrum_absolute, ebv=ebv, R_V = R_V)
    plt.plot(wiyn_wavelength, np.log10(wiyn_spectrum_absolute), 'b-', label='WIYN')
    plt.plot(wiyn_wavelength, np.log10(wiyn_spectrum_absolute_dereddened), 'b-')#, label='WIYN dereddened')
    plt.xlabel("wavelength [$\mathrm{\AA}$]")
    plt.ylabel('log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])')
    plt.legend()
    xlim = [3550.,7585.]
    plt.xlim(xlim)
    #plt.ylim([-10.,-5.9])
    plt.ylim([-10.0,-7.0])
    plt.text(4123., -9.265, 'original')
    #plt.text(4123., -7.0846, 'dereddened with E(B-V)=%.2f'%(ebv))
    plt.text(5550., -8.4, 'dereddened with E(B-V)=%.2f'%(ebv))
    #plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC_absolute_dereddened_R_V=%.1f_ebv=%.3f_x%d-%d' % (R_V, ebv, int(xlim[0]), int(xlim[1]))
    plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+Somme+LeDu_absolute_dereddened_R_V=%.1f_ebv=%.3f_x%d-%d' % (R_V, ebv, int(xlim[0]), int(xlim[1]))
    plotname = plotname.replace('.','_')+'.eps'
    #plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
    print('writing plot to file <'+plotname+'>')
    plt.savefig(plotname, format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
    plt.show()


# observed Flux

gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned_scaled.fits"
#wiyn_file = "/Volumes/work/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned.fits"#Pa30_WN151014_cal_sum_cleaned.fits"

gtc_hdulist = pyfits.open(gtc_file)
#wiyn_hdulist = pyfits.open(wiyn_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum = fits.getdata(gtc_file)
gtc_spectrum_dereddened = pyasl.unred(gtc_wavelength, gtc_spectrum, ebv=ebv, R_V = R_V)
print('gtc_wavelength = ',gtc_wavelength)
print('len(gtc_hdulist) = ',len(gtc_hdulist))
gtc_hdulist[0].data = gtc_spectrum_dereddened
gtc_hdulist.writeto(gtc_file[0:gtc_file.rfind('.')]+'_dereddened.fits', overwrite=False)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1)
wiyn_spectrum = fits.getdata(wiyn_file)
print('wiyn_wavelength = ',wiyn_wavelength)

leDu_spectrum_smoothed_dereddened = pyasl.unred(leDu_wavelength, leDu_spectrum_smoothed, ebv=ebv, R_V = R_V)
plt.plot(leDu_wavelength,np.log10(leDu_spectrum_smoothed),'-', color='#64a030', label='LeDu scaled and smoothed')
plt.plot(leDu_wavelength,np.log10(leDu_spectrum_smoothed_dereddened),'-', color='#64a030')#, label='Le Du scaled, smoothed, and dereddened')

somme_spectrum_smoothed_dereddened = pyasl.unred(somme_wavelength, somme_spectrum_smoothed, ebv=ebv, R_V = R_V)
plt.plot(somme_wavelength,np.log10(somme_spectrum_smoothed),'m-', label='Somme scaled and smoothed')
plt.plot(somme_wavelength,np.log10(somme_spectrum_smoothed_dereddened),'m-')#, label='Somme scaled, smoothed, and dereddened')

plt.plot(gtc_wavelength, np.log10(gtc_spectrum), 'r-', label='GTC scaled')
plt.plot(gtc_wavelength, np.log10(gtc_spectrum_dereddened), 'r-')#, label='GTC scaled')
#plt.xlabel('wavelength [\AA]')
#plt.ylabel('log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])')
#plt.legend()
#plt.show()

wiyn_spectrum_dereddened = pyasl.unred(wiyn_wavelength, wiyn_spectrum, ebv=ebv, R_V = R_V)
plt.plot(wiyn_wavelength, np.log10(wiyn_spectrum), 'b-', label='WIYN')
plt.plot(wiyn_wavelength, np.log10(wiyn_spectrum_dereddened), 'b-')#, label='WIYN dereddened')
plt.xlabel("wavelength [$\mathrm{\AA}$]")
plt.ylabel('log(flux [$\mathrm{erg/cm^2/s/\AA}$])')
plt.legend()
xlim = [3550.,7585.]
plt.xlim(xlim)
#plt.ylim([-10.,-5.9])
plt.ylim([-15.1,-12.2])
plt.text(4123., -14.4, 'original')
#plt.text(4123., -7.0846, 'dereddened with E(B-V)=%.2f'%(ebv))
plt.text(5550., -13.56, 'dereddened with E(B-V)=%.2f'%(ebv))
#plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC_dereddened_R_V=%.1f_ebv=%.3f_x%d-%d' % (R_V, ebv, int(xlim[0]), int(xlim[1]))
plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC+Somme+LeDu_dereddened_R_V=%.1f_ebv=%.3f_x%d-%d' % (R_V, ebv, int(xlim[0]), int(xlim[1]))
plotname = plotname.replace('.','_')+'.eps'
#plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()



