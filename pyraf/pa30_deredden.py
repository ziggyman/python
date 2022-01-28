from astropy.io import fits
import astropy.io.fits as pyfits
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from Pa30_LBT import readLBTFiles,readGvaramadzeFile

from myUtils import getWavelength,calibratedFluxToAbsoluteFlux,boxCarMeanSmooth,smooth
#execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

#ebv = 1.206
ebv = 0.79 # 3D map Bayestar15
#ebv = 0.99 # 3D map Bayestar17
R_V = 3.1  # (Fitzpatrick and Massa 2007)
#R_V = 3.1
A_V = 2.4 #Claire
ebv = A_V / R_V
print('A_V = ',A_V,' => ebv = ',ebv)

lineRangesScalePlot = [[[4281.1,5454.],[4200.,6900.]],
                       [[3704.6,4053.3],[3600.,4053.3]],
                       [[4282.,4430.],[4230.,4430.]],
                       [[4429.,4590.],[4429.,4590.]],
                       [[4282.,4590.],[4200.,4600.]],
                       [[4590.,4736.],[4590.,4736.]],
                       [[4731.,4874.],[4731.,4874.]],
                       [[5112.,5455.],[5100.,5700.]],#5492.
                       [[5912.,6268.],[5800.,6300.]],#5912.,6268.]],
                      ]

leDu_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/_pa30_20181009_941_PLeDu_scaled.fits"
somme_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/_pa30_somme6_scaled.fits"
gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_GT080716_cal_sum_cleaned_scaled.fits"
wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_WN151014_cal_sum_cleaned_t.fits"
gvaramadze_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/gvaramadze.fits"
#    wiyn_file = "/Volumes/work/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned_absoluteFlux.fits"#Pa30_WN151014_cal_sum_cleaned.fits"

gtc_hdulist = pyfits.open(gtc_file)
wiyn_hdulist = pyfits.open(wiyn_file)
leDu_hdulist = pyfits.open(leDu_file)
somme_hdulist = pyfits.open(somme_file)
gvaramadze_hdulist = pyfits.open(gvaramadze_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum = fits.getdata(gtc_file)
print('gtc_wavelength = ',gtc_wavelength)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1)
wiyn_spectrum = fits.getdata(wiyn_file)
print('wiyn_wavelength = ',wiyn_wavelength)

leDu_header = leDu_hdulist[0].header
leDu_wavelength = getWavelength(leDu_header,1)
print('leDu_wavelength = ',leDu_wavelength)
leDu_spectrum = fits.getdata(leDu_file)
#    leDu_spectrum_absolute = calibratedFluxToAbsoluteFlux(leDu_spectrum, 3370.0)
print('leDu_wavelength = ',leDu_wavelength)

somme_header = somme_hdulist[0].header
somme_wavelength = getWavelength(somme_header,1)
somme_spectrum = fits.getdata(somme_file)
#    somme_spectrum_absolute = calibratedFluxToAbsoluteFlux(somme_spectrum, 3370.0)
print('somme_wavelength = ',somme_wavelength)

leDu_spectrum_smoothed = scipy.ndimage.median_filter(leDu_spectrum, 21)#boxCarMeanSmooth(leDu_spectrum, 0, 21)
somme_spectrum_smoothed = scipy.ndimage.median_filter(somme_spectrum, 21)#boxCarMeanSmooth(somme_spectrum, 0, 21)
#leDu_spectrum_absolute_smoothed = scipy.ndimage.median_filter(leDu_spectrumboxCarMeanSmooth(leDu_spectrum_absolute, 0, 21)
#somme_spectrum_absolute_smoothed = boxCarMeanSmooth(somme_spectrum_absolute, 0, 21)

leDu_spectrum_smoothed_dereddened = pyasl.unred(leDu_wavelength, leDu_spectrum_smoothed, ebv=ebv, R_V = R_V)
plt.plot(leDu_wavelength,np.log10(leDu_spectrum_smoothed),'-', color='#64a030', label='LeDu scaled and smoothed')
plt.plot(leDu_wavelength,np.log10(leDu_spectrum_smoothed_dereddened),'-', color='#64a030')#, label='Le Du scaled, smoothed, and dereddened')

somme_spectrum_smoothed_dereddened = pyasl.unred(somme_wavelength, somme_spectrum_smoothed, ebv=ebv, R_V = R_V)
plt.plot(somme_wavelength,np.log10(somme_spectrum_smoothed),'m-', label='Somme scaled and smoothed')
plt.plot(somme_wavelength,np.log10(somme_spectrum_smoothed_dereddened),'m-')#, label='Somme scaled, smoothed, and dereddened')

gtc_spectrum_dereddened = pyasl.unred(gtc_wavelength, gtc_spectrum, ebv=ebv, R_V = R_V)
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
#    plt.xlim(xlim)
#plt.ylim([-10.,-5.9])
#    plt.ylim([-10.0,-7.0])
plt.text(4123., -9.265, 'original')
#plt.text(4123., -7.0846, 'dereddened with E(B-V)=%.2f'%(ebv))
plt.text(5550., -8.4, 'dereddened with E(B-V)=%.2f'%(ebv))
#plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC_dereddened_R_V=%.1f_ebv=%.3f_x%d-%d' % (R_V, ebv, int(xlim[0]), int(xlim[1]))
plotname = '/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_WN+GTC+Somme+LeDu_dereddened_logarithmic_R_V=%.1f_ebv=%.3f_x%d-%d' % (R_V, ebv, int(xlim[0]), int(xlim[1]))
plotname = plotname.replace('.','_')+'.eps'
#plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()


# observed Flux
gvaramadze_header = gvaramadze_hdulist[0].header
gvaramadze_wavelength = getWavelength(gvaramadze_header,1)
gvaramadze_spectrum = fits.getdata(gvaramadze_file)
#gvaramadze_wavelength,gvaramadze_spectrum = readGvaramadzeFile()
gvaramadze_wavelength = np.array(gvaramadze_wavelength)
gvaramadze_spectrum = np.array(gvaramadze_spectrum)
garnavich_wavelength, garnavich_spectrum = readLBTFiles()
garnavich_wavelength = np.array(garnavich_wavelength)
garnavich_spectrum = np.array(garnavich_spectrum)
garnavich_spectrum = garnavich_spectrum[garnavich_wavelength < 5455.]
garnavich_wavelength = garnavich_wavelength[garnavich_wavelength < 5455.]
garnavich_spectrum_smoothed = smooth(garnavich_spectrum,9)#scipy.ndimage.mean_filter(garnavich_spectrum, 7)#boxCarMeanSmooth(somme_spectrum, 0, 21)

#gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_GT080716_cal_sum_cleaned_scaled.fits"
#wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_WN151014_cal_sum_cleaned_t.fits"
#wiyn_file = "/Volumes/work/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned.fits"#Pa30_WN151014_cal_sum_cleaned.fits"

gtc_hdulist = pyfits.open(gtc_file)
wiyn_hdulist = pyfits.open(wiyn_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1)
gtc_spectrum = fits.getdata(gtc_file)
gtc_spectrum_dereddened = pyasl.unred(gtc_wavelength, gtc_spectrum, ebv=ebv, R_V = R_V)
print('gtc_wavelength = ',gtc_wavelength)
print('len(gtc_hdulist) = ',len(gtc_hdulist))
gtc_hdulist[0].data = gtc_spectrum_dereddened
gtc_hdulist.writeto(gtc_file[0:gtc_file.rfind('.')]+'_dereddened.fits', overwrite=True)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1)
wiyn_spectrum = fits.getdata(wiyn_file)
wiyn_spectrum = smooth(wiyn_spectrum,5)#scipy.ndimage.mean_filter(wiyn_spectrum, 7)
print('wiyn_wavelength = ',wiyn_wavelength)

leDu_hdulist = pyfits.open(leDu_file)
somme_hdulist = pyfits.open(somme_file)
leDu_header = leDu_hdulist[0].header
leDu_wavelength = getWavelength(leDu_header,1)
print('leDu_wavelength = ',leDu_wavelength)
leDu_spectrum = fits.getdata(leDu_file)
#    leDu_spectrum_absolute = calibratedFluxToAbsoluteFlux(leDu_spectrum, 3370.0)
print('leDu_wavelength = ',leDu_wavelength)

scaleRangeGvaramadze = [4300.,6800.]
meanGvaramadzeSpectrum = np.mean(gvaramadze_spectrum[(gvaramadze_wavelength > scaleRangeGvaramadze[0]) & (gvaramadze_wavelength < scaleRangeGvaramadze[1])])
print('meanGvaramadzeSpectrum = ',meanGvaramadzeSpectrum)
meanWiynSpectrum = np.mean(wiyn_spectrum[(wiyn_wavelength > scaleRangeGvaramadze[0]) & (wiyn_wavelength < scaleRangeGvaramadze[1])])
print('meanWiynSpectrum = ',meanWiynSpectrum)
gvaramadze_spectrum_scaled = gvaramadze_spectrum * meanWiynSpectrum / meanGvaramadzeSpectrum

scaleRangeGarnavich = [4300.,5500.]
meanGarnavichSpectrum = np.mean(garnavich_spectrum_smoothed[(garnavich_wavelength > scaleRangeGarnavich[0]) & (garnavich_wavelength < scaleRangeGarnavich[1])])
print('meanGarnavichSpectrum = ',meanGarnavichSpectrum)
meanWiynSpectrum = np.mean(wiyn_spectrum[(wiyn_wavelength > scaleRangeGarnavich[0]) & (wiyn_wavelength < scaleRangeGarnavich[1])])
print('meanWiynSpectrum = ',meanWiynSpectrum)
Garnavich_spectrum_smoothed_scaled = garnavich_spectrum_smoothed * meanWiynSpectrum / meanGarnavichSpectrum


#plt.show()
#STOP

leDu_spectrum_smoothed_dereddened = pyasl.unred(leDu_wavelength, leDu_spectrum_smoothed, ebv=ebv, R_V = R_V)
#plt.plot(leDu_wavelength,leDu_spectrum_smoothed,'-', color='#64a030', label='LeDu scaled and smoothed')
#plt.plot(leDu_wavelength,leDu_spectrum_smoothed_dereddened,'-', color='#64a030')#, label='Le Du scaled, smoothed, and dereddened')

somme_spectrum_smoothed_dereddened = pyasl.unred(somme_wavelength, somme_spectrum_smoothed, ebv=ebv, R_V = R_V)
#plt.plot(somme_wavelength,somme_spectrum_smoothed,'m-', label='Somme scaled and smoothed')
#plt.plot(somme_wavelength,somme_spectrum_smoothed_dereddened,'m-')#, label='Somme scaled, smoothed, and dereddened')

#plt.plot(gtc_wavelength, gtc_spectrum_dereddened, 'r-')#, label='GTC scaled')
#plt.xlabel('wavelength [\AA]')
#plt.ylabel('log(absolute flux [$\mathrm{erg/cm^2/s/\AA}$])')
#plt.legend()
#plt.show()

wiyn_spectrum_dereddened = pyasl.unred(wiyn_wavelength, wiyn_spectrum, ebv=ebv, R_V = R_V)
wiyn_hdulist[0].data = wiyn_spectrum_dereddened
wiyn_hdulist.writeto(wiyn_file[0:wiyn_file.rfind('.')]+'_dereddened.fits', overwrite=True)
plt.rcParams.update({'font.size': 14})
plt.plot(garnavich_wavelength,Garnavich_spectrum_smoothed_scaled,'-', color='#64a030', label='Garnavich 11/09/2020')
plt.plot(gvaramadze_wavelength,gvaramadze_spectrum_scaled,'m-', label='Gvaramadze 20/07/2017')
plt.plot(gtc_wavelength, gtc_spectrum, 'r-', label='GTC 08/07/2016')
plt.plot(wiyn_wavelength, wiyn_spectrum, 'b-', label='WIYN 15/10/2014')
#plt.plot(wiyn_wavelength, wiyn_spectrum_dereddened, 'b-')#, label='WIYN dereddened')
plt.xlabel("wavelength [$\mathrm{\AA}$]")
plt.ylabel('flux [$\mathrm{erg/cm^2/s/\AA}$]')
plt.legend(fontsize=12)
xlim = [3550.,7585.]
#plt.xlim(xlim)
#plt.ylim([-10.,-5.9])
#plt.ylim([-15.1,-12.2])
#plt.text(4123., -14.4, 'original')
#plt.text(4123., -7.0846, 'dereddened with E(B-V)=%.2f'%(ebv))
#plt.text(5550., -13.56, 'dereddened with E(B-V)=%.2f'%(ebv))
#plotname = '/Users/azuri/daten/uni/HKU/Pa30/report/images/Pa30_WN+GTC_dereddened_R_V=%.1f_ebv=%.3f_x%d-%d' % (R_V, ebv, int(xlim[0]), int(xlim[1]))
plotname = '/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_WN+GTC+Garnavich+Gvaramadze.fits'#_dereddened_R_V=%.1f_ebv=%.3f_x%d-%d' % (R_V, ebv, int(xlim[0]), int(xlim[1]))
plotname = plotname.replace('.','_')+'.eps'
#plotname = plotname + 'x'+xLimStr[0]+'-'+xLimStr[1]+'_'+yLimStr[0]+'-'+yLimStr[1]+'.pdf'
print('writing plot to file <'+plotname+'>')
plt.savefig(plotname, format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
plt.show()

for line in lineRangesScalePlot:
    scaleRange = line[0]
    plotRange = line[1]
    print('scaleRange = ',scaleRange,', plotRange = ',plotRange)
    plotWiyn = True
    if (wiyn_wavelength[0] > scaleRange[0]) or (wiyn_wavelength[len(wiyn_wavelength)-1] < scaleRange[1]):
        plotWiyn = False
        print('wiyn_wavelength[0] = ',wiyn_wavelength[0],', wiyn_wavelength[',len(wiyn_wavelength)-1,'] = ',wiyn_wavelength[len(wiyn_wavelength)-1],' => plotWiyn set to False')
    plotGTC = True
    if (gtc_wavelength[0] > scaleRange[0]) or (gtc_wavelength[len(gtc_wavelength)-1] < scaleRange[1]):
        plotGTC = False
        print('gtc_wavelength[0] = ',gtc_wavelength[0],', gtc_wavelength[',len(gtc_wavelength)-1,'] = ',gtc_wavelength[len(gtc_wavelength)-1],' => plotGTC set to False')
    plotGvaramadze = True
    if (gvaramadze_wavelength[0] > scaleRange[0]) or (gvaramadze_wavelength[len(gvaramadze_wavelength)-1] < scaleRange[1]):
        plotGvaramadze = False
        print('gvaramadze_wavelength[0] = ',gvaramadze_wavelength[0],', gvaramadze_wavelength[',len(gvaramadze_wavelength)-1,'] = ',gvaramadze_wavelength[len(gvaramadze_wavelength)-1],' => plotGvaramadze set to False')
    plotGarnavich = True
    if (garnavich_wavelength[0] > scaleRange[0]) or (garnavich_wavelength[len(garnavich_wavelength)-1] < scaleRange[1]):
        plotGarnavich = False
        print('garnavich_wavelength[0] = ',garnavich_wavelength[0],', garnavich_wavelength[',len(garnavich_wavelength)-1,'] = ',garnavich_wavelength[len(garnavich_wavelength)-1],' => plotGarnavich set to False')

    if plotGarnavich:
        print('garnavich_wavelength = ',garnavich_wavelength)
        print('plotRange = ',plotRange)
        garnavich_spectrum_line = garnavich_spectrum_smoothed[(garnavich_wavelength >= plotRange[0]) & (garnavich_wavelength <= plotRange[1])] / np.mean(garnavich_spectrum_smoothed[(garnavich_wavelength >= scaleRange[0]) & (garnavich_wavelength <= scaleRange[1])])
        garnavich_wavelength_line = garnavich_wavelength[(garnavich_wavelength >= plotRange[0]) & (garnavich_wavelength <= plotRange[1])]
        plt.plot(garnavich_wavelength_line,garnavich_spectrum_line,'-', color='#64a030', label='Garnavich 11/09/2020')
    if plotGvaramadze:
        gvaramadze_spectrum_line = gvaramadze_spectrum[(gvaramadze_wavelength >= plotRange[0]) & (gvaramadze_wavelength <= plotRange[1])] / np.mean(gvaramadze_spectrum[(gvaramadze_wavelength >= scaleRange[0]) & (gvaramadze_wavelength <= scaleRange[1])])
        gvaramadze_wavelength_line = gvaramadze_wavelength[(gvaramadze_wavelength >= plotRange[0]) & (gvaramadze_wavelength <= plotRange[1])]
        plt.plot(gvaramadze_wavelength_line,gvaramadze_spectrum_line,'m-', label='Gvaramadze 20/07/2017')
    if plotGTC:
        gtc_spectrum_line = gtc_spectrum[(gtc_wavelength >= plotRange[0]) & (gtc_wavelength <= plotRange[1])] / np.mean(gtc_spectrum[(gtc_wavelength >= scaleRange[0]) & (gtc_wavelength <= scaleRange[1])])
        gtc_wavelength_line = gtc_wavelength[(gtc_wavelength >= plotRange[0]) & (gtc_wavelength <= plotRange[1])]
        plt.plot(gtc_wavelength_line,gtc_spectrum_line,'r-', label='GTC 08/07/2016')
    if plotWiyn:
        wiyn_spectrum_line = wiyn_spectrum[(wiyn_wavelength >= plotRange[0]) & (wiyn_wavelength <= plotRange[1])] / np.mean(wiyn_spectrum[(wiyn_wavelength >= scaleRange[0]) & (wiyn_wavelength <= scaleRange[1])])
        wiyn_wavelength_line = wiyn_wavelength[(wiyn_wavelength >= plotRange[0]) & (wiyn_wavelength <= plotRange[1])]
        plt.plot(wiyn_wavelength_line,wiyn_spectrum_line,'b-', label='WIYN 15/10/2014')
    plt.xlabel("wavelength [$\mathrm{\AA}$]")
    plt.ylabel('flux [arbitrary units]')
    plt.legend(fontsize=12)
    plotname = '/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_WN+GTC+Garnavich+Gvaramadze_%d-%d.eps' % (int(plotRange[0]),int(plotRange[1]))
    print('writing plot to file <'+plotname+'>')
    plt.savefig(plotname, format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
    plt.show()
