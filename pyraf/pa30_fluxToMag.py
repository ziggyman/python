import pyphot
from pyphot import unit
import astropy.io.fits as pyfits
from astropy.time import Time
from datetime import datetime as dt
import math
import numpy as np
from matplotlib import pyplot as plt
from myUtils import getWavelength,smooth,toYearFraction
from Pa30_LBT import readLBTFiles

photometry_file = '/Users/azuri/daten/uni/HKU/Pa30/variability/pa30_photometry_visual_1.txt'

gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_GT080716_cal_sum_cleaned.fits"
wiyn_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_WN151014_cal_sum_cleaned_t.fits"

gvaramadze_file = "/Users/azuri/daten/uni/HKU/Pa30/variability/gvaramadze.fits"
gvaramadze_hdulist = pyfits.open(gvaramadze_file)
gvaramadze_header = gvaramadze_hdulist[0].header
gvaramadze_wavelength = np.array(getWavelength(gvaramadze_header,1)) * unit['AA']
gvaramadze_spectrum = np.array(pyfits.getdata(gvaramadze_file)) * unit['erg/s/cm**2/AA']
#gvaramadze_wavelength,gvaramadze_spectrum = readGvaramadzeFile()

garnavich_wavelength, garnavich_spectrum = readLBTFiles()
garnavich_wavelength = np.array(garnavich_wavelength)
garnavich_spectrum = np.array(garnavich_spectrum)
garnavich_spectrum = garnavich_spectrum[garnavich_wavelength < 5455.]
garnavich_wavelength = garnavich_wavelength[garnavich_wavelength < 5455.] * unit['AA']
#garnavich_spectrum_smoothed = smooth(garnavich_spectrum,9)#scipy.ndimage.mean_filter(garnavich_spectrum, 7)#boxCarMeanSmooth(somme_spectrum, 0, 21)
garnavich_spectrum = garnavich_spectrum * unit['erg/s/cm**2/AA']

gtc_hdulist = pyfits.open(gtc_file)
wiyn_hdulist = pyfits.open(wiyn_file)

gtc_header = gtc_hdulist[0].header
gtc_wavelength = getWavelength(gtc_header,1) * unit['AA']
gtc_spectrum = pyfits.getdata(gtc_file) * unit['erg/s/cm**2/AA']
print('gtc_wavelength = ',gtc_wavelength)

wiyn_header = wiyn_hdulist[0].header
wiyn_wavelength = getWavelength(wiyn_header,1) * unit['AA']
wiyn_spectrum = pyfits.getdata(wiyn_file) * unit['erg/s/cm**2/AA']
print('wiyn_wavelength = ',wiyn_wavelength)

with open(photometry_file,'r') as f:
    lines = f.readlines()
photometry = []
for line in lines:
    if line[0] != '#':
        phot = line.strip().split()
        phot[0] = Time(float(phot[0]), format='mjd').decimalyear#toYearFraction(dt.fromtimestamp(float(phot[0])))
        phot[1] = float(phot[1])
        phot[2] = float(phot[2])
        photometry.append(phot)
print('photometry = ',photometry)

lib = pyphot.get_library()
print("Library contains: ", len(lib), " filters")
print("dir(lib) = ",dir(lib))
for i in range(len(lib.content)):
    print('lib.content[',i,'] = ',lib.content[i])

filters = lib.load_filters(['GROUND_JOHNSON_B',
                            'GROUND_JOHNSON_V',
                            'SDSS_u',
                            'SDSS_g',
                            'SDSS_r',
                            'PS1_g',
                            'PS1_r',
                            'PS1_i',
                            'GaiaDR2_BP',
                            'GaiaDR2_RP',
                            'GaiaDR2_G',
#                            '2MASS_H',
#                            '2MASS_J',
#                            '2MASS_Ks',
                            ])

print('dir(filters[0]) = ',dir(filters[0]))
print('filters[0].name = ',filters[0].name)

synphot = []
for specName,date,wlen,spec in [['wiyn',dt(2014,10,15),wiyn_wavelength,wiyn_spectrum],
                                ['gtc',dt(2016,7,8),gtc_wavelength,gtc_spectrum],
                                ['gvaramadze',dt(2017,7,20),gvaramadze_wavelength,gvaramadze_spectrum],
                                ['garnavich',dt(2020,9,11),garnavich_wavelength,garnavich_spectrum]]:
    for filter in filters:
        f = filter.getFlux(wlen, spec)
#    print('f = ',f)
#    print('dir(f) = ',dir(f))
#    print('f.value = ',f.value)
#    print('f.dimensionless = ',f.dimensionless)
#    print('f.unitless = ',f.unitless)
#    print('f.magnitude = ',f.magnitude)
#       
        magV = -2.5 * np.log10(f.value) - filter.Vega_zero_mag
        magAB = -2.5 * np.log10(f.value) - filter.AB_zero_mag
        print(specName,': ',filter.name,': apparent mag Vega = ',magV)
        print(specName,': ',filter.name,': apparent mag AB = ',magAB)
        synphot.append({'magVega': magV, 'magAB': magAB, 'filter_name': filter.name, 'spec_name': specName, 'date':toYearFraction(date)})
#        mag = -2.5 * np.log10(f.value) - filter.AB_zero_mag
#        print('apparent GTC B mag AB = ',mag)

print('synphot = ',len(synphot),': ',synphot)

labels=['WIYN','GTC','Gvaramadze','Garnavich']
for i in range(len(filters)):
    synphots = [x for x in [synphot[i],
                            synphot[len(filters)+i],
                            synphot[(2*len(filters))+i],
                            synphot[(3*len(filters))+i]]]
    dates = [x['date'] for x in synphots]
    Vmags = [x['magVega'] for x in synphots]
    ABmags = [x['magAB'] for x in synphots]
    filter_names = [x['filter_name'] for x in synphots]
    spec_names = [x['spec_name'] for x in synphots]
    print('dates = ',dates)
    print('Vmags = ',Vmags)
    print('ABmags = ',ABmags)
    print('filter_names = ',filter_names)
    print('spec_names = ',spec_names)
    for j in range(4):
        if not math.isinf(ABmags[j]):
            plt.scatter(dates[j],
                        ABmags[j],
                        label=labels[j])
            print('plotted [',dates[j],', ',ABmags[j],'], label = ',labels[j])

    for phot in photometry:
        if (filter_names[0] == 'GROUND_JOHNSON_B') and (phot[3] == 'JohnsonB'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'GROUND_JOHNSON_V') and (phot[3] == 'JohnsonV'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'SDSS_g') and (phot[3] == 'IGAPS.Gunn_g(UVEX)'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'SDSS_r') and (phot[3] == 'IGAPS.Gunn_r(UVEX)'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'PS1_g') and (phot[3] == 'PS1.g'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'PS1_r') and (phot[3] == 'PS1.r'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'PS1_i') and (phot[3] == 'PS1.i'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'GaiaDR2_BP') and (phot[3] == 'phot_bp_mean_mag'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'GaiaDR2_RP') and (phot[3] == 'phot_rp_mean_mag'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
        elif (filter_names[0] == 'GaiaDR2_G') and (phot[3] == 'phot_g_mean_mag'):
            plt.scatter(phot[0],phot[1],label=phot[3])
            #plt.errorbar(phot[0],phot[1], yerr=phot[2], fmt="o")
            print('plotted [',phot[0],', ',phot[1],'], label = ',phot[3])
    plt.title(filter_names[0])
    plt.xlabel('year')
    plt.ylabel('apparent magnitude')
    plt.legend()
    plt.gca().invert_yaxis()
    plotname = '/Users/azuri/daten/uni/HKU/Pa30/variability/Pa30_synphot'+filter_names[0]+'.eps'
    plt.savefig(plotname, format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
    plt.show()

if False:
    filters = lib.load_filters(['GROUND_JOHNSON_B',
                            'GROUND_JOHNSON_V',
                            'SDSS_u',
                            'SDSS_g',
                            'SDSS_r',
                            'PS1_g',
                            'PS1_r',
                            'PS1_i',
                            'GaiaDR2_BP',
                            'GaiaDR2_RP',
                            'GaiaDR2_G',
#                            '2MASS_H',
#                            '2MASS_J',
#                            '2MASS_Ks',
                            ])
    execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

    gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned_scaled.fits"
    wiyn_file = "/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean_52_cal_cleaned.fits"

    gtc_hdulist = pyfits.open(gtc_file)
    gtc_header = gtc_hdulist[0].header
    gtc_wavelength = getWavelength(gtc_header,1)
    gtc_spectrum = fits.getdata(gtc_file)
    gtc_spectrum_absolute = calibratedFluxToAbsoluteFlux(gtc_spectrum, 3370.0)
    print('gtc_wavelength = ',gtc_wavelength.size,': ',gtc_wavelength)
    print('gtc_spectrum = ',gtc_spectrum.size,': ',gtc_spectrum)
    print('gtc_spectrum_absolute = ',gtc_spectrum_absolute.size,': ',gtc_spectrum_absolute)

    wiyn_hdulist = pyfits.open(wiyn_file)
    wiyn_header = wiyn_hdulist[0].header
    wiyn_wavelength = getWavelength(wiyn_header,1)
    wiyn_spectrum = fits.getdata(wiyn_file)
    wiyn_spectrum_absolute = calibratedFluxToAbsoluteFlux(wiyn_spectrum, 3370.0)

    lib = pyphot.get_library()
    print("Library contains: ", len(lib), " filters")
    print("dir(lib) = ",dir(lib))
    for i in range(len(lib.content)):
        print('lib.content[',i,'] = ',lib.content[i])

    filter_JB = lib.load_filters(['GROUND_JOHNSON_B'])
    print('filter_JB = ',filter_JB)
    print('filter_JB[0].info() = ',filter_JB[0].info())

    filter_JV = lib.load_filters(['GROUND_JOHNSON_V'])
    print('filter_JV = ',filter_JV)
    print('filter_JV[0].info() = ',filter_JV[0].info())

    filter_GaiaG_BP = lib.load_filters(['GaiaDR2_BP'])
    print('dir(filter_GaiaG_BP[0]) = ',dir(filter_GaiaG_BP[0]))
    print('filter_GaiaG_BP[0].wavelength = ',filter_GaiaG_BP[0].wavelength)
    print('filter_GaiaG_BP[0].wavelength_unit = ',filter_GaiaG_BP[0].wavelength_unit)
    print('filter_GaiaG_BP[0].info() = ',filter_GaiaG_BP[0].info())

    filter_SDSS_u = lib.load_filters(['SDSS_u'])
    print('filter_SDSS_u = ',filter_SDSS_u[0])
    print('filter_SDSS_u[0].info() = ',filter_SDSS_u[0].info())

    filter_SDSS_g = lib.load_filters(['SDSS_g'])
    print('filter_SDSS_g = ',filter_SDSS_g[0])
    print('filter_SDSS_g[0].info() = ',filter_SDSS_g[0].info())

    filter_SDSS_r = lib.load_filters(['SDSS_r'])
    print('filter_SDSS_r = ',filter_SDSS_r[0])
    print('filter_SDSS_r[0].info() = ',filter_SDSS_r[0].info())

    f = filter_JB[0].getFlux(gtc_wavelength, gtc_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
    print('apparent GTC B mag Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JB[0].AB_zero_mag
    print('apparent GTC B mag AB = ',mag)
    f = filter_JB[0].get_flux(gtc_wavelength, gtc_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
    print('apparent GTC B mag Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JB[0].AB_zero_mag
    print('apparent GTC B mag AB = ',mag)

    f = filter_JB[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
    print('absolute GTC B mag Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JB[0].AB_zero_mag
    print('absolute GTC B mag AB = ',mag)

    f = filter_JV[0].getFlux(gtc_wavelength, gtc_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JV[0].Vega_zero_mag
    print('apparent GTC mag V Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JV[0].AB_zero_mag
    print('apparent GTC mag V AB = ',mag)

    f = filter_JV[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JV[0].Vega_zero_mag
    print('absolute GTC mag V Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JV[0].AB_zero_mag
    print('absolute GTC mag V AB = ',mag)

    f = filter_GaiaG_BP[0].getFlux(gtc_wavelength/10., gtc_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].Vega_zero_mag
    print('apparent GTC mag G_BP Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].AB_zero_mag
    print('apparent GTC mag G_BP AB = ',mag)

    f = filter_GaiaG_BP[0].getFlux(gtc_wavelength/10., gtc_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].Vega_zero_mag
    print('absolute GTC mag G_BP Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].AB_zero_mag
    print('absolute GTC mag G_BP AB = ',mag)

    f = filter_SDSS_u[0].getFlux(gtc_wavelength, gtc_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_u[0].Vega_zero_mag
    print('apparent GTC mag SDSS_u Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_u[0].AB_zero_mag
    print('apparent GTC mag SDSS_u AB = ',mag)

    f = filter_SDSS_u[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_u[0].Vega_zero_mag
    print('absolute GTC mag SDSS_u Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_u[0].AB_zero_mag
    print('absolute GTC mag SDSS_u AB = ',mag)

    f = filter_SDSS_g[0].getFlux(gtc_wavelength, gtc_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_g[0].Vega_zero_mag
    print('apparent GTC mag SDSS_g Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_g[0].AB_zero_mag
    print('apparent GTC mag SDSS_g AB = ',mag)

    f = filter_SDSS_g[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_g[0].Vega_zero_mag
    print('absolute GTC mag SDSS_g Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_g[0].AB_zero_mag
    print('absolute GTC mag SDSS_g AB = ',mag)

    f = filter_SDSS_r[0].getFlux(gtc_wavelength, gtc_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_r[0].Vega_zero_mag
    print('apparent GTC mag SDSS_r Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_r[0].AB_zero_mag
    print('apparent GTC mag SDSS_r AB = ',mag)

    f = filter_SDSS_r[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_r[0].Vega_zero_mag
    print('absolute GTC mag SDSS_r Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_r[0].AB_zero_mag
    print('absolute GTC mag SDSS_r AB = ',mag)
    #J_B = pyphot.helpers.extractPhotometry(gtc_wavelength, gtc_spectrum_absolute, filter_JB)#, absFlux=True, progress=True)
    #print('J_B = ',J_B)

    print(' ')
    #WIYN
    f = filter_JB[0].getFlux(wiyn_wavelength, wiyn_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
    print('apparent WIYN mag B Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JB[0].AB_zero_mag
    print('apparent WIYN mag B AB = ',mag)

    f = filter_JB[0].getFlux(wiyn_wavelength, wiyn_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
    print('absolute WIYN mag B Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JB[0].AB_zero_mag
    print('absolute WIYN mag B AB = ',mag)

    f = filter_JV[0].getFlux(wiyn_wavelength, wiyn_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JV[0].Vega_zero_mag
    print('apparent WIYN mag V Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JV[0].AB_zero_mag
    print('apparent WIYN mag V AB = ',mag)

    f = filter_JV[0].getFlux(wiyn_wavelength, wiyn_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_JV[0].Vega_zero_mag
    print('absolute WIYN mag V Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_JV[0].AB_zero_mag
    print('absolute WIYN mag V AB = ',mag)

    f = filter_GaiaG_BP[0].getFlux(wiyn_wavelength/10., wiyn_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].Vega_zero_mag
    print('apparent WIYN mag G_BP Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].AB_zero_mag
    print('apparent WIYN mag G_BP AB = ',mag)

    f = filter_GaiaG_BP[0].getFlux(wiyn_wavelength/10., wiyn_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].Vega_zero_mag
    print('absolute WIYN mag G_BP Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].AB_zero_mag
    print('absolute WIYN mag G_BP AB = ',mag)

    f = filter_SDSS_g[0].getFlux(wiyn_wavelength, wiyn_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_g[0].Vega_zero_mag
    print('apparent WIYN mag SDSS_g Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_g[0].AB_zero_mag
    print('apparent WIYN mag SDSS_g AB = ',mag)

    f = filter_SDSS_g[0].getFlux(wiyn_wavelength, wiyn_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_g[0].Vega_zero_mag
    print('absolute WIYN mag SDSS_g Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_g[0].AB_zero_mag
    print('absolute WIYN mag SDSS_g AB = ',mag)

    f = filter_SDSS_r[0].getFlux(wiyn_wavelength, wiyn_spectrum)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_r[0].Vega_zero_mag
    print('apparent WIYN mag SDSS_r Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_r[0].AB_zero_mag
    print('apparent WIYN mag SDSS_r AB = ',mag)

    f = filter_SDSS_r[0].getFlux(wiyn_wavelength, wiyn_spectrum_absolute)
    print('f = ',f)
    mag = -2.5 * np.log10(f) - filter_SDSS_r[0].Vega_zero_mag
    print('absolute WIYN mag SDSS_r Vega = ',mag)
    mag = -2.5 * np.log10(f) - filter_SDSS_r[0].AB_zero_mag
    print('absolute WIYN mag SDSS_r AB = ',mag)

