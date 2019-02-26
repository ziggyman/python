import pyphot
execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

gtc_file = "/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned.fits"
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
print('apparent GTC B mag = ',mag)
f = filter_JB[0].get_flux(gtc_wavelength, gtc_spectrum)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
print('apparent GTC B mag = ',mag)

f = filter_JB[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
print('absolute GTC B mag = ',mag)

f = filter_JV[0].getFlux(gtc_wavelength, gtc_spectrum)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_JV[0].Vega_zero_mag
print('apparent GTC mag V = ',mag)

f = filter_JV[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_JV[0].Vega_zero_mag
print('absolute GTC mag V = ',mag)

f = filter_GaiaG_BP[0].getFlux(gtc_wavelength/10., gtc_spectrum)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].Vega_zero_mag
print('apparent GTC mag G_BP = ',mag)

f = filter_GaiaG_BP[0].getFlux(gtc_wavelength/10., gtc_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].Vega_zero_mag
print('absolute GTC mag G_BP = ',mag)

f = filter_SDSS_u[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_SDSS_u[0].Vega_zero_mag
print('absolute GTC mag SDSS_u = ',mag)

f = filter_SDSS_g[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_SDSS_g[0].Vega_zero_mag
print('absolute GTC mag SDSS_g = ',mag)

f = filter_SDSS_r[0].getFlux(gtc_wavelength, gtc_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_SDSS_r[0].Vega_zero_mag
print('absolute GTC mag SDSS_r = ',mag)
#J_B = pyphot.helpers.extractPhotometry(gtc_wavelength, gtc_spectrum_absolute, filter_JB)#, absFlux=True, progress=True)
#print('J_B = ',J_B)

print(' ')
#WIYN
f = filter_JB[0].getFlux(wiyn_wavelength, wiyn_spectrum)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
print('apparent WIYN mag B = ',mag)

f = filter_JB[0].getFlux(wiyn_wavelength, wiyn_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
print('absolute WIYN mag B = ',mag)

f = filter_JV[0].getFlux(wiyn_wavelength, wiyn_spectrum)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_JV[0].Vega_zero_mag
print('apparent WIYN mag V = ',mag)

f = filter_JV[0].getFlux(wiyn_wavelength, wiyn_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_JV[0].Vega_zero_mag
print('absolute WIYN mag V = ',mag)

f = filter_GaiaG_BP[0].getFlux(wiyn_wavelength/10., wiyn_spectrum)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].Vega_zero_mag
print('apparent WIYN mag G_BP = ',mag)

f = filter_GaiaG_BP[0].getFlux(wiyn_wavelength/10., wiyn_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_GaiaG_BP[0].Vega_zero_mag
print('absolute WIYN mag G_BP = ',mag)

f = filter_SDSS_u[0].getFlux(wiyn_wavelength/10., wiyn_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_SDSS_u[0].Vega_zero_mag
print('absolute WIYN mag SDSS_u = ',mag)

f = filter_SDSS_g[0].getFlux(wiyn_wavelength/10., wiyn_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_SDSS_g[0].Vega_zero_mag
print('absolute WIYN mag SDSS_g = ',mag)

f = filter_SDSS_r[0].getFlux(wiyn_wavelength/10., wiyn_spectrum_absolute)
print('f = ',f)
mag = -2.5 * np.log10(f) - filter_SDSS_r[0].Vega_zero_mag
print('absolute WIYN mag SDSS_r = ',mag)


f=1540.1
mag = -2.5 * np.log10(f)# - filter_JB[0].Vega_zero_mag
print('GAIA mag = ',mag)
mag = -2.5 * np.log10(f) - filter_JB[0].Vega_zero_mag
print('GAIA mag = ',mag)

