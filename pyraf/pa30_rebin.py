execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

calibrated = True
wavelengthRange = [4280.0,7092.0]
dlambda = 1.0

inFiles = ["/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_sum_n.fits",
           "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_sum_n.fits"]
if calibrated:
    inFiles = ["/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum.fits",
               "/Users/azuri/daten/uni/HKU/Pa30/Pa30_WN151014_cal_sum.fits"]

for inFile in inFiles:
    hdulist = pyfits.open(inFile)

    header = hdulist[0].header
    spectrum = fits.getdata(inFile)
    wavelength = getWavelength(header,1)

    wavelengthNew, fluxNew = rebin(wavelength, spectrum, wavelengthRange, dlambda)

    # --- write spectrum
    hdulist[0].data = fluxNew
    #        hdulist[0].header['NAXIS1'] = header['NAXIS2']
    hdulist[0].header['CRPIX1'] = 1.
    hdulist[0].header['CRVAL1'] = wavelengthNew[0]
    hdulist[0].header['CDELT1'] = dlambda
    hdulist[0].header['CD1_1'] = dlambda
    hdulist[0].header['WAT1_001'] = 'wtype=linear label=Wavelength units=angstroms'

    specOutName = inFile[0:inFile.rfind('.')] + '_rebinned' + inFile[inFile.rfind('.'):]
    print 'specOutName = <'+specOutName+'>'

    hdulist.writeto(specOutName, clobber=True)
