execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

#wavelengthRange = [4280.0,7092.0]
#nPix = 4000
#dlambda = (wavelengthRange[1] - wavelengthRange[0]) / (nPix-1)
#print('dlambda = ',dlambda)

inFilesList = "/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum.list"#objects_botzfxsEcBld.list"

inFiles = open(inFilesList,'r')
path = inFilesList[0:inFilesList.rfind('/')]
print('path = <'+path+'>')

iFile = 0
for inFile in inFiles:
    print('inFile = <'+inFile+'>')
    inFile = inFile.strip('\n')
    inFile = os.path.join(path,inFile)
    print('inFile = <'+inFile+'>')
    hdulist = pyfits.open(inFile)

    header = hdulist[0].header
#    print('header = ',header)
#    if iFile == 0:
#        wavelengthRange = getMaximumWavelengthRange(header, dispAxis=1)
#        nPix = 4000
#        dlambda = (wavelengthRange[1] - wavelengthRange[0]) / (nPix-1)
#        print('new wavelengthRange = ',wavelengthRange)
#        print('dlambda = ',dlambda)
#        iFile = 1
    spectrum = fits.getdata(inFile)
    print('spectrum.shape = ',spectrum.shape)
    wavelength = getWavelength(header,1)
    print('wavelength.shape = ',wavelength.shape)
    print('wavelength = ',wavelength)

    medianSpec = np.median(spectrum,axis=0)
    print('medianSpec.shape = ',medianSpec.shape)
#    STOP

#    wavelengthNew, fluxNew = rebin(wavelength, spectrum[0,:], wavelengthRange, dlambda)
#    specNew = np.zeros((spectrum.shape[0], fluxNew.shape[0]), dtype=type(spectrum[0,0]))

    for iSpec in range(spectrum.shape[0]):
#        wavelength = getWavelengthMultiSpec(header,iSpec,1)
#        print('iSpec = ',iSpec,': wavelength = ',wavelength.shape,': ',wavelength)
#        wavelengthNew, fluxNew = rebin(wavelength, spectrum[iSpec,:], wavelengthRange, dlambda)
#        print('wavelengthNew = ',wavelengthNew)
#        print('fluxNew = ',fluxNew.shape,': ',fluxNew)

#        specNew[iSpec,:] = fluxNew

        if iSpec == 51:
            plt.plot(wavelength, spectrum[iSpec,:], 'r-', label='object')
#            plt.plot(wavelengthNew, specNew[iSpec,:], 'b-', label='rebinned')
        spectrum[iSpec,:] = spectrum[iSpec,:] - medianSpec
        if iSpec == 51:
            plt.plot(wavelength, spectrum[iSpec,:], 'b-', label='object-sky')
            plt.legend()
            plt.show()

    # --- write spectrum
#    print('wavelengthNew = ',wavelengthNew)
#    print('specNew = ',specNew)
    hdulist[0].data = spectrum
    #        hdulist[0].header['NAXIS1'] = header['NAXIS2']
    #hdulist[0].header['CRPIX1'] = 1.
#    hdulist[0].header['CRVAL1'] = wavelength[0]
#    hdulist[0].header['CDELT1'] = dlambda
#    hdulist[0].header['CD1_1'] = dlambda
#    hdulist[0].header['WAT1_001'] = 'wtype=linear label=Wavelength units=angstroms'

    specOutName = inFile[0:inFile.rfind('.')] + '-skyMedian' + inFile[inFile.rfind('.'):]
    print 'specOutName = <'+specOutName+'>'

    hdulist.writeto(specOutName, clobber=True)
