execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

#tab = [
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0011/gtc_object_av_wl_flt.fits', 'ylim' : [-30000.,50000.], 'OName': 'IPHASXJ190333_Green_Slit', 'SkyLeft':[612,702], 'SkyRight':[1823,1855], 'ObjectAreas':[[1090,1125],[1164,1212],[1238,1252],[1267,1295],[1307,1353],[1367,1444],[1495,1516],[1530,1556],[1600,1650],[1665,1700]]},#faint line at 5005, 6565, 6585, 6829, 6862
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0019/gtc_object_av_wl_flt.fits', 'ylim' : [-50000.,180000], 'OName': 'We2-260', 'SkyLeft':[404,580], 'SkyRight':[1782,1982], 'ObjectAreas':[[800,957],[978,1038],[1090,1240],[1262,1716]]},#line at 6582
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0020/gtc_object_av_wl_flt.fits', 'ylim' : [-3000.,25000.], 'OName': 'IPHASXJ225420', 'SkyLeft':[536,696], 'SkyRight':[1571,1659], 'ObjectAreas':[[1269,1295],[1311,1359]]},#line at 6559, 6580
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0028/gtc_object_av_wl_flt.fits', 'ylim' : [-10000,79000.], 'OName': '', 'SkyLeft':[688,956], 'SkyRight':[1764,1878], 'ObjectAreas':[[974,1038],[1128,1188],[1206,1288],[1308,1514],[1534,1554],[1570,1662]]},#lines at 5005, 6560, 6580
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0043/gtc_object_av_wl_flt_cal.fits', 'ylim' : [,], 'OName': '', 'SkyLeft':[612,0], 'SkyRight':[0,0], 'ObjectAreas':[[0,0],[0,0]]},#
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0049/gtc_object_av_wl_flt.fits', 'ylim' : [-23000.,70000.], 'OName': 'IPHASXJ193849', 'SkyLeft':[428,534], 'SkyRight':[1836,1916], 'ObjectAreas':[[944,998],[1012,1038],[1090,1110],[1138,1150],[1168,1220],[1254,1270],[1284,1354],[1382,1412],[1440,1460],[1486,1510],[1542,1552],[1572,1586],[1602,1610],[1626,1652]]},#lines at 6580
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0067/gtc_object_av_wl_flt.fits', 'ylim' : [,], 'OName': 'ldu18', 'SkyLeft':[1122,1170], 'SkyRight':[1388,1504], 'ObjectAreas':[[1266,1344]]},#lines at 4956, 5005
#]

#tab=[{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0043/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-1000,100000], 'OName': '', 'SkyLeft':[612,0], 'SkyRight':[0,0], 'ObjectAreas':[[0,0],[0,0]]}]
tab=[{'FName': 'Laurence/GTC4-16AMEX/OB0065-1/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1415,1444], 'SkyRight':[1509,1529], 'ObjectAreas':[[1443,1488]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-1a/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1415,1433], 'SkyRight':[1509,1529], 'ObjectAreas':[[1434,1490]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-1b/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1283,1308], 'SkyRight':[1416,1436], 'ObjectAreas':[[1323,1342]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-1c/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1283,1308], 'SkyRight':[1416,1436], 'ObjectAreas':[[1386,1402]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-1d/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1260,1277], 'SkyRight':[1563,1586], 'ObjectAreas':[[1478,1562]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-1e/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1563,1586], 'SkyRight':[1744,1774], 'ObjectAreas':[[1650,1710]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-2/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1264,1320], 'SkyRight':[1501,1539], 'ObjectAreas':[[1406,1421]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-2a/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1264,1320], 'SkyRight':[1501,1539], 'ObjectAreas':[[1372,1461]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-2b/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1200,1264], 'SkyRight':[1369,1382], 'ObjectAreas':[[1278,1296]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-2c/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1264,1320], 'SkyRight':[1369,1382], 'ObjectAreas':[[1337,1355]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-2d/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1300,1320], 'SkyRight':[1568,1600], 'ObjectAreas':[[1433,1537]]},
     {'FName': 'Laurence/GTC4-16AMEX/OB0065-2e/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[1568,1600], 'SkyRight':[1666,1686], 'ObjectAreas':[[1604,1663]]},
    ]

iSpec = 0
spectra = []
for obs in tab:
    directory = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/',obs['FName'][0:obs['FName'].rfind('/')])
    print 'directory = <'+directory+'>'
    print dir(iraf)
#    iraf.user.reduce_gtcmos(indirs=directory)
    image_file = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/',obs['FName'])
    try:
        hdulist = pyfits.open(image_file)
    except:
        print 'file <'+image_file+'> not found'
        image_file = image_file[0:image_file.rfind('_x')]+image_file[image_file.rfind('_x')+2:]
        print 'trying <'+image_file+'>'
        try:
            hdulist = pyfits.open(image_file)
        except:
            print 'file <'+image_file+'> not found'

#    print 'type(hdulist) = ',type(hdulist)
    header = hdulist[0].header
    if header['OBJECT'] != '':
        obs['OName'] = header['OBJECT']
    try:
        print 'airmass = ',header['AIRMASS']
    except:
        """do nothing"""
    print obs
    wavelength = getWavelength(header)
    image_data = fits.getdata(image_file)
#        print 'type(image_data) = ',type(image_data)
#        print 'type(image_data[1000,1000]) = ',type(image_data[1000,1000])
#        print 'len(image_data) = ',len(image_data)
#        print 'image_data.shape = ',image_data.shape
#        print 'image_data[1551,',obs['ObjectAreas'][0][0]+1,'] = ',image_data[1551,obs['ObjectAreas'][0][0]+1]
    imageMinusSky, imageSky = subtractSky(image_data,obs['SkyLeft'],obs['SkyRight'])
#        print 'imageMinusSky[1551,',obs['ObjectAreas'][0][0]+1,':',obs['ObjectAreas'][0][0]+5,'] = ',imageMinusSky[1551,obs['ObjectAreas'][0][0]+1:obs['ObjectAreas'][0][0]+5]

    obsCols = populateObjectArray(imageMinusSky,obs['ObjectAreas'])
#        print 'obsCols[1551,1:5] = ',obsCols[1551,1:5]
#        print 'obsCols[1000,0] = %.3e' % (obsCols[1000,0])
#        print 'sum(obsCols[1551,:]) = ',np.sum(obsCols[1551,:])

    obsColsSmoothed = boxCarMedianSmooth(obsCols, 0, 5)
#        print 'sum(obsColsSmoothed[1551,:]) = ',np.sum(obsColsSmoothed[1551,:])
#        print 'obsColsSmoothed[1551,1:5] = ',obsColsSmoothed[1551,1:5]
#        print 'imageMinusSky[1551,',obs['ObjectAreas'][0][0]+1,':',obs['ObjectAreas'][0][0]+5,'] = ',imageMinusSky[1551,obs['ObjectAreas'][0][0]+1:obs['ObjectAreas'][0][0]+5]
#        print 'obsCols[1551,1:5] = ',obsCols[1551,1:5]
#        print 'obsColsSmoothed[1551,1:5] = ',obsColsSmoothed[1551,1:5]

    spectrum = np.ndarray(imageSky.shape[0], dtype=np.float32)
    for iRow in range(spectrum.shape[0]):
        spectrum[iRow] = np.sum(obsColsSmoothed[iRow,:])
#            print 'spectrum[',iRow,'] = ',spectrum[iRow]

    plt.plot(wavelength,spectrum)
    plt.xlabel('wavelength [A]')
    yLabel = 'calibrated flux [erg/cm2/s/A]'
    if obs['FName'].find("_cal") < 0:
        print '"_cal" not found in obs["FName"] = <',obs['FName'],'> => Setting uLabel to counts [ADU]'
        yLabel = 'counts [ADU]'
    plt.ylabel(yLabel)
    if obs['ylim'][1] > 0:
        plt.ylim(obs['ylim'][0],obs['ylim'][1])
    plt.title(obs['OName'] + obs['FName'][obs['FName'].rfind('-')+1:obs['FName'].rfind('/')])
    plt.savefig(image_file[0:image_file.rfind('.')]+'_spectrum.png')
    plt.show()

    # --- mark sky and object areas in original image
    markAreas(image_data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
    markAreas(imageMinusSky, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])

    hdulist[0].data = obsCols
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_obsCols.fits', clobber=True)

    hdulist[0].data = obsColsSmoothed
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_obsColsSmoothed.fits', clobber=True)

    hdulist[0].data = image_data
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_areasMarked.fits', clobber=True)

    hdulist[0].data = imageSky
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_sky.fits', clobber=True)

    hdulist[0].data = imageMinusSky
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky.fits', clobber=True)

    hdulist[0].data = imageMinusSky
    nCols = 0
    for area in obs['ObjectAreas']:
        hdulist[0].data[:,area[0]:area[1]] = obsCols[:,nCols:nCols+area[1]-area[0]]
        nCols += area[1]-area[0]
    markAreas(hdulist[0].data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky_obs_not_smoothed.fits', clobber=True)

    hdulist[0].data = imageMinusSky
    nCols = 0
    for area in obs['ObjectAreas']:
        hdulist[0].data[:,area[0]:area[1]] = obsColsSmoothed[:,nCols:nCols+area[1]-area[0]]
        nCols += area[1]-area[0]
    markAreas(hdulist[0].data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky_obs_smoothed.fits', clobber=True)

    # --- cut off blue end
    idx = findFirstIdxWithValGT(wavelength, 3740.)
    wavelength = wavelength[idx:]
    spectrum = spectrum[idx:]
#        print 'wavelength.shape = ',wavelength.shape
    print 'spectrum.shape = ',spectrum.shape

    # --- interpolate over 5567-5590 (5577 OII sky line)
    idxL = findFirstIdxWithValGT(wavelength, 5567.)
    idxH = findFirstIdxWithValGT(wavelength, 5590.)
    spectrum, linRegPars = subtractSky(spectrum,[idxL-5,idxL],[idxH,idxH+5])

    # --- interpolate over 5867-5905 (~5885 OII sky line)
    idxL = findFirstIdxWithValGT(wavelength, 5867.)
    idxH = findFirstIdxWithValGT(wavelength, 5905.)
    spectrum, linRegPars = subtractSky(spectrum,[idxL-5,idxL],[idxH,idxH+5])

    # --- interpolate over 6290-6308 (~6299 OII sky line)
    idxL = findFirstIdxWithValGT(wavelength, 6290.)
    idxH = findFirstIdxWithValGT(wavelength, 6308.)
    spectrum, linRegPars = subtractSky(spectrum,[idxL-5,idxL],[idxH,idxH+5])

#    fitsName = image_file[0:image_file.rfind('.')]+'s.fits'
#    write1DFits(spectrum, fitsName)

    # --- continuum normalize spectrum
#    fitsNameOut = fitsName[0:fitsName.rfind('.')]+'c.fits'
#    if os.path.exists(fitsNameOut):
#        os.remove(fitsNameOut)
#    iraf.continuum(fitsName,
#                   fitsNameOut,
#                   ask="yes",
#                   lines="*",
#                   bands="1",
#                   type="ratio",
#                   replace=False,
#                   wavescale=True,
#                   logscale=False,
#                   override=False,
#                   listonly=False,
#                   logfiles="logfile",
#                   interactive=True,
#                   sample="*",
#                   naverage=1,
#                   function="spline3",
#                   order=6,
#                   low_reject=5.,
#                   high_reject=1.,
#                   niterate=10,
#                   grow=1.,
#                   markrej=True)
#    spectrum=read1DFits(fitsNameOut)

#    pymodelfit.fitgui.FitGui(xdata=wavelength, ydata=spectrum, weights=None, model='smoothspline', include_models=None, exclude_models=None, fittype=None, **traits)

    # --- plot cleaned spectrum
    plt.plot(wavelength,spectrum)
    plt.xlabel('wavelength [A]')
    yLabel = 'calibrated flux [erg/cm2/s/A]'
    if obs['FName'].find("_cal") < 0:
        print '"_cal" not found in obs["FName"] = <',obs['FName'],'> => Setting uLabel to counts [ADU]'
        yLabel = 'counts [ADU]'
    plt.ylabel(yLabel)
    if obs['ylim'][1] > 0:
        plt.ylim(obs['ylim'][0],obs['ylim'][1])
    plt.title(obs['OName'])
    plt.savefig(image_file[0:image_file.rfind('.')]+'_spectrum.png')
    plt.show()

    # --- write spectrum
    hdulist[0].data = spectrum
#        hdulist[0].header['NAXIS1'] = header['NAXIS2']
    hdulist[0].header['CRPIX1'] = header['CRPIX2']
    hdulist[0].header['CRVAL1'] = wavelength[0]
    hdulist[0].header['CDELT1'] = header['CDELT2']
    hdulist[0].header['CD1_1'] = header['CD2_2']
    hdulist[0].header['WAT1_001'] = header['WAT2_001']
#        hdulist[0].header[''] = header['']
    obsdate = getDateTime(header['DATE-OBS'])
    datestr = obsdate.strftime('%d%m%y')
    print 'datestr = <'+datestr+'>'
    specOutName = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/'+obs['FName'][0:obs['FName'].rfind('/')],obs['OName']+'_GT'+datestr+'.fits')
    print 'specOutName = <'+specOutName+'>'

    hdulist.writeto(specOutName, clobber=True)

    with open(image_file[0:image_file.rfind('.')]+'_spec.dat','w') as f:
        for i in range(len(wavelength)):
            f.write('%.5e %.5e\n' % (wavelength[i], spectrum[i]))
    spectra.append([wavelength, spectrum])
    if iSpec >= len(tab) / 2:
        plt.plot(spectra[iSpec][0],spectra[iSpec][1],'g-')
        plt.plot(spectra[iSpec-int(len(tab)/2)][0],spectra[iSpec-int(len(tab)/2)][1],'b-')
        plt.xlabel('wavelength [A]')
        yLabel = 'calibrated flux [erg/cm2/s/A]'
        if obs['FName'].find("_cal") < 0:
            print '"_cal" not found in obs["FName"] = <',obs['FName'],'> => Setting uLabel to counts [ADU]'
            yLabel = 'counts [ADU]'
        plt.ylabel(yLabel)
        xMin = 6650.
        xMax = 7200.
        plt.xlim(xMin,xMax)
        yMin = 1.
        yMax = 0.
        specNew = spectrum
        for i in range(len(wavelength)):
            if spectra[iSpec - int(len(tab) / 2)][1][i] < specNew[i]:
                specNew[i] = spectra[iSpec - int(len(tab) / 2)][1][i]
            if np.abs(wavelength[i] - spectra[iSpec - int(len(tab) / 2)][0][i]) > 0.01:
                print 'PROBLEM: np.abs(wavelength[i](=',wavelength[i],') - spectra[iSpec - int(len(tab) / 2)][0][i](=',spectra[iSpec - int(len(tab) / 2)][0][i],'))(=',np.abs(wavelength[i] - spectra[iSpec - int(len(tab) / 2)][0][i]),') > 0.01'
            if (wavelength[i] > xMin) and (wavelength[i] < xMax):
                if spectrum[i] > yMax:
                    yMax = spectrum[i]
                if spectra[iSpec - int(len(tab) / 2)][1][i] > yMax:
                    yMax = spectra[iSpec - int(len(tab) / 2)][1][i]
                if spectra[iSpec - int(len(tab) / 2)][1][i] < yMin:
                    yMin = spectra[iSpec - int(len(tab) / 2)][1][i]
                if spectrum[i] < yMin:
                    yMin = spectrum[i]
        plt.plot(wavelength,specNew,'r-')

        if yMin > 0.:
            yMin = 0.
        plt.ylim(yMin, yMax)
        plt.title(obs['OName'])
        plt.show()

        hdulist[0].data = specNew
        specOutName = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/'+obs['FName'][0:obs['FName'].rfind('/')],obs['OName']+'_GT'+datestr+'_min.fits')
        print 'specOutName = <'+specOutName+'>'

        hdulist.writeto(specOutName, clobber=True)


        # --- plot cleaned spectrum
        plt.plot(wavelength,specNew)
        plt.xlabel('wavelength [A]')
        yLabel = 'calibrated flux [erg/cm2/s/A]'
        if obs['FName'].find("_cal") < 0:
            print '"_cal" not found in obs["FName"] = <',obs['FName'],'> => Setting uLabel to counts [ADU]'
            yLabel = 'counts [ADU]'
        plt.ylabel(yLabel)
        if obs['ylim'][1] > 0:
            plt.ylim(obs['ylim'][0],obs['ylim'][1])
        plt.title(obs['OName'])
        plt.show()

        # --- continuum normalize spectrum
        fitsNameOut = specOutName[0:specOutName.rfind('.')]+'c.fits'
        if os.path.exists(fitsNameOut):
            os.remove(fitsNameOut)
        if False:
            iraf.continuum(specOutName,
                           fitsNameOut,
                           ask="yes",
                           lines="*",
                           bands="1",
                           type="ratio",
                           replace=False,
                           wavescale=True,
                           logscale=False,
                           override=False,
                           listonly=False,
                           logfiles="logfile",
                           interactive=True,
                           sample="*",
                           naverage=1,
                           function="spline3",
                           order=6,
                           low_reject=5.,
                           high_reject=1.,
                           niterate=10,
                           grow=1.,
                           markrej=True)

    iSpec += 1
#    STOP
#    plt.imshow(image_data, cmap='gray')
#    plt.colorbar()
#    i+=1
#    if i == 3:
#        STOP
#print tab
