import pyfits
import astropy.io.fits as fits
import os
import numpy as np
import matplotlib.pyplot as plt
from myUtils import getWavelength, subtractSky, populateObjectArray, boxCarMedianSmooth, markAreas, findFirstIdxWithValGT, getDateTime
#%matplotlib inline

#execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate
tab = [
{'FName': 'Laurence/GTC4-16AMEX/OB0065-1/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[980,1034], 'SkyRight':[1965,2028], 'ObjectAreas':[[1506,1557]]},#very wide, find sky from other imeage nearby
#'FName': 'Laurence/GTC4-16AMEX/OB0065-2/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[980,1034], 'SkyRight':[1965,2028], 'ObjectAreas':[[1387,1397]]},#very wide, find sky from other imeage nearby
#'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0001/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [0,0], 'OName': 'IPHASXJ194645', 'SkyLeft':[1114,1156], 'SkyRight':[1351,1418], 'ObjectAreas':[[1262,1330],[1157,1240]]},#very wide, find sky from other imeage nearby
#{'FName': 'IPHAS_GTC_DATA/DATA2_2017A/GTC12-17AMEX/OB0001/gtc_object_av_x_wl_flt_cal.fits', 'ylim' : [-0.15e-15,0.3e-15], 'OName': 'IPHASXJ194645', 'SkyLeft':[835,960], 'SkyRight':[1570,1760], 'ObjectAreas':[[1262,1330],[1157,1240]]},#very wide, find sky from other imeage nearby
#{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0036/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-0.3e-15,0.8e-15], 'OName': 'IPHASXJ214032', 'SkyLeft':[1291,1303], 'SkyRight':[1359,1380], 'ObjectAreas':[[1312,1356]]},#lines at 4961 and 5005, no central star visible
##{'FName': 'IPHAS_GTC_DATA/DATA4_2016A/GTC4-16AMEX/OB0036/gtc_object_av_wl_flt_cal.fits', 'ylim' : [-0.3e-15,0.8e-15], 'OName': 'IPHASXJ214032', 'SkyLeft':[1090,1250], 'SkyRight':[1552,1661], 'ObjectAreas':[[1312,1356]]},#lines at 4961 and 5005, no central star visible
]
i=1
for obs in tab:
    image_file = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/',obs['FName'])
    hdulist = pyfits.open(image_file)
#    print 'type(hdulist) = ',type(hdulist)
    header = hdulist[0].header
    if header['OBJECT'] != '':
        obs['OName'] = header['OBJECT']
    try:
        print('airmass = ',header['AIRMASS'])
    except:
        """do nothing"""
    print(obs)
    if True:
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

        obsColsSmoothed = boxCarMedianSmooth(obsCols, 0, 15)
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
            print('"_cal" not found in obs["FName"] = <',obs['FName'],'> => Setting uLabel to counts [ADU]')
            yLabel = 'counts [ADU]'
        plt.ylabel(yLabel)
        if obs['ylim'][1] > 0:
            plt.ylim(obs['ylim'][0],obs['ylim'][1])
        plt.title(obs['OName'])
        plt.savefig(image_file[0:image_file.rfind('.')]+'_spectrum.png')
        plt.show()

        # --- mark sky and object areas in original image
        markAreas(image_data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
        markAreas(imageMinusSky, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])

        hdulist[0].data = obsCols
        hdulist.writeto(image_file[0:image_file.rfind('.')]+'_obsCols1.fits', clobber=True)

        hdulist[0].data = obsColsSmoothed
        hdulist.writeto(image_file[0:image_file.rfind('.')]+'_obsColsSmoothed1.fits', clobber=True)

        hdulist[0].data = image_data
        hdulist.writeto(image_file[0:image_file.rfind('.')]+'_areasMarked1.fits', clobber=True)

        hdulist[0].data = imageSky
        hdulist.writeto(image_file[0:image_file.rfind('.')]+'_sky1.fits', clobber=True)

        hdulist[0].data = imageMinusSky
        hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky1.fits', clobber=True)

        hdulist[0].data = imageMinusSky
        nCols = 0
        for area in obs['ObjectAreas']:
            hdulist[0].data[:,area[0]:area[1]] = obsCols[:,nCols:nCols+area[1]-area[0]]
            nCols += area[1]-area[0]
        markAreas(hdulist[0].data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
        hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky_obs_not_smoothed1.fits', clobber=True)

        hdulist[0].data = imageMinusSky
        nCols = 0
        for area in obs['ObjectAreas']:
            hdulist[0].data[:,area[0]:area[1]] = obsColsSmoothed[:,nCols:nCols+area[1]-area[0]]
            nCols += area[1]-area[0]
        markAreas(hdulist[0].data, obs['SkyLeft'], obs['SkyRight'], obs['ObjectAreas'])
        hdulist.writeto(image_file[0:image_file.rfind('.')]+'_mSky_obs_smoothed1.fits', clobber=True)

        # --- cut off blue end
        idx = findFirstIdxWithValGT(wavelength, 3700.)
        wavelength = wavelength[idx:]
        spectrum = spectrum[idx:]
#        print 'wavelength.shape = ',wavelength.shape
#        print 'spectrum.shape = ',spectrum.shape

        # --- interpolate over 5567-5590 (5577 OII sky line)
        idxL = findFirstIdxWithValGT(wavelength, 5567.)
        idxH = findFirstIdxWithValGT(wavelength, 5590.)
        spectrum, linRegPars = subtractSky(spectrum,[idxL-5,idxL],[idxH,idxH+5])

        # --- interpolate over 5567-5590 (5577 OII sky line)
        idxL = findFirstIdxWithValGT(wavelength, 5867.)
        idxH = findFirstIdxWithValGT(wavelength, 5905.)
        spectrum, linRegPars = subtractSky(spectrum,[idxL-5,idxL],[idxH,idxH+5])

        # --- plot cleaned spectrum
        plt.plot(wavelength,spectrum)
        plt.xlabel('wavelength [A]')
        yLabel = 'calibrated flux [erg/cm2/s/A]'
        if obs['FName'].find("_cal") < 0:
            print('"_cal" not found in obs["FName"] = <',obs['FName'],'> => Setting uLabel to counts [ADU]')
            yLabel = 'counts [ADU]'
        plt.ylabel(yLabel)
        if obs['ylim'][1] > 0:
            plt.ylim(obs['ylim'][0],obs['ylim'][1])
        plt.title(obs['OName'])
        plt.savefig(image_file[0:image_file.rfind('.')]+'_spectrum1.png')
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
        print('datestr = <'+datestr+'>')
        specOutName = os.path.join('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/',obs['OName']+'_GT'+datestr+'1.fits')

        hdulist.writeto(specOutName, clobber=True)

        with open(image_file[0:image_file.rfind('.')]+'_spec1.dat','w') as f:
            for i in range(len(wavelength)):
                f.write('%.5e %.5e\n' % (wavelength[i], spectrum[i]))
#        STOP
#    plt.imshow(image_data, cmap='gray')
#    plt.colorbar()
#    i+=1
#    if i == 3:
#        STOP
#print tab