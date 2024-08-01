from drUtils_yichuan import plotSpec,fitEmissionLines,getImageData,getWavelengthArr,fitEmissionLine

oneSpectrum = 'SCIENCE_MPA1757-3547_dbs01046b_otzxfifEcdF_clean.fits'
#plotSpec(oneSpectrum)

#fitsFileList = 'fitsFiles.list'
#plotSpec('@'+fitsFileList)

#fitEmissionLines(oneSpectrum,oneSpectrum[:oneSpectrum.rfind('.')]+'_out.fits')

wLen = getWavelengthArr(oneSpectrum,0)
spec = getImageData(oneSpectrum,0)
xRange = [5000,5014]

with open('results.txt','w') as f:
    area = fitEmissionLine(wLen,spec,xRange,Display=True)
    print('area = ',area)
    f.write('%s: [OIII] 5007 = %.5E\n' % (oneSpectrum,area))
