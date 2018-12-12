execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

#wcsFits = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/gtc_object_av.fits'
#wcsFitsOut = '/Users/azuri/temp/gtc_object_wcsIm_slit_0000956343.fits'
ims = [['/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/temp/0000956343-20160707-OSIRIS-OsirisLongSlitSpectroscopy1.fits',
        '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/temp/0000956343-20160707-OSIRIS-OsirisLongSlitSpectroscopy2.fits',
        '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-2/gtc_object_av.fits',
        '/Users/azuri/temp/gtc_object_wcsIm_slit_'],
       ['/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/temp/0000956342-20160707-OSIRIS-OsirisLongSlitSpectroscopy1.fits',
        '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/temp/0000956342-20160707-OSIRIS-OsirisLongSlitSpectroscopy2.fits',
        '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-2/gtc_object_av.fits',
        '/Users/azuri/temp/gtc_object_wcsIm_sky_'],
       ['/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/temp/0000956339-20160707-OSIRIS-OsirisLongSlitSpectroscopy1.fits',
        '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/temp/0000956339-20160707-OSIRIS-OsirisLongSlitSpectroscopy2.fits',
        '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/gtc_object_av.fits',
        '/Users/azuri/temp/gtc_object_wcsIm_slit_'],
       ['/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/temp/0000956338-20160707-OSIRIS-OsirisLongSlitSpectroscopy1.fits',
        '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/temp/0000956338-20160707-OSIRIS-OsirisLongSlitSpectroscopy2.fits',
        '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/gtc_object_av.fits',
        '/Users/azuri/temp/gtc_object_wcsIm_sky_'],
]
for im in ims:
    tmpStr = im[0][im[0].rfind('/')+1]
    tmpStr = tmpStr[0:tmpStr.find('-')]
    print('tmpStr = <'+tmpStr+'>')
    wcsFitsOut = im[3]+tmpStr+'.fits'
    iraf.imdel(wcsFitsOut)
    iraf.imcopy(im[2], wcsFitsOut)
    headerFull = getHeader(wcsFitsOut,0)
    headerLeft = getHeader(im[0],0)
    headerRight = getHeader(im[1],0)
    nXFull = headerFull['NAXIS1']
    nYFull = headerFull['NAXIS2']
    nXLeft = headerLeft['NAXIS1']
    nYLeft = headerLeft['NAXIS2']
    nXRight = headerRight['NAXIS1']
    nYRight = headerRight['NAXIS2']
    print 'nXFull = ',nXFull,', nXLeft = ',nXLeft,', nXRight = ',nXRight
    print 'nYFull = ',nYFull,', nYLeft = ',nYLeft,', nYRight = ',nYRight

    copyStrLeft = im[0]+'[1:'+str(nXLeft)+',1:'+str(nYLeft)+'], '+wcsFitsOut+'[1:'+str(nXLeft)+',1:'+str(nYLeft)+']'
    copyStrRight = im[1]+'[1:'+str(nXRight)+',1:'+str(nYRight)+'], '+wcsFitsOut+'['+str(nXFull-nXRight)+':'+str(nXFull)+',1:'+str(nYRight)+']'
    print 'copyStrLeft = ',copyStrLeft
    print 'copyStrRight = ',copyStrRight
    iraf.imcopy(im[0]+'[1:'+str(nXLeft)+',1:'+str(nYLeft)+']', wcsFitsOut+'[0][1:'+str(nXLeft)+',1:'+str(nYLeft)+']')
    iraf.imcopy(im[1]+'[1:'+str(nXRight)+',1:'+str(nYRight)+']', wcsFitsOut+'[0]['+str(nXFull-nXRight)+':'+str(nXFull)+',1:'+str(nYRight)+']')
    #iraf.imedit()
