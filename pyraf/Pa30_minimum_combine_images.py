execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

im1 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/gtc_object_av_x_wl_flt_cal.fits'
im2 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-2/gtc_object_av_x_wl_flt_cal.fits'
imOut = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_av_x_wl_flt_cal_minComb.fits'

offset = 49
gap = [1038,1088]
ignoreFirst = 159
ignoreLast = 81

imMinCombine(im1, im2, offset, gap, imOut, True, ignoreFirst, ignoreLast)
