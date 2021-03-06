execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

calibrated = False

im1 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1f/gtc_object_av_x_wl_flt_cal_mSky_obs_not_smoothed.fits'
im2 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-2f/gtc_object_av_x_wl_flt_cal_mSky_obs_not_smoothed.fits'
#im1 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/gtc_object_av_x_wl_flt_cal.fits'
#im2 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-2/gtc_object_av_x_wl_flt_cal.fits'
imOut = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_av_x_wl_flt_cal_mSky_obs_not_smoothed_minComb.fits'

if not calibrated:
    im1 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1f/gtc_object_av_x_wl_flt_mSky_obs_not_smoothed.fits'
    im2 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-2f/gtc_object_av_x_wl_flt_mSky_obs_not_smoothed.fits'
    #im1 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-1/gtc_object_av_x_wl_flt_cal.fits'
    #im2 = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065-2/gtc_object_av_x_wl_flt_cal.fits'
    imOut = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_av_x_wl_flt_mSky_obs_not_smoothed_minComb.fits'

offset = 49
gap = [1038,1088]
ignoreFirst = 159
ignoreLast = 81

imMinCombine(im1, im2, offset, gap, imOut, True, ignoreFirst, ignoreLast)
