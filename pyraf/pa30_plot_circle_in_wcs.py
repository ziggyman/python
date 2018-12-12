import math

execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

for image_file in ['/Users/azuri/temp/gtc_object_wcsIm_slit_0000956343.fits', '/Users/azuri/temp/gtc_object_wcsIm_sky_0000956342.fits']:
    try:
        hdulist = pyfits.open(image_file)
    except:
        print 'file <'+image_file+'> not found'

    #    print 'type(hdulist) = ',type(hdulist)
    header = hdulist[0].header
    radius = math.sqrt((1469.676 - 1289.156)**2) + math.sqrt((986.946 - 674.243)**2)
    x0 = 1469.676
    x0Int = 1470
    y0 = 986.946
    y0Int = 987
    for x in np.arange(x0Int-10, x0Int+10, 1):
        hdulist[0].data[y0Int,x] = 65000.
    for y in np.arange(y0Int - 10, y0Int+10, 1):
        hdulist[0].data[y,x0Int] = 65000.

    for alpha in np.arange(0., 2. * math.pi, 2. * math.pi / 3000.):
        hdulist[0].data[int(y0 + (radius * math.cos(alpha))), int(x0 + (radius * math.sin(alpha)))] = 65000.
    hdulist.writeto(image_file[0:image_file.rfind('.')]+'_marked.fits', clobber=True)

