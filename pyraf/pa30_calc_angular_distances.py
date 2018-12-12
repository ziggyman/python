from __future__ import print_function, division
import math
from PyAstronomy import pyasl
import subprocess

execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

fitsFiles = ['/Users/azuri/temp/gtc_object_wcsIm_sky_0000956338.fits',
             '/Users/azuri/temp/gtc_object_wcsIm_slit_0000956339.fits',
             '/Users/azuri/temp/gtc_object_wcsIm_sky_0000956342.fits',
             '/Users/azuri/temp/gtc_object_wcsIm_slit_0000956343.fits']

#centralStarXY = [1460, 994]
centralStarRaDec = ['00:53:11.20', '+67:30:02.40']
radius = ['00:53:27.20', '+67:30:02.40']

#print('hmsToDeg(0:00:01.0) = ',hmsToDeg('0:00:01.0'))
#print('hmsToDeg(1:00:00.0) = ',hmsToDeg('1:00:00.0'))
#print('hmsToDeg(0:59:59.9) = ',hmsToDeg('0:59:59.9'))
#print('hmsToDeg(23:59:59.9) = ',hmsToDeg('23:59:59.9'))
#print('hmsToDeg(22:59:59.9) = ',hmsToDeg('22:59:59.9'))

for fitsFile in fitsFiles:
    result = subprocess.check_output(['sky2xy', fitsFile, centralStarRaDec[0], centralStarRaDec[1]])
    print(type(result))
    print('sky2xy(',centralStarRaDec[0],', ', centralStarRaDec[1],') = ',result)
    raCenter, decCenter, xCenter, yCenter = getRaDecXY(result)

    print('radius = <',radius,'>')
    print('radius[0] = ',radius[0], ', radius[1] = ',radius[1])
    result = subprocess.check_output(['sky2xy', fitsFile, radius[0], radius[1]])
    print('sky2xy(',radius[0],', ', radius[1],') = ',result)
    raRadius, decRadius, xRadius, yRadius = getRaDecXY(result)

    #result = subprocess.check_output(['xy2sky', '-j', fitsFile, str(centralStarXY[0]), str(centralStarXY[1])])
    #print('xy2sky(',centralStarXY[0],', ',centralStarXY[1],') = ',result)

    try:
        hdulist = pyfits.open(fitsFile)
    except:
        print('file <'+fitsFile+'> not found')

    #    print 'type(hdulist) = ',type(hdulist)
    header = hdulist[0].header
    radiusLen = math.sqrt((xCenter - xRadius)**2) + math.sqrt((yCenter - yRadius)**2)
    x0 = xCenter
    x0Int = int(xCenter)
    y0 = yCenter
    y0Int = int(yCenter)
    for x in np.arange(x0Int-10, x0Int+10, 1):
        hdulist[0].data[y0Int,x] = 65000.
    for y in np.arange(y0Int - 10, y0Int+10, 1):
        hdulist[0].data[y,x0Int] = 65000.

    for alpha in np.arange(0., 2. * math.pi, 2. * math.pi / 3000.):
        hdulist[0].data[int(y0 + (radiusLen * math.cos(alpha))), int(x0 + (radiusLen * math.sin(alpha)))] = 65000.
    hdulist.writeto(fitsFile[0:fitsFile.rfind('.')]+'_marked.fits', clobber=True)


    print('angular distance from [xCenter=',xCenter,', yCenter=',yCenter,'] to [xRadius=',xRadius,', yRadius=',yRadius,'] = ',angularDistanceFromXY(fitsFile, xCenter, yCenter, xRadius, yRadius) * 3600.)
#    print('angular distance from [1068, 1273] to [1459, 994] = ',angularDistanceFromXY(fitsFile, 1068., 1273., 1459., 994.) * 3600.)
#    print('angular distance from [977, 994] to [1471, 988] = ',angularDistanceFromXY(fitsFile, 977., 994., 1471., 988.) * 3600.)
#    print('angular distance from [977, 994] to [1289, 674] = ',angularDistanceFromXY(fitsFile, 977., 994., 1289., 674.) * 3600.)
