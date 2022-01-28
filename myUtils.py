from __future__ import print_function, division

from astropy.coordinates import SkyCoord
import astropy.io.fits as apyfits
import astropy.units as u
from datetime import date, datetime
#import hammer
import itertools
import math
import matplotlib.pyplot as plt
import numpy as np
import os
#import pymodelfit
#from pyraf import iraf
import re
import subprocess
from time import sleep
from hammer import Pixel,XY,LonLat,Hammer

c0 = 299792.458 #km/s

def multikeysort(items, columns):
    from operator import itemgetter
    comparers = [((itemgetter(col[1:].strip()), -1) if col.startswith('-') else
                  (itemgetter(col.strip()), 1)) for col in columns]
    def comparer(left, right):
        for fn, mult in comparers:
            result = cmp(fn(left), fn(right))
            if result:
                return mult * result
        else:
            return 0

def countObjectColumns(areas):
    nCols = 0
    for area in areas:
        nCols += area[1]-area[0]
    return nCols

def findFirstIdxWithValGT(valueArr, val):
    idx = 0
    while valueArr[idx] < val:
        idx += 1
    return idx

def populateSkyArray(imageData, area):
    sky = np.copy(imageData[:,area[0]:area[1]])
    return sky

def populateObjectArray(imageData, areas):
    nCols = countObjectColumns(areas)
    print('nCols = ',nCols)
    obsCols = np.zeros(shape=(imageData.shape[0],nCols), dtype=type(imageData[0,0]))
    nCols = 0
    for area in areas:
#        print 'imageData[1000,',area[0],':',area[1],'] = ',imageData[1000,area[0]:area[1]]
        obsCols[:,nCols:nCols+area[1]-area[0]] = np.copy(imageData[:,area[0]:area[1]])
        nCols += area[1]-area[0]
    return obsCols

# --- fit y = a*x + b
# --- return value w: w[0] = a, w[1] = b
def linReg(x,y):
    oneD = False
    if (len(y.shape) < 2):
        oneD = True
    A = np.array([ x, np.ones(len(x))])
#    print 'linReg: A = ',A
#    print 'linReg: A.T = ',A.T
#    print 'linReg: x.shape = ',x.shape,', A.shape = ',A.shape,', y.shape = ',y.shape
    if not oneD:
        w = []
        for iy in range(y.shape[0]):
            w.append(np.linalg.lstsq(A.T,y[iy,:])[0])
    #        print 'linReg: iy = ',iy,': w[',iy,'] = ',w[iy]
    else:
#        print 'A.shape = ',A.shape
#        print 'y.shape = ',y.shape
        w = np.linalg.lstsq(A.T,y)[0]
    return w

# axis: 0 (columns) or 1 (rows)
# width: odd number
def boxCarMedianSmooth(imageData, axis, width):
    print('boxCarMedianSmooth: imageData.shape = ',imageData.shape)
    print('boxCarMedianSmooth: len(imageData.shape) = ',len(imageData.shape))
    print('boxCarMedianSmooth: imageData = ',imageData)
    newDataArray = None
    if len(imageData.shape) > 1:
        newDataArray = np.zeros(shape=imageData.shape, dtype=type(imageData[0,0]))
        if axis == 0:
            for iRow in range(imageData.shape[0]):
                for iCol in range(imageData.shape[1]):
                    if iCol < int(width/2.0):
                        iColStart = 0
                        iColEnd = iCol+int(width/2.0)+1
                    elif iCol > imageData.shape[1]-int(width/2.0)-1:
                        iColStart = iCol - int(width/2.0)
                        iColEnd = imageData.shape[1]
                    else:
                        iColStart = iCol - int(width/2.0)
                        iColEnd = iCol + int(width/2.0) + 1
                    newDataArray[iRow,iCol] = np.median(imageData[iRow,iColStart:iColEnd])
    #                print 'iRow = ',iRow,', iCol = ',iCol,': iColStart = ',iColStart,', iColEnd = ',iColEnd,': imageData[iRow,iColStart:iColEnd] = ',imageData[iRow,iColStart:iColEnd],': median = ',newDataArray[iRow,iCol]
        elif axis == 1:
            for iCol in range(imageData.shape[1]):
                for iRow in range(imageData.shape[0]):
                    if iRow < int(width/2.0):
                        iRowStart = 0
                        iRowEnd = iRow+int(width/2.0)+1
                    elif iRow > imageData.shape[0]-int(width/2.0)-1:
                        iRowStart = iRow - int(width/2.0)
                        iRowEnd = imageData.shape[1]
                    else:
                        iRowStart = iRow - int(width/2.0)
                        iRowEnd = iRow + int(width/2.0) + 1
                    newDataArray[iRow,iCol] = np.median(imageData[iRowStart:iRowEnd,iCol])
    #                print 'iCol = ',iCol,', iRow = ',iRow,': iRowStart = ',iRowStart,', iRowEnd = ',iRowEnd,': imageData[iRowStart:iRowEnd,iCol] = ',imageData[iRowStart:iRowEnd,iCol],': median = ',newDataArray[iRow,iCol]
        else:
            print('boxCarMedianSmooth: ERROR: axis(=',axis,') out of bounds [0,1]')
    else:
        newDataArray = np.zeros(shape=imageData.shape, dtype=type(imageData[0]))
        for iCol in range(imageData.shape[0]):
            if iCol < int(width/2.0):
                iColStart = 0
                iColEnd = iCol+int(width/2.0)+1
            elif iCol > imageData.shape[0]-int(width/2.0)-1:
                iColStart = iCol - int(width/2.0)
                iColEnd = imageData.shape[0]
            else:
                iColStart = iCol - int(width/2.0)
                iColEnd = iCol + int(width/2.0) + 1
            newDataArray[iCol] = np.median(imageData[iColStart:iColEnd])
    return newDataArray

def boxCarMeanSmooth(imageData, axis, width):
    print('use scipy.ndimage.median_filter instead')
    STOP
    print('imageData.shape = ',imageData.shape)
    print('len(imageData.shape) = ',len(imageData.shape))
    print('imageData = ',imageData)
    newDataArray = None
    if len(imageData.shape) > 1:
        newDataArray = np.zeros(shape=imageData.shape, dtype=type(imageData[0,0]))
        if axis == 0:
            for iRow in range(imageData.shape[0]):
                for iCol in range(imageData.shape[1]):
                    if iCol < int(width/2.0):
                        iColStart = 0
                        iColEnd = iCol+int(width/2.0)+1
                    elif iCol > imageData.shape[1]-int(width/2.0)-1:
                        iColStart = iCol - int(width/2.0)
                        iColEnd = imageData.shape[1]
                    else:
                        iColStart = iCol - int(width/2.0)
                        iColEnd = iCol + int(width/2.0) + 1
                    newDataArray[iRow,iCol] = np.mean(imageData[iRow,iColStart:iColEnd])
    #                print 'iRow = ',iRow,', iCol = ',iCol,': iColStart = ',iColStart,', iColEnd = ',iColEnd,': imageData[iRow,iColStart:iColEnd] = ',imageData[iRow,iColStart:iColEnd],': mean = ',newDataArray[iRow,iCol]
        elif axis == 1:
            for iCol in range(imageData.shape[1]):
                for iRow in range(imageData.shape[0]):
                    if iRow < int(width/2.0):
                        iRowStart = 0
                        iRowEnd = iRow+int(width/2.0)+1
                    elif iRow > imageData.shape[0]-int(width/2.0)-1:
                        iRowStart = iRow - int(width/2.0)
                        iRowEnd = imageData.shape[1]
                    else:
                        iRowStart = iRow - int(width/2.0)
                        iRowEnd = iRow + int(width/2.0) + 1
                    newDataArray[iRow,iCol] = np.mean(imageData[iRowStart:iRowEnd,iCol])
    #                print 'iCol = ',iCol,', iRow = ',iRow,': iRowStart = ',iRowStart,', iRowEnd = ',iRowEnd,': imageData[iRowStart:iRowEnd,iCol] = ',imageData[iRowStart:iRowEnd,iCol],': mean = ',newDataArray[iRow,iCol]
        else:
            print('ERROR: axis(=',axis,') out of bounds [0,1]')
    else:
        newDataArray = np.zeros(shape=imageData.shape, dtype=type(imageData[0]))
        for iCol in range(imageData.shape[0]):
            if iCol < int(width/2.0):
                iColStart = 0
                iColEnd = iCol+int(width/2.0)+1
            elif iCol > imageData.shape[0]-int(width/2.0)-1:
                iColStart = iCol - int(width/2.0)
                iColEnd = imageData.shape[0]
            else:
                iColStart = iCol - int(width/2.0)
                iColEnd = iCol + int(width/2.0) + 1
            newDataArray[iCol] = np.mean(imageData[iColStart:iColEnd])
    return newDataArray

def sigmaRej(values, sigLow, sigHigh, adjustSigLevels, useMean=False):
    ySkyMedian = None
    if useMean:
        ySkyMedian = np.mean(values)
    else:
        ySkyMedian = np.median(values)
    sigma = np.std(values)
    nRej = 0
    if adjustSigLevels:
        if sigma > 3. * ySkyMedian:
            sigLow = 0.5
            sigHigh = 0.05
    for i in range(len(values)):
        if (values[i] < ySkyMedian - (sigLow * sigma)) or (values[i] > ySkyMedian + (sigHigh * sigma)):
            values[i] = ySkyMedian
            nRej += 1
#            print 'subtractSky: median = ',ySkyMedian,', sigma = ',sigma,': removed pixel ',i,' from sky'
#    print 'sigmaRej: rejected ',nRej,' out of ',len(values),' pixels'
    return values

# --- in case of 1D image data replace the area [skyLeftArea[1]:skyRightArea[0]] with the interpolated 'sky'
# --- otherwise subtract the 2D sky image from the 2D image data
def subtractSky(imageData,skyLeftArea,skyRightArea,sigLow=3.0,sigHigh=3.0):
    dtype = None
    axis = None
    oneD = False
    if (len(imageData.shape) > 1):
        dtype = type(imageData[0,0])
        axis = 1
    else:
        dtype = type(imageData[0])
        oneD = True
        axis = 0
#    print 'oneD = ',oneD
    newImageData = np.ndarray(imageData.shape, dtype=dtype)
    skyData = np.ndarray(imageData.shape, dtype=dtype)
    xSky = np.concatenate([np.arange(skyLeftArea[0],skyLeftArea[1],1.0),np.arange(skyRightArea[0],skyRightArea[1],1.0)])
#    print 'subtractSky: xSky = ',xSky
    skyParams = None
    if (oneD):
        skyLeft = imageData[skyLeftArea[0]:skyLeftArea[1]]
        skyLeftSmoothed = boxCarMedianSmooth(skyLeft, 0, 5)
        skyRight = imageData[skyRightArea[0]:skyRightArea[1]]
        skyRightSmoothed = boxCarMedianSmooth(skyRight, 0, 5)
        ySky = np.concatenate([skyLeftSmoothed, skyRightSmoothed])
        ySky = sigmaRej(ySky, sigLow, sigHigh, True)
        skyParams = linReg(xSky,ySky)
    else:
        skyLeft = populateSkyArray(imageData, skyLeftArea)
        skyLeftSmoothed = boxCarMedianSmooth(skyLeft, 0, 9)
        skyRight = populateSkyArray(imageData, skyRightArea)
        skyRightSmoothed = boxCarMedianSmooth(skyRight, 0, 9)
        ySky = np.concatenate([skyLeftSmoothed, skyRightSmoothed],axis=1)
#    print 'subtractSky ySky = ',len(ySky),': ',ySky
#    print 'xSky.shape = ',xSky.shape
#    print 'ySky.shape = ',ySky.shape
#    skyParams = linReg(xSky,ySky)
#    print 'type(skyParams) = ',type(skyParams)
#    print 'subtractSky: imageData.shape = ',imageData.shape,', newImageData.shape = ',newImageData.shape,', xSky.shape = ',xSky.shape,', len(skyParams) = ',len(skyParams),', len(skyParams[0]) = ',len(skyParams[0])
    xObs = np.arange(0.0,imageData.shape[axis],1.0,dtype=dtype)
#    print 'subtractSky: xObs = ',xObs.shape,': ',xObs

    # in case of oneD only replace the "object area" with the sky, otherwise subtract the sky from each row
    if oneD:
        skyData[:] = (skyParams[0] * xObs) + skyParams[1]
        newImageData[:] = imageData[:]
        newImageData[skyLeftArea[1]:skyRightArea[0]] = skyData[skyLeftArea[1]:skyRightArea[0]]
    else:
        for iRow in range(imageData.shape[0]):
#            print 'subtractSky: iRow = ',iRow
            ySkyRow = sigmaRej(ySky[iRow,:], sigLow, sigHigh, True)
            skyParams = linReg(xSky,ySkyRow)
    #        print 'subtractSky: imageData[',iRow,',:] = ',imageData[iRow,:],', newImageData[',iRow,',:].shape = ',newImageData[iRow,:].shape,', skyParams[',iRow,'] = ',len(skyParams[iRow]),': ',skyParams[iRow]
            skyData[iRow,:] = (skyParams[0] * xObs) + skyParams[1]
    #        print 'subtractSky: skyData[',iRow,',:] = ',skyData[iRow,:]
            newImageData[iRow,:] = imageData[iRow,:] - skyData[iRow,:]
    #        print 'subtractSky: newImageData[',iRow,',:] = ',newImageData[iRow,:]
    return [newImageData,skyData]

#def subtractSky(obsAreas,obsData,skyLeftArea,skyLeft,skyRightArea,skyRight):
#    newObsData = np.ndarray(obsData.shape, dtype=type(obsData[0,0]))
##    xSkyLeft = np.arange(skyLeftArea[0],skyLeftArea[1],1.0)
##    xSkyRight = np.arange(skyRightArea[0],skyRightArea[1],1.0)
#    xSky = np.concatenate([np.arange(skyLeftArea[0],skyLeftArea[1],1.0),np.arange(skyRightArea[0],skyRightArea[1],1.0)])
#    print 'subtractSky: xSky = ',xSky
#    nCols = 0
#    for area in obsAreas:
#        nCols += area[1]-area[0]
#    print 'subtractSky: nCols = ',nCols
#
#    xObs = np.ndarray(nCols, dtype=type(obsData[0,0]))
#    nCols = 0
#    for area in obsAreas:
#        xObs[nCols:nCols+area[1]-area[0]] = np.arange(area[0],area[1],1.0)
#        print 'subtractSky: area = ',area,': nCols = ',nCols,': xObs[',nCols,':',nCols+area[1]-area[0],'] = ',xObs[nCols:nCols+area[1]-area[0]]
#        nCols += area[1]-area[0]
#    print 'subtractSky: nCols = ',nCols
#    ySky = np.concatenate([skyLeft, skyRight],axis=1)
#    print 'subtractSky = ',len(ySky),': ',ySky
#    skyParams = linReg(xSky,ySky)
#    print 'type(skyParams) = ',type(skyParams)
#    print 'subtractSky: obsData.shape = ',obsData.shape,', newObsData.shape = ',newObsData.shape,', xObs.shape = ',xObs.shape,', len(skyParams) = ',len(skyParams),', len(skyParams[0]) = ',len(skyParams[0])
#    for iRow in range(obsData.shape[0]):
#        print 'subtractSky: obsData.shape = ',obsData.shape,', newObsData.shape = ',newObsData.shape,', obsData[',iRow,',:].shape = ',obsData[iRow,:].shape,', xObs.shape = ',xObs.shape,', len(skyParams) = ',len(skyParams)
#        newObsData[iRow,:] = obsData[iRow,:] - ((skyParams[iRow][1] * xObs) + skyParams[iRow][0])
#        print 'subtractSky: newObsData[',iRow,',:] = ',newObsData[iRow,:]
#    STOP
#    return newObsData

def getWavelength(header, axis=2):
    nPix = int(header['NAXIS'+str(axis)])
    crPix = int(header['CRPIX'+str(axis)])
    crVal = float(header['CRVAL'+str(axis)])
    cDelt = float(header['CDELT'+str(axis)])
#    lam = np.ndarray(nPix, dtype=np.float32)
#    lam[0] =
#    for i in np.arange(1,nPix):
    lam = ((np.arange(nPix)+1.0)-crPix)*cDelt+crVal
#    lam = np.arange(crVal, crVal + ((nPix-0.5)*cDelt), cDelt, dtype=np.float32)
#    print 'getWavelength: lam = ',len(lam),': ',lam
    return lam

def markAreas(imageData, skyAreaLeft, skyAreaRight, objectAreas):
    maxImageData = np.amax(imageData)
    minImageData = np.amin(imageData)
    imageData[:,skyAreaLeft[0]] = minImageData
    imageData[:,skyAreaLeft[1]] = minImageData
    imageData[:,skyAreaRight[0]] = minImageData
    imageData[:,skyAreaRight[1]] = minImageData
    for area in objectAreas:
        imageData[:,area[0]] = maxImageData
        imageData[:,area[1]] = maxImageData

def writeFits(imData, fname, clobber=True):
    hdu = apyfits.PrimaryHDU(imData)
    hdu.writeto(fname, clobber=clobber)

def getImageData(fname,hduNum=1):
    hdulist = apyfits.open(fname)
    scidata = hdulist[hduNum].data
    hdulist.close()
    return scidata

def specCombineMinimum(imA, imB):
    imOut = np.ndarray(shape = imA.shape, dtype = type(imA[0]))
    for i in range(len(imA)):
        imOut[i] = imA[i]
        if imB[i] < imA[i]:
            imOut[i] = imB[i]
    return imOut

def getHeader(imName, hdu=1):
    hdulist = apyfits.open(imName, ignore_missing_end=True)
    #print('len(hdulist) = ',len(hdulist))
    hdulist[hdu].verify('fix')
    header = hdulist[hdu].header
    return header

def setHeaderKeyWord(fName, keyWord, value, hdu=0):
    apyfits.setval(fName, keyWord, value=value, ext=hdu, ignore_missing_end=True, output_verify='ignore')

# --- get date from string of the form yyyy-mm-ddThh:mm:ss.mmm
def getDate(dateStr):
    dateStr=dateStr[:dateStr.find('T')]
    day=date(int(dateStr[:dateStr.find('-')]),
             int(dateStr[dateStr.find('-')+1:dateStr.rfind('-')]),
             int(dateStr[dateStr.rfind('-')+1:]))
    print('day = ',day)
    return day

def getDateTime(dateTimeStr):
    dateStr = dateTimeStr[:dateTimeStr.find('T')]
    timeStr = dateTimeStr[dateTimeStr.find('T')+1:]
    time=datetime(int(dateStr[:dateStr.find('-')]),
             int(dateStr[dateStr.find('-')+1:dateStr.rfind('-')]),
             int(dateStr[dateStr.rfind('-')+1:]),
             int(timeStr[0:timeStr.find(':')]),
             int(timeStr[timeStr.find(':')+1:timeStr.rfind(':')]),
             int(timeStr[timeStr.rfind(':')+1:timeStr.find('.')]))
    print('time = ',time)
    return time

def findClosestInTime(time, times):
    diffMin = 100000000.
    timeOut = times[0]
    timeIdx = 0
    idx = 0
    for timeA in times:
        diff = time - timeA
        timediff = abs(diff.total_seconds())
        print('diff = ',diff,', timediff = ',type(timediff),': ',timediff)
        if timediff < diffMin:
            diffMin = timediff
            timeOut = timeA
            timeIdx = idx
            idx += 1
    return [timeOut, timeIdx, diffMin]

def getRaDecXY(string):
    strs = re.sub( '\s+', ' ', string ).strip()
#    print('getRaDecXY: strs = ',strs)
    strs = strs.rstrip().split(' ')
#    print('getRaDecXY: strs = ',strs)
    return [strs[0], strs[1], float(strs[len(strs)-2]), float(strs[len(strs)-1])]

#@brief convert RA and DEC to Galactic Longitude and Latitude
#@param ra: RA in degrees
#@param dec: DEC in degrees
#@return [lon, lat]
def raDecToLonLat(ra, dec):
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
#    print('c.galactic = ',type(c.galactic),': ',c.galactic)
#    print('dir(c) = ',type(c),': ',dir(c))
#    print('dir(c.galactic) = ',type(c.galactic),': ',dir(c.galactic))
#    print('c.galactic.l = ',type(c.galactic.l),': ',dir(c.galactic.l))
#    print('c.galactic.l.value = ',type(c.galactic.l.value),': ',dir(c.galactic.l.value))
    return [float(c.galactic.l.value),float(c.galactic.b.value)]

def lonLatToRaDec(l, b):
    c = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
    return [float(c.icrs.ra.value),float(c.icrs.dec.value)]

#def lonLatToRaDec(lon, lat):
#    c = SkyCoord(frame="galactic", l="1h12m43.2s", b="+1d12m43s")
#    return [float(c.ra.value),float(c.dec.value)]

#def getPixel(hammerX, hammerY):
#    ham = hammer.Hammer()
#    for pix in ham.getPixels():
#        if ham.isInside(pix, hammerX, hammerY):
#            return pix
#    return None

# string = xx:yy:zz.zzz
def hmsToDeg(string):
    h, m, s = [float(i) for i in string.split(':')]
    return (15. * s / 3600.) + (15. * m / 60.) + (h * 15.)

def degToHMS(degrees):
    h = int(degrees / 15.)
    m = int((degrees - (h * 15.)) * 4.)
    s = (degrees - (m/4.) - (h*15.)) * 240.
    sStr = '%.3f' % (s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (h,m,sStr)

# string = xx:yy:zz.zzz
def dmsToDeg(string):
    d, m, s = [float(i) for i in string.split(':')]
#    print('dmsToDeg: string = <'+string+'>: d = ',d,', m = ',m,', s = ',s)
    if string[0] == '-':
        d = 0. - d
        return 0. - (s / 3600. + m / 60. + d)
    return s / 3600. + m / 60. + d

def degToDMS(degrees):
    d = int(degrees)
    m = int((degrees - d) * 60)
    s = (degrees - d - (m/60.)) * 3600.
    sStr = '%.3f' % abs(s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (d,abs(m),sStr)

def hmsToCoordStr(string):
    h, m, s = [i for i in string.split(':')]
    return h+'h'+m+'m'+s+'s'

def dmsToCoordStr(string):
    d, m, s = [i for i in string.split(':')]
    return d+'d'+m+'m'+s+'s'

def getPNG(lon,lat):
    return str(lon)[:str(lon).find('.')+2].zfill(5)+('+'if lat >= 0. else '-')+str(abs(lat))[:str(abs(lat)).find('.')+2].zfill(4)

# all angles must be in degrees
#@return: angular distance in degrees
def angularDistancePyAsl(ra1, dec1, ra2, dec2):
    from PyAstronomy import pyasl
#    print('ra1 = ',ra1,', dec1 = ',dec1,', ra2 = ',ra2,', dec2 = ',dec2)
    return pyasl.getAngDist(ra1, dec1, ra2, dec2)

def angularDistanceFromXYPyAsl(fitsName, x1, y1, x2, y2):
    result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x1), str(y1)])
    result2 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x2), str(y2)])
#    print('xy2sky(',str(x1),', ',str(y1),') = ',result1)
#    print('xy2sky(',str(x2),', ',str(y2),') = ',result2)
    raHMS1, decHMS1, x1, y1 = getRaDecXY(result1.decode('utf-8'))
    raHMS2, decHMS2, x2, y2 = getRaDecXY(result2.decode('utf-8'))
#    print('raHMS1 = ',raHMS1,', decHMS1 = ',decHMS1)
#    print('raHMS2 = ',raHMS2,', decHMS2 = ',decHMS2)
    ra1 = hmsToDeg(raHMS1)
    dec1 = dmsToDeg(decHMS1)
    ra2 = hmsToDeg(raHMS2)
    dec2 = dmsToDeg(decHMS2)
#    print('ra1 = ',ra1,', dec1 = ',dec1)
#    print('ra2 = ',ra2,', dec2 = ',dec2)
    return angularDistance(ra1, dec1, ra2, dec2)

#@param ra1: RA in degrees
#@param dec1: DEC in degrees
#@param ra2: RA in degrees
#@param dec2: DEC in degrees
#@return: angular distance in arc seconds
def angularDistance(ra1, dec1, ra2, dec2):
    if ':' in str(ra1):
        print('values must be in degrees!')
        STOP
#    print('ra1 = ',ra1,', dec1 = ',dec1,', ra2 = ',ra2,', dec2 = ',dec2)
    mm1 = SkyCoord(ra=ra1, dec=dec1, frame='icrs', unit="deg")
    mm2 = SkyCoord(ra=ra2, dec=dec2, frame='icrs', unit="deg")
    return mm1.separation(mm2).arcsecond

def angularDistanceFromXY(fitsName, x1, y1, x2, y2):
    result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x1), str(y1)])
    result2 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x2), str(y2)])
#    print('xy2sky(',str(x1),', ',str(y1),') = ',result1)
#    print('xy2sky(',str(x2),', ',str(y2),') = ',result2)
    raHMS1, decHMS1, x1, y1 = getRaDecXY(result1.decode('utf-8'))
    raHMS2, decHMS2, x2, y2 = getRaDecXY(result2.decode('utf-8'))
#    print('raHMS1 = ',raHMS1,', decHMS1 = ',decHMS1)
#    print('raHMS2 = ',raHMS2,', decHMS2 = ',decHMS2)
    ra1 = hmsToCoordStr(raHMS1)
    dec1 = dmsToCoordStr(decHMS1)
    ra2 = hmsToCoordStr(raHMS2)
    dec2 = dmsToCoordStr(decHMS2)
#    print('ra1 = ',ra1,', dec1 = ',dec1)
#    print('ra2 = ',ra2,', dec2 = ',dec2)
    return angularDistance(ra1, dec1, ra2, dec2)


def getArcsecDistance(fitsName, x1, y1, x2, y2):
    result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x1), str(y1)])
    result2 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x2), str(y2)])
    #print('xy2sky(',str(x1),', ',str(y1),') = ',result1)
    #print('xy2sky(',str(x2),', ',str(y2),') = ',result2)
    raHMS1, decHMS1, x1, y1 = getRaDecXY(result1)
    raHMS2, decHMS2, x2, y2 = getRaDecXY(result2)
    #print('raHMS1 = ',raHMS1,', decHMS1 = ',decHMS1)
    #print('raHMS2 = ',raHMS2,', decHMS2 = ',decHMS2)
    mm1 = SkyCoord(ra=raHMS1, dec=decHMS1, unit=(u.hourangle, u.deg))
    mm2 = SkyCoord(ra=raHMS2, dec=decHMS2, unit=(u.hourangle, u.deg))
    return mm1.separation(mm2).arcsecond

def degToArcsec(deg):
    return deg * 3600.

# offset < 0: fitsName2 shifted left relative to fitsName1
# offset > 0: fitsName2 shifted right relative to fitsName1 (spectrum appears further left)
# gap: tuple: [x0,x1]
#def imMinCombine(fitsName1,
#                 fitsName2,
#                 offset,
#                 gap=None,
#                 outFileName=None,
#                 scale=None,
#                 ignoreFirst=None,
#                 ignoreLast=None,
#                ):
if False:
    def setRowNaNs(row):
#        print('setRowNaNs: row = ',row)
        row[:offset] = 'nan'
        row[nColsOut-offset:] = 'nan'
        if gap is not None:
            row[gap[0]:gap[1]] = 'nan'
            row[gap[0]+offset:gap[1]+offset] = 'nan'
        if ignoreFirst is not None:
            row[0:ignoreFirst+offset] = 'nan'
        if ignoreLast is not None:
            row[nColsOut-ignoreLast-offset:] = 'nan'
#        print('setRowNaNs: row with nans = ',row)
        return row

        hdulist1 = None
        hdulist2 = None
        data1 = None
        data2 = None
        if offset > 0:
            hdulist1 = apyfits.open(fitsName1)
            hdulist2 = apyfits.open(fitsName2)
        else:
            hdulist2 = apyfits.open(fitsName1)
            hdulist1 = apyfits.open(fitsName2)
            offset = 0 - offset
            print('changed offset to ',offset)
        data1 = hdulist1[0].data
        data2 = hdulist2[0].data
        data1 = np.nan_to_num(data1)
        data2 = np.nan_to_num(data2)
        print('data1.shape = ',data1.shape,', data2.shape = ',data2.shape)
    #    if scale:
    #        data1median = np.median(data1)
    #        data2median = np.median(data2)
    #        print('data1median = ',data1median)
    #        print('data2median = ',data2median)
    #        print('data1median / data2median = ',data1median / data2median,', data2median / data1median = ',data2median / data1median)
    #        if data1median < data2median:
    #            data1 = data1 * data2median / data1median
    #        else:
    #            data2 = data2 * data1median / data2median

        nColsIn = data1.shape[1]
        nRows = data1.shape[0]
        print('nColsIn = ',nColsIn,', nRows = ',nRows)

        outArr = np.ndarray(shape=(nRows, nColsIn+offset), dtype=type(data1[0][0]))
        print('outArr.shape = ',outArr.shape)
        nColsOut = outArr.shape[1]
        print('nColsOut = ',nColsOut)

        if scale:
            medians1 = np.zeros(shape=(data1[:,0].shape), dtype=type(data1[int(nRows/2),int(nColsIn/2)]))
            medians2 = np.zeros(shape=(data1[:,0].shape), dtype=type(data1[int(nRows/2),int(nColsIn/2)]))
            dataToMedian1 = np.zeros(shape=outArr.shape, dtype=type(data1[int(nRows/2),int(nColsIn/2)]))
            dataToMedian2 = np.zeros(shape=outArr.shape, dtype=type(data1[int(nRows/2),int(nColsIn/2)]))
            print('dataToMedian1.shape = ',dataToMedian1.shape)
            for row in range(nRows):
                dataToMedian1[row,:nColsOut-offset] = data1[row, :].copy()
                dataToMedian1[row,:] = setRowNaNs(dataToMedian1[row,:])
                if row == 0:
                    print('row =',row,': after setRowNaNs: nNaNs = ',countNaNs(dataToMedian1[row,:]))
                medians1[row] = np.nanmedian(dataToMedian1[row,:])

                dataToMedian2[row,offset:] = data2[row,:].copy()
                dataToMedian2[row,:] = setRowNaNs(dataToMedian2[row,:])
                medians2[row] = np.nanmedian(dataToMedian2[row,:])

            nNaNs = countNaNs(dataToMedian1[0,:])
            x = np.ndarray(shape=(dataToMedian2[0,:].shape[0] - nNaNs), dtype=type(data1[500,500]))
            iX = 0
            for i in range(dataToMedian1[0,:].shape[0]):
    #            print('dataToMedian1[0,',i,'] = ',dataToMedian1[0,i])
                if not math.isnan(dataToMedian1[0,i]):
                    x[iX] = i
                    iX += 1
    #                print('iX = ',iX)
            y = np.arange(nRows)
            print('medians1 = ',medians1)
            print('medians2 = ',medians2)
            ratio = None
            medians1mean = np.mean(medians1)
            medians2mean = np.mean(medians2)
            print('medians1mean = ',medians1mean,', medians2mean = ',medians2mean)
            if medians1mean < medians2mean:
                ratio = dataToMedian2 / dataToMedian1
            else:
                ratio = dataToMedian1 / dataToMedian2
            z = np.ndarray(shape=(len(y), x.shape[0]), dtype=type(x[0]))
            print('z.shape = ',z.shape)
            iCol = 0
            for col in x:
    #            print('type(col) = ',type(col),', type(iCol) = ',type(iCol))
    #            print('z.shape = ',z.shape,', ratio.shape = ',ratio.shape)
                z[:,iCol] = ratio[:,int(col)]
                iCol += 1
            print('x.shape = ',x.shape,', y.shape = ',np.asarray(y).shape,', z.shape = ',z.shape,', ratio.shape = ',ratio.shape)
            #deg = [1,1]
            pathOut = fitsName1[:fitsName1.rfind('/')]
            if outFileName is not None:
                pathOut = outFileName[:outFileName.rfind('/')]
            polyFitIn = os.path.join(pathOut,'polyFitIn.fits')
            polyFitOut = os.path.join(pathOut,'polyFitOut.fits')
            z = np.nan_to_num(z)
            hdulist1[0].data = ratio
            hdulist1.writeto(polyFitIn, clobber=True)
    #        coeffs = nppolyfit2d(y, x, z, deg)
    #        print('coeffs = ',coeffs)
    #        xFit = np.arange(nColsOut)
    #        zFit = nppolyval2d(xFit, y, coeffs, deg)
    #        hdulist1[0].data = zFit
    #        hdulist1.writeto('/Users/azuri/daten/uni/HKU/Pa30/ratioFit.fits', clobber=True)
            if os.path.exists(polyFitOut):
                os.remove(polyFitOut)
            iraf.imfit.imsurfit(polyFitIn,
                                polyFitOut,
                                3,
                                5,
                                type_ou='fit',
                                function='leg',
                                cross_t='yes',
                                xmedian=10,
                                ymedian=10,
                                median_=50.,
                                lower=2.0,
                                upper=2.0,
                                ngrow=1,
                                niter=3,
                                regions='all',
                                rows='[132:2051]',
                                columns='[209:1038,1138:2029]',
                                border=50,
                                section='',
                                circle='',
                                div_min='INDEF')
            hdulistPolyFitOut = apyfits.open(polyFitOut)
            zFit = hdulistPolyFitOut[0].data
            print('data1.shape = ',data1.shape)
            print('data1.shape = ',data1.shape)
            print('zFit.shape = ',zFit.shape)
            for row in range(nRows):
                for col in range(nColsIn):
    #                print('row = ',row,', col = ',col)
                    if medians1mean < medians2mean:
                        data1[row, col] = data1[row, col] * zFit[row, col]
                    else:
                        data2[row, col] = data2[row, col] * zFit[row, col+offset]

        for col in range(nColsOut):
            if col < offset:
                outArr[:,col] = data1[:,col]
                print('col(=',col,') < offset(=',offset,')')
    #        elif col > (nCols - offset):
    #            print('outArr[:,',col,'].shape = ',outArr[:,col].shape,', data2[:,',col,'] = ',data2[:,col].shape)
    #            outArr[:,col] = data2[:,col]
    #            print('col(=',col,') > (nCols(=',nCols,') - offset(=',offset,')) = ',nCols-offset)
            elif (gap is not None) and (col >= gap[0]) and (col <= gap[1]):
                outArr[:,col] = data2[:,col-offset]
                print('col(=',col,') >= gap[0](=',gap[0],') and col <= gap[1](=',gap[1],'), setting outArr[:,',col,'] to data2[:,',col-offset,']')
            elif (gap is not None) and (col >= (gap[0]+offset)) and (col <= (gap[1]+offset)):
                outArr[:,col] = data1[:,col]
                print('col(=',col,') >= gap[0]+offset(=',gap[0]+offset,') and col <= gap[1]+offset(=',gap[1]+offset,'), setting outArr[:,',col,'] to data1[:,',col,']')
            elif col >= nColsIn:
                outArr[:,col] = data2[:,col-offset]
            else:
                print('col = ',col,': taking minimum of each pixel')
                outArr[:,col] = data1[:,col]
                for row in range(nRows):
                    if data2[row,col-offset] < outArr[row,col]:
                        outArr[row,col] = data2[row,col-offset]
        hdulist1[0].data = outArr
        if outFileName:
            hdulist1.writeto(outFileName, clobber=True)
        return hdulist1

def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def nppolyfit2d(x, y, f, deg):
    from numpy.polynomial import polynomial
    import numpy as np
    X = np.array(np.meshgrid(x,y))
    xx = np.asarray(X[0])
    yy = np.asarray(X[1])
    ff = np.asarray(f)
    degdeg = np.asarray(deg)
    vander = polynomial.polyvander2d(xx, yy, degdeg)
    vander = vander.reshape((-1,vander.shape[-1]))
    ff = ff.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, ff)[0]
    return c.reshape(degdeg+1)

def nppolyval2d(x, y, c, deg):
    from numpy.polynomial import polynomial
    import numpy as np
    X = np.array(np.meshgrid(x,y))
    f = polynomial.polyval2d(X[0], X[1], c)
    #c1 = nppolyfit2d(X[0], X[1], f, deg)
    return f

def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

def countNaNs(x):
    nNaNs = 0
    for i in range(len(x)):
        if math.isnan(x[i]):
            nNaNs += 1
    return nNaNs

def rebin(wavelength, flux, wavelengthRange, dlambda):
    print('rebin: wavelength = ',wavelength)
    print('rebin: flux = ',flux)
    print('rebin: wavelengthRange = ',wavelengthRange,', dlambda = ',dlambda)

    """create new wavelength array and empty flux array"""
    wavelengthNew = np.arange(wavelengthRange[0], wavelengthRange[1], dlambda)
    print('wavelengthNew = ',wavelengthNew)
    print('type(wavelengthNew) = ',type(wavelengthNew))
    print('len(flux) = ',len(flux))
    print('len(wavelength) = ',len(wavelength))
    fluxNew = np.zeros(len(wavelengthNew))
    print('len(fluxNew) = ',len(fluxNew))
    print('len(wavelengthNew) = ',len(wavelengthNew))

    for wavelengthNewIndex in range(len(wavelengthNew)):
        if wavelengthNew[wavelengthNewIndex] < wavelengthRange[0]:
            fluxNew[wavelengthNewIndex] = 0.0
        elif wavelengthNew[wavelengthNewIndex] > wavelengthRange[1]:
            fluxNew[wavelengthNewIndex] = 0.0
        else:
            for i in range(len(wavelength)-1):
    #            print('rebin: wavelength[',i,'] = ',wavelength[i])
                if (wavelength[i] <= wavelengthNew[wavelengthNewIndex]) and (wavelength[i+1] > wavelengthNew[wavelengthNewIndex]):
                    fluxNew[wavelengthNewIndex] = flux[i] + ((wavelengthNew[wavelengthNewIndex] - wavelength[i]) * (flux[i+1] - flux[i]) / (wavelength[i+1] - wavelength[i]))
#                    print('rebin: wavelengthNew[',wavelengthNewIndex,'] = ',wavelengthNew[wavelengthNewIndex],': fluxNew[',wavelengthNewIndex,'] set to ',fluxNew[wavelengthNewIndex])
#    print('rebin: wavelengthNew = ',wavelengthNew)
    return [wavelengthNew, fluxNew]

# --- iSpec=[0...nSpec-1]
def getWavelengthMultiSpec(header, iSpec, axis=1):
    dispStr = header['WAT2_001']
    iWat = 2
    while(True):
        try:
            watStr = 'WAT2_%03i' % (iWat)
#            print('getWavelengthMultiSpec: watStr = <'+watStr+'>')
            dispStr += header[watStr]
            iWat += 1
        except:
            break
#    print('getWavelengthMultiSpec: dispStr = <'+dispStr+'>')
#    a='82 1 0 4275.4858499095 2.1134804589752 1333 0. 49.83 56.22'
#    b='2 0 0 4275.4858499095 2.1134804589752 1333 0. 603.94 610.34'
    startStr = 'spec2'
    endStr = '"'
    specDispStr = dispStr[dispStr.find(startStr):]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDispStr = specDispStr[specDispStr.find(endStr)+1:]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDispStr = specDispStr[:specDispStr.find(endStr)]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDisp = specDispStr.split(' ')
#    print('getWavelengthMultiSpec: specDisp = ',specDisp)
    nDigitsBehindDot = []
    for dispVal in specDisp:
#        print(dispVal+'.find(".") = ',dispVal.find('.'))
        if dispVal.find('.') >= 0:
#            print('getWavelengthMultiSpec: dispVal = <'+dispVal+'>: found .')
            nDigitsBehindDot.append(len(dispVal) - dispVal.find('.') - 1)
        else:
#            print('getWavelengthMultiSpec: dispVal = <'+dispVal+'>: did not find .')
            if len(nDigitsBehindDot) == 0:
                nDigitsBehindDot.append(-2)
            else:
                nDigitsBehindDot.append(0-len(dispVal))
#    print('getWavelengthMultiSpec: nDigitsBehindDot = ',nDigitsBehindDot)

    startStr = 'spec'+str(iSpec+1)
    endStr = '"'
    specDispStr = dispStr[dispStr.find(startStr):]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDispStr = specDispStr[specDispStr.find(endStr)+1:]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDispStr = specDispStr[:specDispStr.find(endStr)]
#    print('getWavelengthMultiSpec: specDispStr = <'+specDispStr+'>')
    specDisp = specDispStr.split(' ')
#    print('getWavelengthMultiSpec: specDisp = ',specDisp)
    if len(specDisp) != len(nDigitsBehindDot):
        specDispNew = []
        iNDig = 0
        for iDispVal in range(len(specDisp)):
            dispVal = specDisp[iDispVal]
#            print('getWavelengthMultiSpec: iNDig = ',iNDig)
#            print('dispVal.find(".") = ',dispVal.find('.'))
#            print('nDigitsBehindDot[',iNDig,'] < 0 = ',nDigitsBehindDot[iNDig] < 0)
            if (dispVal.find('.') < 0) or ((dispVal.find('.') > 0) and (nDigitsBehindDot[iNDig] < 0)):
                if len(dispVal) > (0-nDigitsBehindDot[iNDig]):
                    specDispNew.append(dispVal[0:0-nDigitsBehindDot[iNDig]])
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
                    specDispNew.append(dispVal[0-nDigitsBehindDot[iNDig]:])
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
                    iNDig += 1
#                    print('getWavelengthMultiSpec: iNDig = ',iNDig)
                else:
                    specDispNew.append(dispVal)
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
            else:
#                print('getWavelengthMultiSpec: nDigitsBehindDot[',iNDig,'] = ',nDigitsBehindDot[iNDig])
                if nDigitsBehindDot[iNDig] < (len(dispVal) - dispVal.find('.') - 1):
#                    print('getWavelengthMultiSpec: wrong length: nDigitsBehindDot[iNDig] = ',nDigitsBehindDot[iNDig],', (len(dispVal)(=',len(dispVal),') - dispVal.find(".")(=',dispVal.find('.'),' - 1) = ',len(dispVal) - dispVal.find('.') - 1)
                    specDispNew.append(dispVal[0:dispVal.find('.')+1+nDigitsBehindDot[iNDig]])
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
                    specDispNew.append(dispVal[dispVal.find('.')+1+nDigitsBehindDot[iNDig]:])
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
                    iNDig += 1
#                    print('getWavelengthMultiSpec: iNDig = ',iNDig)
                else:
                    specDispNew.append(dispVal)
#                    print('getWavelengthMultiSpec: appended specDispNew[',len(specDispNew)-1,'] = ',specDispNew[len(specDispNew)-1])
            iNDig += 1
        specDisp = specDispNew
#    print('getWavelengthMultiSpec: specDisp = ',specDisp)
    nPix = int(header['NAXIS'+str(axis)])
    crVal = float(specDisp[3])
    cDelt = float(specDisp[4])
    lam = np.arange(crVal, crVal + (nPix*cDelt), cDelt, dtype=np.float32)
    if lam.shape[0] > nPix:
        lam = lam[0:nPix]
#    print 'getWavelength: lam = ',len(lam),': ',lam
    return lam

def getMaximumWavelengthRange(header, dispAxis=1):
    if dispAxis == 1:
        nSpec = header['NAXIS2']
    else:
        nSpec = header['NAXIS1']
#    print('getMaximumWavelengthRange: nSpec = ',nSpec)
    specRangeMax = [0.,100000.]
#    print('getMaximumWavelengthRange: specRangeMax = ',specRangeMax)
    for iSpec in range(nSpec):
        wavelength = getWavelengthMultiSpec(header, iSpec, axis=dispAxis)
        if wavelength[0] > specRangeMax[0]:
            specRangeMax[0] = wavelength[0]
#            print('getMaximumWavelengthRange: wavelength[0] = ',wavelength[0],' => specRangeMax[0] = ',specRangeMax[0])
        if wavelength[wavelength.shape[0]-1] < specRangeMax[1]:
            specRangeMax[1] = wavelength[wavelength.shape[0]-1]
#            print('getMaximumWavelengthRange: wavelength[',wavelength.shape[0]-1,'] = ',wavelength[wavelength.shape[0]-1],' => specRangeMax[1] = ',specRangeMax[1])
    print('getMaximumWavelengthRange: specRangeMax = ',specRangeMax)
    return specRangeMax

def calibratedFluxToAbsoluteFlux(calibratedFlux, distanceIn_pc):
    return calibratedFlux * distanceIn_pc * distanceIn_pc / 100.0

def rgb2grey(rgb):

    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray

def trace(image, searchRadius=5):
    axis = 0
    print('image.shape = ', image.shape)

    minPos = []#np.argmin(image,axis=1)
    startFound = False
    minIndex = 0
    for row in range(image.shape[axis]):
        addValue = True
        if not startFound:
            minIndex = np.argmin(image[row,:])
            print('row = ',row,': minIndex = ',minIndex,': image[',row,',',minIndex,'] = ',image[row,minIndex])
            if image[row,minIndex] < 10:
                startFound = True
                minPos.append(minIndex)
        else:
            searchStart = minPos[len(minPos)-1]-searchRadius
            if searchStart < 0:
                searchStart = 0

            searchEnd = minPos[len(minPos)-1]+searchRadius
            if searchEnd >= image.shape[axis]:
                searchEnd = image.shape[axis]-1
            print('row ',row,': searchStart = ',searchStart,', searchEnd = ',searchEnd,': ',image[row,searchStart:searchEnd])
            minIndex = np.argmin(image[row,searchStart:searchEnd])
            print('minIndex = ',minIndex)
            searchRadiusTemp = searchRadius
            iRun = 0
            while (minIndex == 0) and (image[row,searchStart + minIndex] == image[row,searchStart + minIndex + int(searchRadius/2)]):
                searchRadiusTemp = 2 * searchRadiusTemp
                searchStart = minPos[len(minPos)-1]-searchRadiusTemp
                if searchStart < 0:
                    searchStart = 0

                searchEnd = minPos[len(minPos)-1]+searchRadiusTemp
                if searchEnd >= image.shape[axis]:
                    searchEnd = image.shape[axis]-1
                print('row ',row,': searchStart = ',searchStart,', searchEnd = ',searchEnd,': ',image[row,searchStart:searchEnd])
                minIndex = np.argmin(image[row,searchStart:searchEnd])
                print('minIndex = ',minIndex)
                iRun += 1
                if iRun > 10:
                    addValue = False
                    break
            if addValue:
                minPos.append(searchStart + minIndex)
        if len(minPos) > 1:
            print('minPos[',len(minPos)-1,'] = ',minPos[len(minPos)-1])

    print('minPos = ',len(minPos),': ',minPos)
    minPos = np.array(minPos)
    minPos = image.shape[1] - np.flip(minPos,0)
    print('minPos = ',minPos.shape,': ',minPos)
    plt.plot(minPos)
    plt.show()
#    for row in range(image.shape[0]):
#        minPos.append()
#    print(image[:,100])
    return minPos

# gaiaData.ra, gaiaData.dec, ra, and dec are in decimal degrees
def getStarWithMinDist(gaiaData, ra, dec, iStar=0):
    dist = None
    index = None
    print('gaiaData.header = ',gaiaData.header)
    for i in range(gaiaData.size()):
#        print('gaiaData.getData(',i,') = ',gaiaData.getData(i))
#        print('gaiaData.getData(parallax,',i,') = ',gaiaData.getData('parallax',i))
        failed = False
        try:
            if gaiaData.getData('parallax',i) != '':
                parallax = float(gaiaData.getData('parallax',i))
                if parallax < 0.:
                    parallax = 0.1
                distance = Distance(parallax=parallax * u.mas)
    #            print('distance = ',distance)
    #            print("gaiaData.getData('pmra',i) = ",gaiaData.getData('pmra',i))
    #            print("gaiaData.getData('pmdec',i) = ",gaiaData.getData('pmdec',i))
                time = Time(float(gaiaData.getData('ref_epoch',i)), format='decimalyear')
    #            print('time = ',time)
                c = SkyCoord(ra=float(gaiaData.getData('ra',i))*u.degree,
                             dec=float(gaiaData.getData('dec',i))*u.degree,
                             distance=distance,
                             pm_ra_cosdec=float(gaiaData.getData('pmra',i)) * u.mas/u.yr,
                             pm_dec=float(gaiaData.getData('pmdec',i)) * u.mas/u.yr,
                             obstime=time)
    #            print('c = ',c)
                c_epoch2000 = c.apply_space_motion(Time('2000-01-01'))
    #            print('c_epoch2000 = ',c_epoch2000)
    #            print('c_epoch2000 = ',type(c_epoch2000),': ',dir(c_epoch2000),': ',c_epoch2000)
    #            print('c_epoch2000.ra = ',type(c_epoch2000.ra),': ',dir(c_epoch2000.ra),': ',c_epoch2000.ra)
    #            print('c_epoch2000.ra.deg = ',c_epoch2000.ra.deg)
                thisDist = angularDistance(ra, dec, c_epoch2000.ra.deg, c_epoch2000.dec.deg) * 3600.
        except:
            failed= True
        if failed:
            thisDist = angularDistance(ra, dec, float(gaiaData.getData('ra',i)), float(gaiaData.getData('dec',i))) * 3600.
        if (dist is None):
            dist = thisDist
            index = i
        else:
            if dist > thisDist:
                dist = thisDist
                index = i
                print('star ',iStar,': closest star index: ',index,': distance = ',dist)
        if (dist) < 1.:
            return [index, dist]
    return [index, dist]

def isLineEmpty(line):
    return len(line.strip()) == 0

def readFileToArr(fname):
    text_file = open(fname, "r")
    lines = text_file.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut

def lockAndWriteToFile(filename,flag,lockName,textToWrite):
    while os.path.isfile(lockName):
        print('lockName <'+lockName+'> exists, waiting 0.01 seconds')
        sleep(0.010)#ms
    open(lockName,'w').close()
    with open(filename, flag) as f:
        f.write(textToWrite)
    os.remove(lockName)

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

def getPNGName(lon,lat):
    png = '%010.6f' % lon
    png = png[:png.find('.')+2]
    #png = png.zfill(3)
    if lat > 0:
        png = png+'+'
    else:
        png = png+'-'
    if abs(lat) < 10.:
        png += '0'
    latStr = '%08.6f' % abs(lat)
    png = png + latStr
#    print('png = ',png)
    png = png[:png.rfind('.')+2]
#    print('lon = ',lon,', lat = ',lat,', png = <'+png+'>')
    return png


# @brief plot l and b every x degrees
def plotLBMarks(x):
    ham = Hammer()
    pixels = ham.getPixels()
    lArr = np.arange(0,360.1,0.1)
    bArr = np.arange(-90, 90.1, 0.1)
    xArr = []
    yArr = []
    for l in lArr:
        for b in np.arange(-90,91,x):
            xy = ham.lonLatToXY(l,b)
            xArr.append(xy.x)
            yArr.append(xy.y)
    for b in bArr:
        for l in np.arange(0,361,x):
            xy = ham.lonLatToXY(l,b)
            xArr.append(xy.x)
            yArr.append(xy.y)
        l=180.0001
        xy = ham.lonLatToXY(l,b)
        xArr.append(xy.x)
        yArr.append(xy.y)
    plt.scatter(xArr,yArr,s=0.1)
#    xy = ham.lonLatToXY(181.,0.)
#    plt.scatter([xy.x],[xy.y],s=100)

def applyVRadCorrection(wavelength, vRad):
    return np.array([wLen - (vRad * wLen / c0) for wLen in wavelength])

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[int(window_len/2):-int(window_len/2)]
