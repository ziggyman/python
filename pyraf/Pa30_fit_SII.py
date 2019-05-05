from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
import astropy.units as u
import numpy as np
#from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import re
from scipy.optimize import curve_fit
import subprocess

c0 = 299792.458
mean1 = 6716.44
mean2 = 6730.815
dMean = mean2 - mean1
sigma = 7.5 / 2.355
doPlot = False
SIIRatio = 2. / 3.

minRow = 1573
maxRow = 1613
minCol = 1000
maxCol = 1920
starCol = 1458
xCenter = starCol - minCol
yCenter = 994
sigmaVRad = 185.

maxDMean = 150.#50.
maxDSigma = 200.
maxRes = 3.e-37#0.8e-36#0.4e-36
minMean = 0.3e-19
minSum = 0.5e-21
maxMeanSigma = 0.1#1.e-1
maxSigSigma = 800.#1.e-1

colPlot = []
colRejecta = [104,110,127,279,313,468,490,529,536,543,557,559,613,629,651,652,656,857,864]
colRejectb = [22,53,348,362,416,418,606,613,621,819,822,824,825,856]

twoDSpecFile = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_av_x_wl_flt_cal_mSky_obs_not_smoothed_minComb.fits'
imFile = '/Users/azuri/daten/uni/HKU/Pa30/gtc_object_wcsIm_sky_0000956342.fits'


def removeBackground(hdulist, minRow, maxRow, minCol, maxCol, twoDSpec):
    upper = hdulist[0].data[maxRow:maxRow+20,minCol:maxCol]
    lower = hdulist[0].data[minRow-20:minRow,minCol:maxCol]
    for col in np.arange(0,twoDSpec.shape[1],1):
        twoDSpec[:,col] = twoDSpec[:,col] - np.mean([np.mean(upper[:,col]), np.mean(lower[:,col])])
    return twoDSpec

def getWavelength(header, axis=2):
    nPix = int(header['NAXIS'+str(axis)])
    print('getWavelength: nPix = ',nPix)
    crPix = int(header['CRPIX'+str(axis)])
    print('getWavelength: crPix = ',crPix)
    crVal = float(header['CRVAL'+str(axis)])
    print('getWavelength: crVal = ',crVal)
    cDelt = float(header['CDELT'+str(axis)])
    print('getWavelength: cDelt = ',cDelt)
  #    lam = np.ndarray(nPix, dtype=np.float32)
  #    lam[0] =
  #    for i in np.arange(1,nPix):
    lam = np.arange(crVal, crVal + ((nPix-0.5)*cDelt), cDelt, dtype=np.float32)
    print('getWavelength: lam = ',len(lam),': ',lam)
    return lam

hdulist = pyfits.open(twoDSpecFile)
header = hdulist[0].header
twoDSpec = hdulist[0].data[minRow:maxRow,minCol:maxCol]
twoDSpec = removeBackground(hdulist, minRow, maxRow, minCol, maxCol, twoDSpec)
wavelength = getWavelength(header)[minRow:maxRow]

def norm(x, mean, sd):
    norm = []
    for i in np.arange(0,x.size,1):
        norm += [1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x[i] - mean)**2/(2*sd**2))]
    return np.array(norm)

#x = np.linspace(-20, 20, 500)
#y_real = norm(x, mean1, std1) + norm(x, mean2, std2)

def res(p, y, x):
    m, dm, sd1, sd2 = p
    m1 = m
    m2 = m1 + dm
    y_fit = norm(x, m1, sd1) + norm(x, m2, sd2)
    err = y - y_fit
    return err

######################################
# Solving
#m, dm, sd1, sd2 = [5, 10, 1, 1]
#p = [m, dm, sd1, sd2] # Initial guesses for leastsq
#y_init = norm(x, m, sd1) + norm(x, m + dm, sd2) # For final comparison plot

#plsq = leastsq(res, p, args = (y_real, x))

#y_est = norm(x, plsq[0][0], plsq[0][2]) + norm(x, plsq[0][0] + plsq[0][1], plsq[0][3])


def gauss(x, mean, sigma, amp, y0=0.):
    return (amp * np.exp(0.-(0.5*(((x-mean)/(sigma))**2.)))) + y0

def getRaDecXY(string):
    strs = re.sub( '\s+', ' ', string ).strip()
    #  print('getRaDecXY: strs = ',strs)
    strs = strs.rstrip().split(' ')
    #print('getRaDecXY: strs = ',strs)
    return [strs[0], strs[1], float(strs[len(strs)-2]), float(strs[len(strs)-1])]

def getArcsecDistance(fitsName, x1, y1, x2, y2):
    result1 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x1), str(y1)])
    result2 = subprocess.check_output(['xy2sky', '-j', fitsName, str(x2), str(y2)])
    #print('xy2sky(',str(x1),', ',str(y1),') = ',result1)
    #print('xy2sky(',str(x2),', ',str(y2),') = ',result2)
    raHMS1, decHMS1, x1, y1 = getRaDecXY(result1.decode('utf-8'))
    raHMS2, decHMS2, x2, y2 = getRaDecXY(result2.decode('utf-8'))
    #print('raHMS1 = ',raHMS1,', decHMS1 = ',decHMS1)
    #print('raHMS2 = ',raHMS2,', decHMS2 = ',decHMS2)
    mm1 = SkyCoord(ra=raHMS1, dec=decHMS1, unit=(u.hourangle, u.deg))
    mm2 = SkyCoord(ra=raHMS2, dec=decHMS2, unit=(u.hourangle, u.deg))
    return mm1.separation(mm2).arcsecond

def lambdaToY(lam, lambdaArr):
#  return len(lambdaArr) - ((lam-lambdaArr[0]) * len(lambdaArr) / (lambdaArr[len(lambdaArr)-1] - lambdaArr[0]))# + 0.5
  return ((lam-lambdaArr[0]) * len(lambdaArr) / (lambdaArr[len(lambdaArr)-1] - lambdaArr[0])) - 0.5

# pixRangeMove[0] < 0
# pixRangeMove[1] > 0
#def crossCorrelate(yStatic, yMove, pixRangeMove):
#    print('yStatic.shape = ',yStatic.shape)
#    print('yMove.shape = ',yMove.shape)
#    print('pixRangeMove = ',pixRangeMove)
#    yStaticTemp = yStatic[pixRangeMove[1]:len(yStatic) + pixRangeMove[0]]
#    print('yStaticTemp.shape = ',yStaticTemp.shape)
#    res = []
#    for shift in np.arange(pixRangeMove[0], pixRangeMove[1]+1, 1):
#        x0 = pixRangeMove[1]-shift
#        x1 = len(yMove) + pixRangeMove[0] - shift
#        yMoveTemp = yMove[x0:x1]
#        print('shift = ',shift,': x0 = ',x0,', x1 = ',x1,': yMoveTemp.shape = ',yMoveTemp.shape)
#        diffSquare = (yStaticTemp - yMoveTemp)**2
#        print('shift = ',shift,': diffSquare = ',diffSquare)
#        res.append(np.sum(diffSquare))
#    print('res = ',res)
#    minimum = min(res)
#    maximum = max(res)
#    minPos = res.index(minimum)
#    print('minPos = ',minPos,', minimum = ',minimum)
#    popt,pcov = curve_fit(gauss,range(len(res)),res,p0=[minimum - maximum,minPos,0.5,maximum])
#    print('popt = ',popt,', pcov = ',pcov)
#    plt.plot(gauss(range(len(res)), popt[0], popt[1], popt[2], popt[3]), 'g-')
#    plt.show()
#    return [minPos,minimum]

def run():
    rows = np.arange(0,twoDSpec.shape[0],1)
    nRows = len(rows)
    cols = np.arange(0,twoDSpec.shape[1],1)
    nCols = len(cols)

    goodVRads = []
    guesses = []
    for col in np.arange(0,nCols,1):
        print(' ')
        print('col = ',col)
        doPlot = False
        if col in colPlot:
            doPlot = True
        if doPlot:
            plt.plot(wavelength, twoDSpec[:,col],'g-')
            plt.show()
        # for negative and positive radial velocities (near side and far side)
        iVRad = 0
        for vradRange in [[-1200.,300.],[300.,1200]]:
            if ((iVRad == 0) and (col not in colRejecta)) or ((iVRad > 0.) and (col not in colRejectb)):
    #        print('iVRad = ',iVRad,': vradRange = ',vradRange)
                res = []
                vradArr = []
                specs = []
                synSpecs = []
                wavelengths = []
                ranges = []
                # for each radial velocity in range in steps of 5
                for vrad in np.arange(vradRange[0],vradRange[1],5.):
        #            print('vrad = ',vrad)
                    if vrad < 0.:
                        valRange = [0,int(twoDSpec.shape[0]*4/5)]
                    else:
                        valRange = [int(twoDSpec.shape[0]*1/5),twoDSpec.shape[0]]

    #                print('valRange = ',valRange)
                    ranges.append(valRange)
                    vradArr.append(vrad)
                    wavelengths.append(wavelength[valRange[0]:valRange[1]])

                    synSpec = np.zeros(valRange[1]-valRange[0])
                    iRun = 0
                    # add peak for both [SII] species
                    amp = 0.
                    for mean in [((vrad * mean1) / c0) + mean1, ((vrad * mean2) / c0) + mean2]:
        #                print('mean = ',mean)
                        amp = np.max(twoDSpec[valRange[0]:valRange[1],col])
                        if iRun == 1:
                            amp = amp * SIIRatio
        #                print('amp = ',amp)
        #                print('sigma = ',sigma)
                        synSpec += gauss(wavelengths[len(wavelengths)-1], mean, sigma, amp)
                        iRun += 1
        #            print('synSpec = ',synSpec)
                    synSpecs.append(synSpec)
            #        plt.plot(wavelength, synSpec,'b-')
                    specs.append(twoDSpec[valRange[0]:valRange[1],col])

                    # calculate redisual
                    res.append(np.sum((synSpec - specs[len(specs)-1]) ** 2.))

    #        print('vradArr = ',len(vradArr),': ',vradArr,': res = ',len(res),': ',res)

                # determine minimum residuum and position
                minimum = min(res)
                #maximum = max(res)
                minPos = res.index(minimum)
                print('valRange = ',valRange)
                print('minPos = ',minPos,', minimum = ',minimum)
                print('len(res) = ',len(res))
    #            print('mean(spec) = ',np.mean(specs[minPos]))
                print('mean(specs[minPos]) = ',np.mean(specs[minPos]))
                print('sum(specs[minPos]) = ',np.sum(specs[minPos]))

                if (minPos > 3) and (minPos < (len(res)-4)):
                    # plot spectrum and best fitting synthetic spectrum
                    if doPlot:
                        plt.plot(wavelengths[minPos],specs[minPos],'b-')
                        plt.plot(wavelengths[minPos],synSpecs[minPos],'g-')
                        plt.xlabel('wavelength')
                        plt.ylabel('spectrum')
                        plt.show()

                    # fit Gaussian to residuum curve as function of radial velocity
                    # for continuum level take minimum of maxima on boths sides of best vrad
                    y0 = min([max(res[0:minPos]), max(res[minPos:])])
        #            print('y0 = ',y0)

                    # only fit minPos +/- 20 radial velocities
                    x1 = minPos-20
                    if x1 < 0:
                        x1 = 0
                    x2 = minPos+20
                    if x2 >= len(vradArr):
                        x2 = len(vradArr)-1
        #            print('x1 = ',x1,', x2 = ',x2)
        #            print('vradArr[x1 = ',x1,'] = ',vradArr[x1],', vradArr[x2 = ',x2,'] = ',vradArr[x2])
                    xFit = vradArr[x1:x2]
        #            print('xFit = ',xFit)
                    yFit = res[x1:x2]
        #            print('yFit = ',yFit)

                    goodFit = False
                    try:
                        # fit Gaussian
                        guess = [vradArr[minPos], sigmaVRad, minimum-y0, y0]
                        print('guess = [mean=',guess[0],', sigma=',guess[1],', amplitude=',guess[2],', y0=',guess[3],']')
                        popt,pcov = curve_fit(gauss,xFit,yFit,p0=guess)
                        print('popt = ',popt,', pcov = ',pcov)
                        if doPlot:
                            plt.plot(vradArr, gauss(vradArr, popt[0], popt[1], popt[2], popt[3]), 'g-')
                        # fit was successful
                        goodFit = True
                    except:
                        print('FAILED TO FIT GAUSSIAN')
    #                    if col in colPlot:
    #                        STOP

                    # check if fit was good
                    if (goodFit
                        and (abs(popt[0] - vradArr[minPos]) < maxDMean)
                        and (abs(popt[1] - sigmaVRad) < maxDSigma)
                        and (popt[2] < 0.)
                        and (minimum < maxRes)
                        and (np.mean(specs[minPos]) > minMean)
                        and (np.sum(specs[minPos]) > minSum)
                        and (pcov[0,0] < maxMeanSigma)
                        and (pcov[1,1] < maxSigSigma)):
                        print('minMean = ',minMean,', mean(specs[',minPos,']) = ',np.mean(specs[minPos]))

                        dist = getArcsecDistance(imFile, starCol, yCenter, minCol+col, yCenter)
                        # if spectrum left of star
                        if col < (starCol - minCol):
                            dist = 0. - dist
                        goodVRads.append([dist, popt, pcov, guess, col, res, minPos, specs, synSpecs])
                        if doPlot:
                            plt.plot(vradArr,res,'b+')
                            plt.xlabel('v_rad')
                            plt.ylabel('residual')
                        print('Gaussian considered good')
                    else:
                        if not goodFit:
                            print('fit failed')
                        else:
                            if abs(popt[0] - vradArr[minPos]) >= maxDMean:
                                print('abs(popt[0] - vradArr[minPos])=',abs(popt[0] - vradArr[minPos]),' >= maxDMean=',maxDMean)
                            if abs(popt[1] - sigmaVRad) >= maxDSigma:
                                print('abs(popt[1] - sigmaVRad)=',abs(popt[1] - sigmaVRad),' >= maxDSigma=',maxDSigma)
                            if popt[2] >= 0.:
                                print('popt[2]=',popt[2],' >= 0.')
                            if minimum >= maxRes:
                                print('minimum=',minimum,' >= maxRes=',maxRes)
                            if np.mean(specs[minPos]) <= minMean:
                                print('np.mean(specs[minPos])=',np.mean(specs[minPos]),' <= minMean=',minMean)
                            if np.sum(specs[minPos]) <= minSum:
                                print('np.sum(specs[minPos])=',np.sum(specs[minPos]),' <= minSum=',minSum)
                            if pcov[0,0] >= maxMeanSigma:
                                print('pcov[0,0]=',pcov[0,0],' >= maxMeanSigma=',maxMeanSigma)
                            if pcov[1,1] >= maxSigSigma:
                                print('pcov[1,1]=',pcov[1,1],' >= maxSigSigma=',maxSigSigma)
                        print('Gaussian fit not good')
                    if doPlot:
                        plt.show()

                #        minimum, minPos = crossCorrelate(twoDSpec[0:int(twoDSpec.shape[0]*2/3),col], synSpec[0:int(twoDSpec.shape[0]*2/3)], [-,0])
                #        minima.append(minimum)
                #        minPositions.append(minPos)
                #        print('minima = ',minima)
                #        print('minPositions = ',minPositions)
                else:
                    print('minPos outside range')

                    if doPlot:
                        plt.plot(wavelengths[minPos],specs[minPos],'b-')
                        plt.show()

                        plt.plot(vradArr,res,'b+')
                        plt.xlabel('v_rad1')
                        plt.ylabel('residual1')
                        plt.show()

                if False:#col == 1390-minCol:
                    plt.plot(wavelengths[len(wavelengths)-1], specs[len(specs)-1],'g-')
                    plt.plot(wavelengths[len(wavelengths)-1], synSpecs[len(synSpecs)-1],'b-')
                    plt.show()
                    distPlot = [x[0] for x in goodVRads]
                    vradPlot = [x[1] for x in goodVRads]
                    plt.plot(distPlot,vradPlot,'g+')
                    plt.show()
                    STOP
            iVRad += 1
    print('goodVRads = ',len(goodVRads),': ',goodVRads)

    if False:
        distPlot = [x[4] for x in goodVRads]
        vradPlot = [x[1][0] for x in goodVRads]
        plt.plot(distPlot,vradPlot,'g+')
        plt.xlabel('column')#center distance [arcsec]')
        plt.ylabel('radial velocity [km/s]')
        plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_vrad_SII_fit_maxDM=%d_maxSig=%d.eps' % (maxDMean, maxDSigma), format='eps', frameon=False, bbox_inches='tight', pad_inches=0.1)
        plt.show()

    return goodVRads


# find x in Rows from [xFrom, yFrom] for column yTo, do it the crude way...
def findXForArcSecDistanceFrom(xFrom, yFrom, yTo, dist, xRange, fitsName):
  #print('findXForArcSecDistanceFrom: xFrom = ',xFrom,', yFrom =',yFrom,', yTo = ',yTo,', dist = ',dist,', xRange = ',xRange,', fitsName = ',fitsName)
  xRangeTemp = xRange
  distTemp = dist
  if dist < 0:
    distTemp = 0. - dist
    xRangeTemp = [xRange[0], xFrom]
  else:
    if xFrom > xRange[0]:
      xRangeTemp = [xFrom, xRange[1]]

#print('findXForArcSecDistanceFrom: dist = ',dist,', distTemp = ',distTemp)
  hdulist = pyfits.open(fitsName)
  nCols = hdulist[0].data.shape[1]
#print('findXForArcSecDistanceFrom: nCols = ',nCols)
  previousDist = 100000000.
  for x in np.arange(xRangeTemp[0], xRangeTemp[1]):
    #   print('findXForArcSecDistanceFrom: xFrom = ',xFrom,', yFrom =',yFrom,', yTo = ',yTo,', distTemp = ',distTemp,', xRange = ',xRange,', fitsName = ',fitsName)
    tempDist = getArcsecDistance(fitsName, xFrom, yFrom, x, yTo)
    #print('findXForArcSecDistanceFrom): distTemp = ',distTemp,': x = ',x,': tempDist = ',tempDist)
    if (((dist > 0.) and ((previousDist < distTemp) and (tempDist >= distTemp)))
        or ((dist < 0.) and ((previousDist > distTemp) and (tempDist <= distTemp)))
        or ((distTemp == 0) and (tempDist == 0.))):
      #print('findXForArcSecDistanceFrom: returning x = ',x)
      return x

    previousDist = tempDist
  raise('findXForArcSecDistanceFrom: ERROR: could not determine x')

def findYForWLen(wLenArr, wLen):
    print('findYForWLen: wLen = ',wLen,', wLenArr = ',wLenArr)
    previousWLen = 10000000.
    pos = 0
    for w in wLenArr:
        #    print('findYForWLen: previousWLen = ',previousWLen,', w = ',w,', wLen = ',wLen)
        if (previousWLen < wLen) and (w >= wLen):
          #  print('findYForWLen: returning pos = ',pos)
          return pos
        pos += 1
        previousWLen = w
    raise('findYForWLen: ERROR: wLen ',wLen,' not found')

def lambda0AndVradToLambda(lambda0, vrad):
    return lambda0 + (vrad * lambda0 / c0) # vrad = c * (lambda0 - lambda) / lambda0 -> lambda = lambda0 - (vrad * lambda0 / c)

def plotGoodVRads():
    goodVRads = run()#                    goodVRads.append([dist, popt, pcov, guess, col, res, minPos, specs, synSpecs])
                     #                    guess = [vradArr[minPos], sigmaVRad, minimum-y0, y0]

    rows = range(twoDSpec.shape[0])
    cols = range(twoDSpec.shape[1])

    xlim = [-0.5,len(cols)-0.5]
    ylim = [-0.5,len(rows)-0.5]

    # show 2D spectrum and reverse y (because of how matrixes are stored / interpreted)
    twoDSpecPlot = np.ndarray(twoDSpec.shape)
    for row in range(twoDSpec.shape[0]):
        twoDSpecPlot[row,:] = twoDSpec[twoDSpec.shape[0]-1-row,:]
    plt.imshow(twoDSpecPlot, cmap='Greys',vmin=0., vmax=2.0e-19, extent=(xlim[0],xlim[1],ylim[0],ylim[1]), aspect='auto')

    # mark CS with X
    plt.plot([len(cols)/2.],[len(rows)/2.],'rx',markersize=10)
    limits = plt.axis()

    # plot measured lines with vrad as color
    cm = plt.cm.get_cmap('rainbow')

    vrads = [x[1][0] for x in goodVRads]

    positiveVRads = []
    positiveIndices = []
    negativeVRads = []
    negativeIndices = []
    for i in np.arange(0,len(vrads),1):
        if vrads[i] < 0:
            negativeVRads.append(vrads[i])
            negativeIndices.append(i)
        else:
            positiveVRads.append(vrads[i])
            positiveIndices.append(i)
    print('positiveVRads = ',len(positiveVRads),': ',positiveVRads)
    print('positiveIndices = ',len(positiveIndices),': ',positiveIndices)
    print('negativeVRads = ',len(negativeVRads),': ',negativeVRads)
    print('negativeIndices = ',len(negativeIndices),': ',negativeIndices)

    SII6716a = [lambda0AndVradToLambda(mean1, vrad) for vrad in negativeVRads]
    print('mean1 = ',mean1,', negativeVRads[0] = ',negativeVRads[0],': SII6716a[0] = ',SII6716a[0])
    SII6716b = [lambda0AndVradToLambda(mean1, vrad) for vrad in positiveVRads]
    SII6731a = [lambda0AndVradToLambda(mean2, vrad) for vrad in negativeVRads]
    SII6731b = [lambda0AndVradToLambda(mean2, vrad) for vrad in positiveVRads]
    print('SII6716a = ',len(SII6716a),': ',SII6716a)
    print('SII6716b = ',len(SII6716b),': ',SII6716b)
    print('SII6731a = ',len(SII6731a),': ',SII6731a)
    print('SII6731b = ',len(SII6731b),': ',SII6731b)
    vradSII6716a = negativeVRads
    vradSII6716b = positiveVRads
    vradSII6731a = negativeVRads
    vradSII6731b = positiveVRads
    print('vradSII6716a = ',len(vradSII6716a),': ',vradSII6716a)
    print('vradSII6716b = ',len(vradSII6716b),': ',vradSII6716b)
    print('vradSII6731a = ',len(vradSII6731a),': ',vradSII6731a)
    print('vradSII6731b = ',len(vradSII6731b),': ',vradSII6731b)

    for colorCode in ['vrad','res','mu','sum','meanSigma','sigSigma']:
        i = 0
        sc = None
        for lam in [SII6716a,SII6716b,SII6731a,SII6731b]:
            plotVRadY = lambdaToY(lam, wavelength)
            print('lam = ',lam)
            print('lambdaToY(lam, wavelength) = ',lambdaToY(lam, wavelength))
            #plotVRadY = len(wavelength) - lambdaToY(lam, wavelength)
            print('plotVRadY = ',len(plotVRadY),': ',plotVRadY)
            plotVRadX = None
            plotRes = None
            plotMean = None
            plotSum = None
            plotMeanSigma = None
            plotSigSigma = None
              #    if i < 2:
              #  marker = 'r+'
              #else:
              #  marker = 'b+'
                  #    plt.plot(plotVRadX,
                  #   plotVRadY,
                  #   marker)
            marker = None
            label = None
            if i == 0:
                vrad = vradSII6716a
                plotVRadX = np.array([goodVRads[a][4] for a in negativeIndices])
                plotRes = np.array([goodVRads[a][5][goodVRads[a][6]] for a in negativeIndices])
                plotMean = np.array([np.mean(goodVRads[a][7][goodVRads[a][6]]) for a in negativeIndices])
                plotSum = np.array([np.sum(goodVRads[a][7][goodVRads[a][6]]) for a in negativeIndices])
                plotMeanSigma = np.array([goodVRads[a][2][0,0] for a in negativeIndices])
                plotSigSigma = np.array([goodVRads[a][2][1,1] for a in negativeIndices])
                marker = 'o'
                label = '[SII] 6716'
            elif i == 1:
                vrad = vradSII6716b
                plotVRadX = np.array([goodVRads[a][4] for a in positiveIndices])
                plotRes = np.array([goodVRads[a][5][goodVRads[a][6]] for a in positiveIndices])
                plotMean = np.array([np.mean(goodVRads[a][7][goodVRads[a][6]]) for a in positiveIndices])
                plotSum = np.array([np.sum(goodVRads[a][7][goodVRads[a][6]]) for a in positiveIndices])
                plotMeanSigma = np.array([goodVRads[a][2][0,0] for a in positiveIndices])
                plotSigSigma = np.array([goodVRads[a][2][1,1] for a in positiveIndices])
                marker = 'o'
#                label = '[SII] 6716b'
            elif i == 2:
                vrad = vradSII6731a
                plotVRadX = np.array([goodVRads[a][4] for a in negativeIndices])
                plotRes = np.array([goodVRads[a][5][goodVRads[a][6]] for a in negativeIndices])
                plotMean = np.array([np.mean(goodVRads[a][7][goodVRads[a][6]]) for a in negativeIndices])
                plotSum = np.array([np.sum(goodVRads[a][7][goodVRads[a][6]]) for a in negativeIndices])
                plotMeanSigma = np.array([goodVRads[a][2][0,0] for a in negativeIndices])
                plotSigSigma = np.array([goodVRads[a][2][1,1] for a in negativeIndices])
                marker = '^'
                label = '[SII] 6731'
            elif i == 3:
                vrad = vradSII6731b
                plotVRadX = np.array([goodVRads[a][4] for a in positiveIndices])
                plotRes = np.array([goodVRads[a][5][goodVRads[a][6]] for a in positiveIndices])
                plotMean = np.array([np.mean(goodVRads[a][7][goodVRads[a][6]]) for a in positiveIndices])
                plotSum = np.array([np.sum(goodVRads[a][7][goodVRads[a][6]]) for a in positiveIndices])
                plotMeanSigma = np.array([goodVRads[a][2][0,0] for a in positiveIndices])
                plotSigSigma = np.array([goodVRads[a][2][1,1] for a in positiveIndices])
                marker = '^'
#                label = '[SII] 6731b'
            else:
                print('i = ',i,' outside range')
                STOP
        #    label = '[SII] 6731'
            print('plotVRadX = ',plotVRadX)
            print('plotVRadX = ',len(plotVRadX),': ',plotVRadX)

            if colorCode == 'vrad':
                cPlot = vrad
                colLabel = 'radial velocity [km/s]'
                vrange=[-1200.,1200.]
            elif colorCode == 'res':
                cPlot = plotRes
                colLabel = 'residuals'
                residuals = np.array([goodVRads[a][5][goodVRads[a][6]] for a in np.arange(0,len(goodVRads))])
                vrange=[min(residuals),max(residuals)]
            elif colorCode == 'mu':
                cPlot = plotMean
                colLabel = 'mean of spectrum'
                mus = np.array([np.mean(goodVRads[a][7][goodVRads[a][6]]) for a in np.arange(0,len(goodVRads))])
                vrange=[min(mus),max(mus)]
            elif colorCode == 'sum':
                cPlot = plotSum
                colLabel = 'sum of spectrum'
                mus = np.array([np.sum(goodVRads[a][7][goodVRads[a][6]]) for a in np.arange(0,len(goodVRads))])
                vrange=[min(mus),max(mus)]
            elif colorCode == 'meanSigma':
                cPlot = plotMeanSigma
                colLabel = 'sigma(mu)'
                sigs = np.array([goodVRads[a][2][0,0] for a in np.arange(0,len(goodVRads))])
                vrange=[min(sigs),max(sigs)]
            elif colorCode == 'sigSigma':
                cPlot = plotSigSigma
                colLabel = 'sigma(sigma)'
                sigs = np.array([goodVRads[a][2][1,1] for a in np.arange(0,len(goodVRads))])
                vrange=[min(sigs),max(sigs)]

            for j in np.arange(0,len(plotVRadX),1):
                print('x = ',plotVRadX[j],', y = ',plotVRadY[j],': value = ',cPlot[j])

            sc = plt.scatter(plotVRadX,
                             plotVRadY,
                             c=cPlot,
                             vmin=vrange[0],#-1200.,
                             vmax=vrange[1],#1200.,
                             s=15,
                             cmap=cm,
                             marker = marker,
                             label = label
                            )
            i += 1
        plt.xlim = xlim
        plt.ylim = ylim
        plt.axis(limits)
        plt.autoscale(False)
        plt.colorbar(sc,label=colLabel)
        plt.legend()

        if colorCode == 'vrad':
            # change x-Axis
            xTicksDist = [-100., -80., -60., -40., -20., 0., 20., 40., 60.,80.,100.]
            xTicksCols = []
            xRange = [minCol, maxCol]
            for xTick in xTicksDist:
              if xTick != xTicksDist[0]:
                xRange = [xTicksCols[len(xTicksCols)-1], xRange[1]]
              xTicksCols.append(findXForArcSecDistanceFrom(xCenter+minCol, yCenter, yCenter, xTick, xRange, imFile))
            print('xTicksCols = ',xTicksCols)

            xTicks = []
            for x in xTicksCols:
              xTicks.append(x - minCol)
            print('xTicks = ',xTicks)

            xTicksDistStr =()
            for x in xTicksDist:
              xTicksDistStr = xTicksDistStr + ('%.0d' % (x),)
            print('findXForArcSecDistanceFrom: xTicksDistStr = ',xTicksDistStr)

            locks, labels = plt.xticks()
            print('findXForArcSecDistanceFrom: locks = ',locks,', labels[0] = ',labels[0])

            plt.xticks(xTicks,
                       xTicksDistStr)

            plt.xlabel('center distance [arcsec]')

        if colorCode == 'vrad':
            # change y axis
            yTicksWLen = ('%.0d' % ((int(getWavelength(header)[minRow] / 10.) + 1.) * 10.),)
            #yTicksRow = [findYForWLen(wavelength, int(yTicksWLen[len(yTicksWLen)-1]))]
            yTicksRow = [lambdaToY(int(yTicksWLen[len(yTicksWLen)-1]), wavelength)]
            while True:
              if int(yTicksWLen[len(yTicksWLen)-1]) + 10. > getWavelength(header)[maxRow]:
                break
              yTicksWLen = yTicksWLen + ('%.0d' % (int(yTicksWLen[len(yTicksWLen)-1]) + 10.),)
#              yTicksRow.append(findYForWLen(wavelength, int(yTicksWLen[len(yTicksWLen)-1])))
              yTicksRow.append(lambdaToY(int(yTicksWLen[len(yTicksWLen)-1]),wavelength))
            print('findXForArcSecDistanceFrom: yTicksRow = ',yTicksRow)
            print('findXForArcSecDistanceFrom: yTicksWLen = ',yTicksWLen)
            plt.yticks(yTicksRow, yTicksWLen)
            plt.ylabel('wavelength [$\mathrm{\AA}$]')

        # insert magnification of most striking features
        if False:
            a = plt.axes([.135, .67, .15, .18])
            magMinRow = 4
            magMaxRow = twoDSpecPlot.shape[0] - 4
            magMinCol = 300
            magMaxCol = 400
            plt.imshow(twoDSpecPlot[magMinRow:magMaxRow,magMinCol:magMaxCol], cmap='Greys',vmin=0., vmax=2.0e-19, extent=(xlim[0],xlim[1],ylim[0],ylim[1]), aspect='auto')
            #plt.plot([1.,2.], [1.,2.], 'b+')
            #plt.title('Probability')
            plt.xticks([])
            plt.yticks([])


        plt.savefig('/Users/azuri/daten/uni/HKU/Pa30/report/images/pa30_'+colorCode+'_map_on_2dspec_new_maxDMean=%d_maxDSigma=%d_maxRes='
                    % (maxDMean, maxDSigma) + str(maxRes)+'_minMean='+str(minMean)+'_minSum='+str(minSum)+'_maxMeanSigma=%.3f_maxSigSigma=%.3f' % (maxMeanSigma, maxSigSigma) +'.eps',
                    format='eps',
                    frameon=False,
                    bbox_inches='tight',
                    pad_inches=0.1)

        plt.show()


plotGoodVRads()