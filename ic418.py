import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from drUtils import getImageData
from myUtils import angularDistanceFromXYPyAsl,hmsToDeg,dmsToDeg,angularDistancePyAsl

centerPos = [1763.47,1545.43]#[ratioMap.shape[0]/2.,ratioMap.shape[1]/2.]
xC = centerPos[0]
yC = centerPos[1]

""" check if a pixel with center coordinates [xC,yC] is in between the 2 lines each defined by 2 points [xyi0[0],xyi0[1]] and [xyi1[0],xyi1[1]]"""
def pixelCenterIsInArea(xC,yC,xy00,xy01,xy10,xy11):
    if (xy00[1] == xy01[1]):
        return ((yC > xy00[1]) and (yC < xy10[1])) or ((yC < xy00[1]) and (yC > xy10[1]))
    if (xy00[0] == xy01[0]):
        return ((xC > xy00[0]) and (xC < xy10[0])) or ((xC < xy00[0]) and (xC > xy10[0]))
    m0 = (xy01[1]-xy00[1]) / (xy01[0]-xy00[0])
    n0 = xy00[1] - (m0*xy00[0])
    m1 = (xy11[1]-xy10[1]) / (xy11[0]-xy10[0])
    n1 = xy10[1] - (m1*xy10[0])
    y0 = m0 * xC + n0
    y1 = m1 * xC + n1
    if ((yC - y0 > 0.) and (yC - y1 < 0.)) or ((yC - y0 < 0.) and (yC - y1 > 0.)):
        return True

def integrateBetween(xyL,xy00,xy01,xy10,xy11):
    fig,ax = plt.subplots(1)
    ax.plot(xyL.shape[0],xyL.shape[0],marker='v',color='white')
    ax.imshow(np.transpose(xyL),origin='lower')
    ax.plot([xy00[0],xy01[0]],[xy00[1],xy01[1]],'r-')
    ax.plot([xy10[0],xy11[0]],[xy10[1],xy11[1]],'r-')
    left, right = plt.xlim()
    if (left < 0.) or (right > 2.*xC):
        ax.set_xlim([0.,2. * xC])
#    plt.show()
#    return [0.0,1]
    squares = []
    integral = 0.
    nSquares = 0
    for ix in range(xyL.shape[0]):
        print('checking pixel at [ix=',ix,']')
        for iy in range(xyL.shape[1]):
            if pixelCenterIsInArea(ix+0.5,iy+0.5,xy00,xy01,xy10,xy11):
                #print('[',ix,',',iy,'] is inside')
#                STOP
                if xyL[ix,iy] > 0:
                    integral += xyL[ix,iy]
                    squares.append(Rectangle((ix,iy),1,1))
                    nSquares += 1
#            if ix > xy00[0] and ix < xy10[0]:
#                STOP
    pc = PatchCollection(squares,facecolor='r',alpha=0.9,edgecolor='b')
    ax.add_collection(pc)
    plt.show()
    return [integral / nSquares,nSquares]


def integrateBetween(xyLOIII,xyLHbeta,xy00,xy01,xy10,xy11):
    fig,ax = plt.subplots(1)
    ax.plot(xyLOIII.shape[0],xyLOIII.shape[0],marker='v',color='white')
    ax.imshow(np.transpose(xyLOIII),origin='lower')
    ax.plot([xy00[0],xy01[0]],[xy00[1],xy01[1]],'r-')
    ax.plot([xy10[0],xy11[0]],[xy10[1],xy11[1]],'r-')
    left, right = plt.xlim()
    if (left < 0.) or (right > 2.*xC):
        ax.set_xlim([0.,2. * xC])
#    plt.show()
#    return [0.0,1]
    squares = []
    integralOIII = 0.
    integralHbeta = 0.
    nSquares = 0
    for ix in range(xyLOIII.shape[0]):
        print('checking pixel at [ix=',ix,']')
        for iy in range(xyLOIII.shape[1]):
            if pixelCenterIsInArea(ix+0.5,iy+0.5,xy00,xy01,xy10,xy11):
                #print('[',ix,',',iy,'] is inside')
#                STOP
                if xyLOIII[ix,iy] > 0:
                    integralOIII += xyLOIII[ix,iy]
                    squares.append(Rectangle((ix,iy),1,1))
                    nSquares += 1
                if xyLHbeta[ix,iy] > 0:
                    integralHbeta += xyLHbeta[ix,iy]
#            if ix > xy00[0] and ix < xy10[0]:
#                STOP
    pc = PatchCollection(squares,facecolor='r',alpha=0.9,edgecolor='b')
    ax.add_collection(pc)
    plt.show()
    return [integralOIII,integralHbeta,nSquares]

def getAverageRatioFromRatioImage():

    fitsName = '/Users/azuri/daten/uni/HKU/IC418/IC418_bluebeam_ratiomap2_without_cs.fits'#IC418_4861_sum2_without_cs_8_v2.fits'
    ratioMap = getImageData(fitsName,0)

    dist1Pix = angularDistancePyAsl(hmsToDeg('5:27:28.421'),dmsToDeg('-12:41:59.03'),hmsToDeg('5:27:28.3735'),dmsToDeg('-12:41:59.029')) * 3600.
    print('dist1Pix = ',dist1Pix)

    dist = angularDistanceFromXYPyAsl(fitsName,centerPos[0],centerPos[1],centerPos[0]+1.,centerPos[1])
    print('dist = ',dist)

    slitWidthInArcSec = 4.
    slitWidthInPixels = slitWidthInArcSec * 100. / dist
    print('slitWidthInPixels = ',slitWidthInPixels)

    print('centerPos = ',centerPos)

    angles = np.arange(0,181,15)
    print('angles = ',angles)

    plt.plot(ratioMap.shape[1],ratioMap.shape[0],marker='v',color='white')
    plt.imshow(ratioMap,origin='lower')
    plt.show()

    xl = np.arange(0,ratioMap.shape[1]*100,1)
    yl = np.arange(0,ratioMap.shape[0]*100,1)
    xyL = np.zeros((xl.shape[0],yl.shape[0]))
    print('xyL.shape = ',xyL.shape)
    for ix in range(xl.shape[0]):
        for iy in range(yl.shape[0]):
            xyL[ix,iy] = ratioMap[int(iy/100),int(ix/100)]
    plt.imshow(np.transpose(xyL),origin='lower')
    plt.show()

    if False:
        xy00 = [centerPos[0] - (slitWidthInPixels / 2.),0]
        print('xy00 = ',xy00)
        xy01 = [centerPos[0] - (slitWidthInPixels / 2.),centerPos[1]*2.]
        print('xy01 = ',xy01)

        xy10 = [centerPos[0] + (slitWidthInPixels / 2.),0]
        print('xy10 = ',xy10)
        xy11 = [centerPos[0] + (slitWidthInPixels / 2.),centerPos[1]*2.]
        print('xy11 = ',xy11)

        integral = integrateBetween(xyL,xy00,xy01,xy10,xy11)
        print('integral = ',integral)
        #STOP

    resultsFileName = '/Users/azuri/daten/uni/HKU/IC418/results.csv'
    results = []
    with open(resultsFileName,'w') as f:
        f.write('angle,averge ratio\n')
        for angle in angles:
            ang = np.radians(angle)
            print('angle = ',angle,', ang = ',ang)
            if angle != 90.:

                r = yC / (np.cos(ang) + 1e-10)
                print('r = ',r)

                xx = xC + (r * np.sin(ang))
                print('xx = ',xx)
                xy00 = [xx - (slitWidthInPixels / (2. * np.cos(ang))),0.]
                print('xy00 = ',xy00)
                xy10 = [xx + (slitWidthInPixels / (2. * np.cos(ang))),0.]
                print('xy10 = ',xy10)

                xy = xC - (r * np.sin(ang))
                print('xy = ',xy)
                xy01 = [xy - (slitWidthInPixels / (2. * np.cos(ang))),2. * yC]
                print('xy01 = ',xy01)
                xy11 = [xy + (slitWidthInPixels / (2. * np.cos(ang))),2. * yC]
                print('xy11 = ',xy11)
                #STOP
                integral,nPix = integrateBetween(xyL,xy00,xy01,xy10,xy11)
                results.append([angle,integral,nPix])
                print('first integral = 1.2713')
                print('this integral = ',integral)
                #STOP
                f.write('%d,%.5f\n' % (angle,integral))
            else:
                xy00 = [0.,yC-(slitWidthInPixels / 2.)]
                print('xy00 = ',xy00)
                xy01 = [2.*xC,yC-(slitWidthInPixels/2.)]
                print('xy01 = ',xy01)
                xy10 = [0.,yC+(slitWidthInPixels / 2.)]
                print('xy10 = ',xy10)
                xy11 = [2.*xC,yC+(slitWidthInPixels/2.)]
                print('xy11 = ',xy11)
                integral,nPix = integrateBetween(xyL,xy00,xy01,xy10,xy11)
                results.append([angle,integral,nPix])
                print('first integral = 1.2713')
                print('this integral = ',integral)
                #STOP
                f.write('%d,%.5f\n' % (angle,integral))

    plt.plot(angles,[r[1] for r in results])
    plt.xlabel('slit orientation angle [deg]')
    plt.ylabel('average [OIII]/H_beta')
    plt.savefig(resultsFileName[:resultsFileName.rfind('.')]+'_IC418_slit%d.png' % (slitWidthInArcSec))
    plt.show()

def getRatioFromIndividualImages():

    fitsNameOIII = '/Users/azuri/daten/uni/HKU/IC418/IC418_5007_sum2_without_cs_10_v2.fits'#IC418_4861_sum2_without_cs_8_v2.fits'
    OIIIMap = getImageData(fitsNameOIII,0)

    fitsNameHbeta = '/Users/azuri/daten/uni/HKU/IC418/IC418_4861_sum2_without_cs_8_v2.fits'
    HbetaMap = getImageData(fitsNameHbeta,0)

    if OIIIMap.shape != HbetaMap.shape:
        print('shapes are different!')
        STOP

    dist1Pix = angularDistancePyAsl(hmsToDeg('5:27:28.421'),dmsToDeg('-12:41:59.03'),hmsToDeg('5:27:28.3735'),dmsToDeg('-12:41:59.029')) * 3600.
    print('dist1Pix = ',dist1Pix)

    dist = angularDistanceFromXYPyAsl(fitsNameOIII,centerPos[0],centerPos[1],centerPos[0]+1.,centerPos[1])
    print('dist = ',dist)
    STOP

    slitWidthInArcSec = 4.
    slitWidthInPixels = slitWidthInArcSec * 100. / dist
    print('slitWidthInPixels = ',slitWidthInPixels)

    print('centerPos = ',centerPos)

    angles = np.arange(0,181,15)
    print('angles = ',angles)

    plt.plot(OIIIMap.shape[1],OIIIMap.shape[0],marker='v',color='white')
    plt.imshow(OIIIMap,origin='lower')
    plt.show()

    xl = np.arange(0,OIIIMap.shape[1]*100,1)
    yl = np.arange(0,OIIIMap.shape[0]*100,1)
    xyLOIII = np.zeros((xl.shape[0],yl.shape[0]))
    xyLHbeta = np.zeros((xl.shape[0],yl.shape[0]))
    print('xyLOIII.shape = ',xyLOIII.shape)
    for ix in range(xl.shape[0]):
        for iy in range(yl.shape[0]):
            xyLOIII[ix,iy] = OIIIMap[int(iy/100),int(ix/100)] / 10000.
            xyLHbeta[ix,iy] = HbetaMap[int(iy/100),int(ix/100)] / 10000.

    resultsFileName = '/Users/azuri/daten/uni/HKU/IC418/resultsFrom2Images_slit=%darcsec.csv' % (slitWidthInArcSec)
    results = []
    with open(resultsFileName,'w') as f:
        f.write('angle,averge ratio\n')
        for angle in angles:
            ang = np.radians(angle)
            print('angle = ',angle,', ang = ',ang)
            if angle != 90.:

                r = yC / (np.cos(ang) + 1e-10)
                print('r = ',r)

                xx = xC + (r * np.sin(ang))
                print('xx = ',xx)
                xy00 = [xx - (slitWidthInPixels / (2. * np.cos(ang))),0.]
                print('xy00 = ',xy00)
                xy10 = [xx + (slitWidthInPixels / (2. * np.cos(ang))),0.]
                print('xy10 = ',xy10)

                xy = xC - (r * np.sin(ang))
                print('xy = ',xy)
                xy01 = [xy - (slitWidthInPixels / (2. * np.cos(ang))),2. * yC]
                print('xy01 = ',xy01)
                xy11 = [xy + (slitWidthInPixels / (2. * np.cos(ang))),2. * yC]
                print('xy11 = ',xy11)
            else:
                xy00 = [0.,yC-(slitWidthInPixels / 2.)]
                print('xy00 = ',xy00)
                xy01 = [2.*xC,yC-(slitWidthInPixels/2.)]
                print('xy01 = ',xy01)
                xy10 = [0.,yC+(slitWidthInPixels / 2.)]
                print('xy10 = ',xy10)
                xy11 = [2.*xC,yC+(slitWidthInPixels/2.)]
                print('xy11 = ',xy11)
                #STOP
            integralOIII,integralHbeta,nPix = integrateBetween(xyLOIII,xyLHbeta,xy00,xy01,xy10,xy11)
            results.append([angle,integralOIII,integralHbeta,integralOIII/integralHbeta,nPix])
            #STOP
            f.write('%d,%.5f,%.5f,%.5f\n' % (angle,integralOIII,integralHbeta,integralOIII/integralHbeta))

    plt.plot(angles,[r[3] for r in results])
    plt.xlabel('slit orientation angle [deg]')
    plt.ylabel('[OIII]/H_beta')
    plt.savefig(resultsFileName[:resultsFileName.rfind('.')]+'_IC418_slit%d_OIII_and_Hbeta.png' % (slitWidthInArcSec))
    plt.show()


if __name__ == '__main__':
    getRatioFromIndividualImages()
    fitsNameOIII = '/Users/azuri/daten/uni/HKU/IC418/IC418_5007_sum2_without_cs_10_v2.fits'#IC418_4861_sum2_without_cs_8_v2.fits'
    OIIIMap = getImageData(fitsNameOIII,1)
    plt.imshow(OIIIMap)
    plt.show()

    fHbeta = getImageData('/Users/azuri/daten/uni/HKU/IC418/TExpCombi/MyFits/fHb.fits',1)
    plt.imshow(fHbeta)
    plt.show()
    eHbeta = getImageData('/Users/azuri/daten/uni/HKU/IC418/TExpCombi/MyFits/efHb.fits',1)
    eHbeta[np.isnan(eHbeta)] = 0.
    plt.imshow(eHbeta)
    plt.show()
    fOIII = getImageData('/Users/azuri/daten/uni/HKU/IC418/TExpCombi/MyFits/fOIII4959.fits',1)
    plt.imshow(fOIII)
    plt.show()
    eOIII = getImageData('/Users/azuri/daten/uni/HKU/IC418/TExpCombi/MyFits/efOIII4959.fits',1)
    plt.imshow(eOIII)
    plt.show()

    sumHbeta = np.sum(fHbeta)
    print('sumHbeta = ',sumHbeta)
    sumOIII = np.sum(fOIII)
    print('sumOIII = ',sumOIII)
    sumeHbeta = np.sum(eHbeta)
    print('sumeHbeta = ',sumeHbeta)
    sumeOIII = np.sum(eOIII)
    print('sumeOIII = ',sumeOIII)
    ratio = sumOIII / sumHbeta
    ratio5007 = 3. * ratio
    print('ratio = ',ratio)
    print('ratio5007 = ',ratio5007)
    eRatio = np.sqrt((sumeOIII / sumOIII)**2 + (sumeHbeta / sumHbeta)**2)
    print('eRatio = ',eRatio)

    ratio = fOIII / fHbeta
    e = ratio * np.sqrt((eOIII/fOIII)**2 + (eHbeta/fHbeta)**2)
    print('ratio = ',ratio)
    print('e = ',e)
    plt.imshow(e)
    plt.show()
