import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from drUtils import getImageData,getWavelengthArr,getHeader
from myUtils import angularDistanceFromXYPyAsl,hmsToDeg,dmsToDeg,angularDistancePyAsl
import astropy.io.fits as pyfits
from matplotlib.patches import Polygon
#from pywifes_utils import extract_stellar_spectra_fits

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


def integrateBetween(xyLOIII,xyLHbeta,xy00,xy01,xy10,xy11,fName,slitAngle):
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
    pc = PatchCollection(squares,facecolor='r',alpha=0.4,edgecolor='b')
    ax.add_collection(pc)
    plt.savefig(fName[:fName.rfind('.')]+'_%ddeg.png' % (slitAngle))
    plt.close()
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

    dist = 0.5#angularDistanceFromXYPyAsl(fitsNameOIII,centerPos[0],centerPos[1],centerPos[0]+1.,centerPos[1])
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
    if False:
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

    # WiFeS
    if False:
        fitsFileNames = ['/Users/azuri/daten/uni/HKU/IC418/WiFes/T2m3wb-20190905.190910-0020.p11.fits',
                    '/Users/azuri/daten/uni/HKU/IC418/WiFes/T2m3wb-20190907.185630-0044.p11.fits',
                    #'/Users/azuri/daten/uni/HKU/IC418/WiFes/T2m3wb-20190907.190323-0045.p11.fits',#saturated
                    ]
        for fitsFileName in fitsFileNames:

    #        p08 = fitsFileName[:fitsFileName.rfind('.p')]+'.p08.fits'
    #        hdul = pyfits.open(p08)
    #        print('len(hdul) = ',len(hdul))
    #        print('p08 = '+p08)
    #        sigFrac = np.zeros([2840,76,25])
    #        for i in np.arange(1,26,1):
    #            sig = np.nan_to_num(np.array(getImageData(p08,i)))
    #            dat = np.nan_to_num(np.array(getImageData(p08,i+25)))
    #            print('i=',i,': dat = ',dat.shape,': ',dat)
    #            sigFrac[:,:,25-i] = np.transpose(np.nan_to_num(sig/dat))
    #        print('sigFrac = ',sigFrac)
    #        STOP
    #        sig = np.nan_to_num(np.array(getImageData(p08,1),dtype=float))
    #        print('sig = ',type(sig),': ',sig.shape,': ',sig)
    #        data = np.nan_to_num(np.array(getImageData(p08,0),dtype=float))
    #        print('data = ',type(data),': ',data.shape,': ',data)
    #        sigFrac = np.nan_to_num(sig / data)
            hdul = pyfits.open(fitsFileName)
            print('len(hdul) = ',len(hdul))
            scidata = np.nan_to_num(getImageData(fitsFileName,0))
    #        sigdata = scidata * sigFrac
    #        sigHeader = getHeader(fitsFileName,1)
    #        print('sigHeader = ',sigHeader)
    #        obsdata = getImageData(fitsFileName,2)
            header = getHeader(fitsFileName,0)
            print('header = ',header)
            wLen = getWavelengthArr(fitsFileName,0,dim=3)
            print('wLen = ',wLen.shape,': ',wLen)
            print('scidata = ',type(scidata),': ',scidata.shape,': ',scidata)
    #        print('sigdata = ',type(sigdata),': ',sigdata.shape,': ',sigdata)
    #        print('obsdata = ',type(obsdata),': ',obsdata.shape,': ',obsdata)
    #        hdul[0].data = sigFrac
    #        hdul.writeto(fitsFileName[:fitsFileName.rfind('.fits')]+'_sigFrac.fits',overwrite=True)
            #STOP

            if fitsFileName == '/Users/azuri/daten/uni/HKU/IC418/WiFes/T2m3wb-20190905.190910-0020.p11.fits':
                objArea = [[3,8],[18,42]]
                ignorePixVal = 1e-13
                if False:
                    ignorePix = [[0,0],[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,12],[0,13],[0,14],[0,15],[0,16],
                                [1,0],[1,1],[1,2],[1,3],[1,4],[1,13],[1,14],[1,15],[1,16],
                                [2,0],[2,1],[2,2],[2,3],[2,14],[2,15],[2,16],
                                [3,0],[3,1],[3,2],[3,3],[3,14],[3,15],[3,16],
                                [4,0],[4,1],[4,2],[4,3],[4,14],[4,15],[4,16],
                                [5,],
                                [6,],
                                [7,],
                                [8,],
                                [9,],
                                [10,],
                                [11,],
                                [12,],
                                [13,],
                                [14,],
                                [15,],
                                [16,],
                                [17,],
                                [18,],
                                [19,],
                                [20,],
                                [21,],
                                [22,],
                                [23,],
                                [24,],
                                [25,],
                                [26,],
                                [27,],
                                [28,],
                                [29,],
                                [30,],
                                [31,],
                                [32,],
                                [33,],
                                [34,],
                                ]
            else:
                objArea = [[7,2],[22,34]]
                ignorePixVal = 2e-14
                if False:
                    ignorePix = [[0,0],[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,12],[0,13],[0,14],[0,15],[0,16],
                                [1,0],[1,1],[1,2],[1,3],[1,4],[1,13],[1,14],[1,15],[1,16],
                                [2,0],[2,1],[2,2],[2,3],[2,14],[2,15],[2,16],
                                [3,0],[3,1],[3,2],[3,3],[3,14],[3,15],[3,16],
                                [4,0],[4,1],[4,2],[4,3],[4,14],[4,15],[4,16],
                                [5,0],[5,1],[5,2],[5,15],[5,16],
                                [6,0],[6,1],[6,15],[6,16],
                                [7,0],[7,1],[7,15],[7,16],
                                [8,0],[8,16],
                                [9,0],[9,16],
                                [10,0],[10,16],
                                [11,0],[11,16],
                                [12,0],[12,16],
                                [13,],
                                [14,],
                                [15,],
                                [16,],
                                [17,],
                                [18,],
                                [19,],
                                [20,],
                                [21,],
                                [22,],
                                [23,],
                                [24,],
                                [25,],
                                [26,],
                                [27,],
                                [28,],
                                [29,],
                                [30,],
                                [31,],
                                [32,],
                                [33,],
                                [34,],
                                ]

    #        sig2D = np.nan_to_num(sigdata[:,objArea[0][1]:objArea[1][1],objArea[0][0]:objArea[1][0]])
    #        sig1D = np.sum(sig2D,axis=2)
    #        print('sig1D = ',sig1D.shape,': ',sig1D)
    #        sig = np.sum(sig1D,1)
    #        print('sig = ',sig.shape,': ',sig)

    #        for i in flux:
    #            print(i)

    #        nans = np.where(flux is np.nan)
    #        print('nans = ',nans)

            background = np.median(scidata[np.where((wLen > 4810.) & (wLen < 5077))[0],:,:])
            print('background = ',background)

            flux2D = scidata[:,objArea[0][1]:objArea[1][1],objArea[0][0]:objArea[1][0]] - background
            print('flux2D = ',flux2D.shape,': ',flux2D)
            print('mean(flux2D) = ',np.mean(flux2D))
            idx = np.where(flux2D < 0.)
            print('idx = ',idx)
            flux2D[idx] = 0.
            print('mean(flux2D) = ',np.mean(flux2D))
            hdul[0].data = flux2D
            hdul.writeto(fitsFileName[:fitsFileName.rfind('.fits')]+'_flux2D.fits',overwrite=True)

            flux1D = np.sum(flux2D,axis=2)
            print('flux1D = ',flux1D.shape,': ',flux1D)
            flux = np.sum(flux1D,1)
            print('flux = ',flux.shape,': ',flux)

            fig, ax = plt.subplots()
            background = np.median(flux[np.where((wLen > 4810.) & (wLen < 5077))[0]])
            flux = flux - background
            ax.plot(wLen,flux)

            idxGT = np.where(wLen >= 4858)[0]
            idxLT = np.where(wLen[idxGT] <= 4867)[0]
            Hbeta2D = np.sum(flux2D[idxGT[idxLT],:,:],axis=0)
            Hbeta = np.trapz(flux[idxGT[idxLT]],wLen[idxGT[idxLT]])
            sumHbeta = np.sum(flux[idxGT[idxLT]])
    #        sigHbeta = np.trapz(sig[idxGT[idxLT]],wLen[idxGT[idxLT]])
    #        sumSigHbeta = np.sum(sig[idxGT[idxLT]])
            print('Hbeta = ',Hbeta,', sumHbeta = ',sumHbeta)
    #        print('sigHbeta = ',sigHbeta,', sumSigHbeta = ',sumSigHbeta)

            for i in np.arange(idxGT[idxLT[0]],idxGT[idxLT[len(idxLT)-1]],1):
                coordsPoly = [(wLen[i],0.),(wLen[i],flux[i]),(wLen[i+1],flux[i+1]),(wLen[i+1],0.)]
                p=Polygon(coordsPoly,facecolor='g')
                ax.add_patch(p)


            idxGT = np.where(wLen >= 5002)[0]
            idxLT = np.where(wLen[idxGT] <= 5012)[0]
            OIII2D = np.sum(flux2D[idxGT[idxLT],:,:],axis=0)
            goodPix = np.where(Hbeta2D > ignorePixVal)
            ignorePix = np.where(Hbeta2D < ignorePixVal)
            OIII_over_Hbeta2D = OIII2D / Hbeta2D

            print('OIII_over_Hbeta2D = ',OIII_over_Hbeta2D.shape,': ',OIII_over_Hbeta2D)
            mean_OIII_over_Hbeta2D = np.mean(OIII_over_Hbeta2D[goodPix])
            sigma_OIII_over_Hbeta2D = np.std(OIII_over_Hbeta2D[goodPix])
            print('mean_OIII_over_Hbeta2D = ',mean_OIII_over_Hbeta2D)
            print('sigma_OIII_over_Hbeta2D = ',sigma_OIII_over_Hbeta2D)
            OIII5007 = np.trapz(flux[idxGT[idxLT]],wLen[idxGT[idxLT]])
            sumOIII = np.sum(flux[idxGT[idxLT]])
    #        sigOIII5007 = np.trapz(sig[idxGT[idxLT]],wLen[idxGT[idxLT]])
    #        sumSigOIII = np.sum(sig[idxGT[idxLT]])
    #        print('OIII5007 = ',OIII5007,', sumOIII = ',sumOIII)
    #        print('sigOIII5007 = ',sigOIII5007,', sumSigOIII = ',sumSigOIII)

            ratio = OIII5007 / Hbeta
            print('OIII5007 / Hbeta = ',ratio)
    #        error = ratio * ((sigOIII5007/OIII5007) + (sigHbeta/Hbeta))
    #        print('error OIII / Hbeta = ',error)
    #        sigma = ratio * np.sqrt((sigOIII5007/OIII5007)**2 + (sigHbeta/Hbeta)**2)
    #        print('sigma OIII / Hbeta = ',sigma)

            ratio = sumOIII / sumHbeta
            print('sumOIII / sumHbeta = ',ratio)
    #        error = ratio * ((sigOIII5007/OIII5007) + (sigHbeta/Hbeta))
    #        print('error OIII / Hbeta = ',error)
    #        sigma = ratio * np.sqrt((sigOIII5007/OIII5007)**2 + (sigHbeta/Hbeta)**2)
    #        print('sigma OIII / Hbeta = ',sigma)

    #        ratio = Hbeta / OIII5007
    #        print('Hbeta / OIII5007 = ',ratio)

            for i in np.arange(idxGT[idxLT[0]],idxGT[idxLT[len(idxLT)-1]],1):
                coordsPoly = [(wLen[i],0.),(wLen[i],flux[i]),(wLen[i+1],flux[i+1]),(wLen[i+1],0.)]
                p=Polygon(coordsPoly,facecolor='g')
                ax.add_patch(p)

            plt.show()

            fig,axs = plt.subplots(1,3)
            fig.set_figheight(9)
            fig.set_figwidth(15)
            OIII2D[ignorePix] = 0.
            Hbeta2D[ignorePix] = 0.
            OIII_over_Hbeta2D[ignorePix] = 0.
            oiii = axs[0].imshow(OIII2D)
            fig.colorbar(oiii,ax=axs[0])
            axs[0].title.set_text('[OIII] 5007')

            hbet = axs[1].imshow(Hbeta2D)
            fig.colorbar(hbet,ax=axs[1])
            axs[1].title.set_text('H_beta')

            rat = axs[2].imshow(OIII_over_Hbeta2D)
            fig.colorbar(rat,ax=axs[2])
            axs[0].title.set_text('[OIII] 5007 / H_beta')

            plt.show()

    #SAAO Sep 2024
    if False:
        fitsFileNames = ['/Users/azuri/daten/uni/HKU/IC418/SAAO/HW1/20240922/IC418_combined.fits',
                        '/Users/azuri/daten/uni/HKU/IC418/SAAO/HW1/20240922/SCIENCE_IC418_a5851090_otzxfif-skyEcdF.fits',
                        '/Users/azuri/daten/uni/HKU/IC418/SAAO/HW1/20240922/SCIENCE_IC418_a5851091_otzxfif-skyEcdF.fits',
                        '/Users/azuri/daten/uni/HKU/IC418/SAAO/HW1/20240922/SCIENCE_IC418_a5851092_otzxfif-skyEcdF.fits',
        ]
        ratios = []
        for fitsFileName in fitsFileNames:
            wLenRange = [4700,5200]
            wLen = getWavelengthArr(fitsFileName,0)
            spec = getImageData(fitsFileName,0)

            idxL = np.where(wLen > wLenRange[0])[0]
            idxG = np.where(wLen[idxL] < wLenRange[1])[0]

            wLen = wLen[idxL[idxG]]
            spec = spec[idxL[idxG]]
            background = np.median(spec[np.where(spec < 2e-13)[0]])
            spec -= background
            fig, ax = plt.subplots()
            ax.plot(wLen,spec)

            wLenRangeHbeta = [4831.4,4902.3]
            idxGT = np.where(wLen > wLenRangeHbeta[0])[0]
            idxLT = np.where(wLen[idxGT] < wLenRangeHbeta[1])[0]
            Hbeta = np.trapz(spec[idxGT[idxLT]],wLen[idxGT[idxLT]])
            for i in np.arange(idxGT[idxLT[0]],idxGT[idxLT[len(idxLT)-1]],1):
                coordsPoly = [(wLen[i],0.),(wLen[i],spec[i]),(wLen[i+1],spec[i+1]),(wLen[i+1],0.)]
                p=Polygon(coordsPoly,facecolor='g')
                ax.add_patch(p)


            wLenRangeOIII = [4984,5052.1]
            idxGT = np.where(wLen > wLenRangeOIII[0])[0]
            idxLT = np.where(wLen[idxGT] < wLenRangeOIII[1])[0]
            OIII = np.trapz(spec[idxGT[idxLT]],wLen[idxGT[idxLT]])
            for i in np.arange(idxGT[idxLT[0]],idxGT[idxLT[len(idxLT)-1]],1):
                coordsPoly = [(wLen[i],0.),(wLen[i],spec[i]),(wLen[i+1],spec[i+1]),(wLen[i+1],0.)]
                p=Polygon(coordsPoly,facecolor='g')
                ax.add_patch(p)

            plt.show()
            ratio = OIII / Hbeta
            print('ratio = ',ratio)
            ratios.append(ratio)
        ratios = np.array(ratios)
        mean = np.mean(ratios[1:])
        stddev = np.std(ratios[1:])
        print('mean = ',mean,', stddev = ',stddev)


    #WiFeS
    fitsFileNames = ['/Users/azuri/daten/uni/HKU/IC418/WiFes/T2m3wb-20190905.190910-0020.p11_flux2D.fits',
                '/Users/azuri/daten/uni/HKU/IC418/WiFes/T2m3wb-20190907.185630-0044.p11_flux2D.fits',
                #'/Users/azuri/daten/uni/HKU/IC418/WiFes/T2m3wb-20190907.190323-0045.p11.fits',#saturated
                ]
    #slit width in arcsec:
    slitWidth=2.69
    for fitsFileName in fitsFileNames:
        if fitsFileName == '/Users/azuri/daten/uni/HKU/IC418/WiFes/T2m3wb-20190905.190910-0020.p11_flux2D.fits':
            center = [7.6,14.7]
            ignorePixVal = 1e-14
        else:
            center = [7.5,17.5]
            ignorePixVal = 2e-14
        print('center = ',center)

        hdul = pyfits.open(fitsFileName)
        print('len(hdul) = ',len(hdul))
        flux2D = np.nan_to_num(getImageData(fitsFileName,0))
        flux2D[np.where(flux2D < ignorePixVal)] = 0.
#        sigdata = scidata * sigFrac
#        sigHeader = getHeader(fitsFileName,1)
#        print('sigHeader = ',sigHeader)
#        obsdata = getImageData(fitsFileName,2)
        plt.imshow(np.sum(flux2D,axis=0))
        plt.show()
        header = getHeader(fitsFileName,0)
        print('header = ',header)
        wLen = getWavelengthArr(fitsFileName,0,dim=3)
        print('wLen = ',wLen.shape,': ',wLen)
        print('scidata = ',type(flux2D),': ',flux2D.shape,': ',flux2D)


        idxGT = np.where(wLen >= 4858)[0]
        idxLT = np.where(wLen[idxGT] <= 4867)[0]
        Hbeta2D = np.sum(flux2D[idxGT[idxLT],:,:],axis=0)
        Hbeta1D = np.sum(Hbeta2D,axis=1)
        plt.plot(wLen[idxGT[idxLT]],np.sum(np.sum(flux2D[idxGT[idxLT],:,:],axis=2),axis=1),'b-',label='Hbeta')

        idxGT = np.where(wLen >= 5002)[0]
        idxLT = np.where(wLen[idxGT] <= 5012)[0]
        OIII2D = np.sum(flux2D[idxGT[idxLT],:,:],axis=0)
        OIII1D = np.sum(OIII2D,axis=1)
        print('OIII1D = ',OIII1D.shape,': ',OIII1D)
        plt.plot(wLen[idxGT[idxLT]],np.sum(np.sum(flux2D[idxGT[idxLT],:,:],axis=2),axis=1),'r-',label='OIII')
        plt.show()
        STOP

        lowLim = center[1]-slitWidth
        lowWeight = 1. - (lowLim - int(lowLim))
        print('lowLim = ',lowLim,', lowWeight = ',lowWeight)

        highLim = center[1]+slitWidth
        highWeight = highLim - int(highLim)
        print('highLim = ',highLim,', highWeight = ',highWeight)

        sumHbeta = lowWeight * Hbeta1D[int(lowLim)]
        sumOIII = lowWeight * OIII1D[int(lowLim)]
        print('sumHbeta = ',sumHbeta,', sumOIII = ',sumOIII)

        sumHbeta += highWeight * Hbeta1D[int(highLim)]
        sumOIII += highWeight * OIII1D[int(highLim)]
        print('sumHbeta = ',sumHbeta,', sumOIII = ',sumOIII)

        for i in range(int(lowLim)+1,int(highLim),1):
            sumHbeta += Hbeta1D[i]
            sumOIII += OIII1D[i]
            print('i = ',i,': sumHbeta = ',sumHbeta,', sumOIII = ',sumOIII)

        ratio = sumOIII/sumHbeta
        print(fitsFileName,': ratio = ',ratio)

        plt.plot(Hbeta1D)
        plt.show()
        plt.plot(OIII1D)
        plt.show()







        OIIIMap = np.zeros((OIII2D.shape[0]*100,200*OIII2D.shape[1]))
        HbetaMap = np.zeros(OIIIMap.shape)
        for ix in range(OIIIMap.shape[0]):
            for iy in range(OIIIMap.shape[1]):
                OIIIMap[ix,iy] = OIII2D[int(ix/100),int(iy/200)] / 20000.
                HbetaMap[ix,iy] = Hbeta2D[int(ix/100),int(iy/200)] / 20000.
        #plt.imshow(OIIIMap)
        #plt.show()
        #plt.imshow(HbetaMap)
        #plt.show()

        if OIIIMap.shape != HbetaMap.shape:
            print('shapes are different!')
            STOP

        dist = 0.5#angularDistanceFromXYPyAsl(fitsNameOIII,centerPos[0],centerPos[1],centerPos[0]+1.,centerPos[1])
        print('dist = ',dist)

        slitWidthInArcSec = 2.#69
        slitWidthInPixels = slitWidthInArcSec * 100. / dist
        print('slitWidthInPixels = ',slitWidthInPixels)

        centerPos = center*100
        print('centerPos = ',centerPos)

        angles = np.arange(0,181,35)
        print('angles = ',angles)

        #plt.plot(OIIIMap.shape[1],OIIIMap.shape[0],marker='v',color='white')
        #plt.imshow(OIIIMap,origin='lower')
        #plt.show()

        xyLOIII = OIIIMap
        xyLHbeta = HbetaMap
        fig,axs = plt.subplots(1,3)
        plt1 = axs[0].imshow(xyLOIII)
        axs[0].set_title('[OIII]')
        plt2 = axs[1].imshow(xyLHbeta)
        axs[1].set_title('Hbeta')
        plt3 = axs[2].imshow(xyLOIII/xyLHbeta,vmin=1.,vmax = 4.)
        axs[2].set_title('[OIII] / Hbeta')
        plt.colorbar(plt1,ax=axs[0])
        plt.colorbar(plt2,ax=axs[1])
        plt.colorbar(plt3,ax=axs[2])
        plt.show()
        print('OIIIMap.shape = ',OIIIMap.shape)
        print('xyLOIII.shape = ',xyLOIII.shape)
#        for ix in range(xl.shape[0]):
#            for iy in range(yl.shape[0]):
#                xyLOIII[ix,iy] = OIIIMap[int(iy/100),int(ix/100)] / 10000.
#                xyLHbeta[ix,iy] = HbetaMap[int(iy/100),int(ix/100)] / 10000.

        resultsFileName = fitsFileName[:fitsFileName.rfind('.')]+'_slit=%.2farcsec.csv' % (slitWidthInArcSec)
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
                integralOIII,integralHbeta,nPix = integrateBetween(xyLOIII,xyLHbeta,xy00,xy01,xy10,xy11, fitsFileName,angle)
                results.append([angle,integralOIII,integralHbeta,integralOIII/integralHbeta,nPix])
                #STOP
                f.write('%d,%.5f,%.5f,%.5f\n' % (angle,integralOIII,integralHbeta,integralOIII/integralHbeta))

        plt.plot(angles,[r[3] for r in results])
        plt.xlabel('slit orientation angle [deg]')
        plt.ylabel('[OIII]/H_beta')
        plt.savefig(resultsFileName[:resultsFileName.rfind('.')]+'_IC418_slit%d_OIII_and_Hbeta.png' % (slitWidthInArcSec))
        plt.show()
