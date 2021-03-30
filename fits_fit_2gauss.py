import astropy.io.fits as pyfits
#from astropy import units as u
import math
import numpy as np
from scipy import exp
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

imagename = '/Users/azuri/spectra/GTC/IPHASXJ055242.8+262116_GT030118.fits'

# calculate one Gaussian
def gauss(x,a,x0,sigma,yBackground=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground

#calculate 2 Gaussians
def gauss2(x,a1,a2,x01,x02,sigma1,sigma2,yBackground=0.):
    return gauss(x,a1,x01,sigma1,yBackground) + gauss(x,a2,x02,sigma2,yBackground)

#calculate 3 Gaussians
def gauss3(x,a1,a2,a3,x01,x02,x03,sigma1,sigma2,sigma3,yBackground=0.):
    return gauss(x,a1,x01,sigma1,yBackground) + gauss(x,a2,x02,sigma2,yBackground) + gauss(x,a3,x03,sigma3,yBackground)

# read image data from fits file
def getImageData(fname, hduNum=1):
    hdulist = pyfits.open(fname)
    scidata = hdulist[hduNum].data
    hdulist.close()
    return scidata

def getWavelength(header, axis=2):
    nPix = int(header['NAXIS'+str(axis)])
    crPix = int(header['CRPIX'+str(axis)])
    crVal = float(header['CRVAL'+str(axis)])
    cDelt = float(header['CDELT'+str(axis)])
#    lam = np.ndarray(nPix, dtype=np.float32)
#    lam[0] =
#    for i in np.arange(1,nPix):
    lam = np.arange(crVal, crVal + ((nPix-0.5)*cDelt), cDelt, dtype=np.float32)
#    print 'getWavelength: lam = ',len(lam),': ',lam
    return lam

def getAreas2Gauss(x,imageData,a1,a2,x01,x02,sigma1,sigma2,addOnBothSidesOfX=0.,show=True,save=None):
    # fit 2 Gaussians
    try:
        popt,pcov = curve_fit(gauss2,x,imageData,p0=[a1,
                                                     a2,
                                                     x01,
                                                     x02,
                                                     sigma1,
                                                     sigma2,
                                                     ])
    except Exception as e:
        print(e)
        STOP
#    print('popt = ',popt)

    print('amplitude a1 = ',popt[0])
    print('amplitude a2 = ',popt[1])
    print('position x01 = ',popt[2])
    print('position x02 = ',popt[3])
    print('sigma1 = ',popt[4])
    print('sigma2 = ',popt[5])

    #calculate fitted gaussians and overplot them

    xFit = x
    if addOnBothSidesOfX != 0.:
        xFit = np.arange(x[0]-addOnBothSidesOfX,x[len(x)-1]+addOnBothSidesOfX,x[1]-x[0])

    yGauss1 = gauss(xFit,popt[0],popt[2],popt[4])
    yGauss2 = gauss(xFit,popt[1],popt[3],popt[5])

    areaUnderCurve1 = np.trapz(yGauss1,x=xFit)
    print('areaUnderCurve1 = ',areaUnderCurve1)
    areaUnderCurve2 = np.trapz(yGauss2,x=xFit)
    print('areaUnderCurve2 = ',areaUnderCurve2)

    gauss12 = gauss2(xFit,*popt)

    if show or (save is not None):
        plt.plot(x,imageData)
        plt.plot(xFit,yGauss1)
        plt.plot(xFit,yGauss2)
        plt.plot(xFit,gauss12)
    if save is not None:
        plt.savefig(save,bbox_inches='tight')
    if show:
        plt.show()
    if save is not None:
        plt.close()
    print('plot finished, areaUnderCurve1 = ',areaUnderCurve1,', areaUnderCurve2 = ',areaUnderCurve2)
    return [areaUnderCurve1,areaUnderCurve2,popt]

def getAreaGauss(x,imageData,a1,x01,sigma1,addOnBothSidesOfX=0.,show=True,save=None):
    # fit 2 Gaussians
    try:
        popt,pcov = curve_fit(gauss,x,imageData,p0=[a1,
                                                    x01,
                                                    sigma1,
                                                   ])
    except Exception as e:
        print(e)
        return [0.,[0.,0.,0.]]
#    print('popt = ',popt)

    print('amplitude a1 = ',popt[0])
    print('position x01 = ',popt[1])
    print('sigma1 = ',popt[2])

    #calculate fitted gaussians and overplot them
    xFit = x
    if addOnBothSidesOfX != 0.:
        xFit = np.arange(x[0]-addOnBothSidesOfX,x[len(x)-1]+addOnBothSidesOfX,x[1]-x[0])
    yGauss1 = gauss(xFit,popt[0],popt[1],popt[2])

    areaUnderCurve1 = np.trapz(yGauss1,x=xFit)
    print('areaUnderCurve1 = ',areaUnderCurve1)
    if show or (save is not None):
        plt.plot(x,imageData)
        plt.plot(xFit,yGauss1)
    if save is not None:
        plt.savefig(save,bbox_inches='tight')
    if show:
        plt.show()
    if save is not None:
        plt.close()

    return [areaUnderCurve1,popt]

def getAreas3Gauss(x,imageData,a1,a2,a3,x01,x02,x03,sigma1,sigma2,sigma3,show=True,save=None):
    # fit 2 Gaussians
    try:
        popt,pcov = curve_fit(gauss3,x,imageData,p0=[a1,
                                                     a2,
                                                     a3,
                                                     x01,
                                                     x02,
                                                     x03,
                                                     sigma1,
                                                     sigma2,
                                                     sigma3,
                                                     ])
    except Exception as e:
        print(e)
        STOP
#    print('popt = ',popt)

    print('amplitude a1 = ',popt[0])
    print('amplitude a2 = ',popt[1])
    print('amplitude a3 = ',popt[2])
    print('position x01 = ',popt[3])
    print('position x02 = ',popt[4])
    print('position x03 = ',popt[5])
    print('sigma1 = ',popt[6])
    print('sigma2 = ',popt[7])
    print('sigma3 = ',popt[8])

    #calculate fitted gaussians and overplot them
    yGauss1 = gauss(x,popt[0],popt[3],popt[6])
    yGauss2 = gauss(x,popt[1],popt[4],popt[7])
    yGauss3 = gauss(x,popt[2],popt[5],popt[8])

    areaUnderCurve1 = np.trapz(yGauss1,x=x)
    print('areaUnderCurve1 = ',areaUnderCurve1)
    areaUnderCurve2 = np.trapz(yGauss2,x=x)
    print('areaUnderCurve2 = ',areaUnderCurve2)
    areaUnderCurve3 = np.trapz(yGauss3,x=x)
    print('areaUnderCurve3 = ',areaUnderCurve3)

    gauss123 = gauss3(x,*popt)

    if show or (save is not None):
        plt.plot(x,imageData)
        plt.plot(x,yGauss1)
        plt.plot(x,yGauss2)
        plt.plot(x,yGauss3)
        plt.plot(x,gauss123)
    if save is not None:
        plt.savefig(save,bbox_inches='tight')
    if show:
        plt.show()
    if save is not None:
        plt.close()

    return [areaUnderCurve1,areaUnderCurve2,areaUnderCurve3,popt]

# main function
if __name__ == "__main__":
    #read image data from fits file
    imageData = getImageData(imagename, 0)
#    print('imageData = ',imageData)
#    print('type(imageData) = ',type(imageData))
#    print('imageData.shape = ',imageData.shape)

    #create x array
    header = pyfits.getheader(imagename,0)
#    print('header = ',header)
    x = getWavelength(header,1)#np.arange(0,imageData.shape[0],1)
#    print('x = ',x.shape,': ',x)

    #plot the spectrum
    plt.plot(x,imageData)
    plt.show()

    areas = getAreas3Gauss(x,imageData,6.83E-17,5.74E-15,2.13E-16,6548.7,6563.7,6583.,3.,3.,3.)
