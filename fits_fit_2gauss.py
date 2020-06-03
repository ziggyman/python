import astropy.io.fits as pyfits
#from astropy import units as u
import math
import numpy as np
from scipy import exp
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

imagename = '/Users/azuri/daten/uni/HKU/HASH/StenholmAcker_pn_g024_1+03_8_id316.fits'

# calculate one Gaussian
def gauss(x,a,x0,sigma,yBackground=0.):
    return a*exp(-(x-x0)**2/(2*sigma**2))+yBackground

#calculate 2 Gaussians
def gauss2(x,a1,a2,x01,x02,sigma1,sigma2,yBackground=0.):
    return gauss(x,a1,x01,sigma1,yBackground) + gauss(x,a2,x02,sigma2,yBackground) + yBackground

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

    areas = getAreas2Gauss(x,imageData,565,2891,6548.7,6563.7,4.5,4.5)

def getAreas2Gauss(x,imageData,a1,a2,x01,x02,sigma1,sigma2,show=True):
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
    yGauss1 = gauss(x,popt[0],popt[2],popt[4])
    yGauss2 = gauss(x,popt[1],popt[3],popt[5])

    areaUnderCurve1 = np.trapz(yGauss1,x=x)
    print('areaUnderCurve1 = ',areaUnderCurve1)
    areaUnderCurve2 = np.trapz(yGauss2,x=x)
    print('areaUnderCurve2 = ',areaUnderCurve2)

    if show:
        plt.plot(x,imageData)
#    plt.plot(x,yGauss1)
#    plt.plot(x,yGauss2)

    gauss12 = gauss2(x,*popt)
    plt.plot(x,gauss12)
    if show:
        plt.show()

    return [areaUnderCurve1,areaUnderCurve2,popt]