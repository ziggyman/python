import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
import os

from scipy.optimize import curve_fit

from drUtils import getImageData,getWavelengthArr,getHeader,getHeaderValue,getInsideRanges,gauss,gauss_lin

def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

if __name__ == '__main__':
    path = '/Users/azuri/daten/uni/HKU/Kamila'
    fName = os.path.join(path,'Knot_Spectra','spectrumB_J1.fits')
