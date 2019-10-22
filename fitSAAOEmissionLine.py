import numpy as np
import matplotlib.pyplot as plt

#from hammer import Hammer
import MpFit# import MPFitTwoGaussLim

specFileNames = ['/Volumes/work/azuri/spectra/saao/saao_sep2019/20190905/redux/a1171070_1D.txt',
                '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190905/redux/a1171071_1D.txt',
                '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190905/redux/a1171072_1D.txt',
                '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190905/redux/a1171073_1D.txt',
                '/Volumes/work/azuri/spectra/saao/saao_sep2019/20190905/redux/a1171074_1D.txt',
]

def findValue(arr, val):
    for i in np.arange(1,arr.shape[0],1):
        if (arr[i-1] < val) and (arr[i] >= val):
            return i
    return -1

for specFileName in specFileNames:
    with open(specFileName,'r') as f:
        lines = f.readlines()

    wavelength = np.zeros(len(lines)-1)
    flux = np.zeros(len(lines)-1)

    for iLine in np.arange(1,len(lines),1):
    #    print('lines[',iLine,'] = ',lines[iLine])
        lines[iLine] = lines[iLine].split(' ')
        wavelength[iLine-1] = lines[iLine][0]
        flux[iLine-1] = lines[iLine][1]

    plt.plot(wavelength, flux)
    plt.show()

    xStart = findValue(wavelength,7875.)
    xEnd = findValue(wavelength,8020.)
    x = wavelength[xStart:xEnd]
    y = flux[xStart:xEnd]
    eY = np.sqrt(y)

    guess = np.ndarray(shape=6, dtype=np.float)
    guess[0] = 850. #     p[0] = constant offset
    guess[1] = 1600. #     p[1] = peak y value 1st Gauss
    guess[2] = 7929. #     p[2] = x centroid position 1st Gauss
    guess[3] = 40. #     p[3] = gaussian sigma width
    guess[4] = 540. #     p[4] = peak y value 2nd Gauss
    guess[5] = 7963. #     p[5] = x centroid position 2nd Gauss

    limited = np.ndarray(shape=(6,2), dtype=np.int)
    limited[0,0] = 1
    limited[0,1] = 1
    limited[1,0] = 1
    limited[1,1] = 1
    limited[2,0] = 1
    limited[2,1] = 1
    limited[3,0] = 1
    limited[3,1] = 1
    limited[4,0] = 1
    limited[4,1] = 1
    limited[5,0] = 1
    limited[5,1] = 1

    limits = np.ndarray(shape=(6,2), dtype=np.float)
    limits[0,0] = 800.
    limits[0,1] = 900.
    limits[1,0] = 1400.
    limits[1,1] = 1800.
    limits[2,0] = 7910.
    limits[2,1] = 7940.
    limits[3,0] = 20.
    limits[3,1] = 60.
    limits[4,0] = 400.
    limits[4,1] = 650.
    limits[5,0] = 7940.
    limits[5,1] = 7980.

    coeffs = np.zeros(shape=(6,2), dtype=float)
    eCoeffs = np.zeros(shape=(6,2), dtype=float)

    if False:
        result = MpFit.MPFitTwoGaussLim(x,
                                  y,
                                  eY,
                                  guess,
                                  limited,
                                  limits,
                                  True,
                                  False,
                                  coeffs,
                                  eCoeffs)
        print('result = ',result)

        print('coeffs = ',coeffs)
        print('eCoeffs = ',eCoeffs)
