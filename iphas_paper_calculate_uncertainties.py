import matplotlib.pyplot as plt
import numpy as np
import os
import pyneb as pn

import csvFree, csvData

from drUtils import sfit#(x, y, fittingFunction, solveFunction, order, nIterReject, nIterFit, lowReject, highReject, adjustSigLevels=False, useMean=False, display=False)
from drUtils import continuum#(spectrumFileNameIn, spectrumFileNameOut, fittingFunction, evalFunction, order, nIterReject, nIterFit, lowReject, highReject, type='difference', adjustSigLevels=False, useMean=False, display = False)
from drUtils import writeFits1D#(spec,spectrumFileNameOut,wavelength=None,header=spectrumFileNameIn,CRVAL1=getHeaderValue(spectrumFileNameIn,'CRVAL1'),CRPIX1=getHeaderValue(spectrumFileNameIn,'CRPIX1'),CDELT1=getHeaderValue(spectrumFileNameIn,'CDELT1'),)
from drUtils import getHeaderValue#(fname, keyword, hduNum=0)
from drUtils import getWavelengthArr#(fname, hduNum=0)
from drUtils import getImageData#(fname,hduNum=1)
from drUtils import getHeader#(fName, hduNum=0)
#from myUtils import 
#from myUtils import 
#from myUtils import 
#from myUtils import 

linesOfInterest = {'H1r_6563A':6562.801,
                   'H1r_4861A':4861.363,
                   'S2_6716A':6716.44,
                   'S2_6731A':6730.82,
                   'N2_5755A':5754.59,
                   'N2_6548A':6548.05,
                   'N2_6584A':6583.45,
                   'O3_4363A':4363.209,
                   'O3_5007A':5006.843	,
                    }

def getSigmaVsWLen(spectrumName,nPix=20,show=False):
    wLen = getWavelengthArr(spectrumName,0)
    print(spectrumName,': wLen = ',wLen)
    arrayLength = wLen.shape[0]
    spec = getImageData(spectrumName,0)
    if show:
        plt.plot(wLen,spec)
        plt.show()
    sigmas = np.zeros(arrayLength)
    for i in range(arrayLength):
        if i < int(nPix/2.):
            sigmas[i] = np.std(spec[:nPix])
        elif i > arrayLength - int(nPix/2.):
            sigmas[i] = np.std(spec[:-nPix])
        else:
            sigmas[i] = np.std(spec[i-int(nPix/2.):i+int(nPix/2.)+1])
    sigmaSpecFileName = os.path.join(spectrumName[:spectrumName.rfind('/')],'test','sigma'+spectrumName[spectrumName.rfind('/')+1:])
    sigmaFitSpecFileName = os.path.join(spectrumName[:spectrumName.rfind('/')],'test','sigmaFit'+spectrumName[spectrumName.rfind('/')+1:])
    writeFits1D(sigmas,
                sigmaSpecFileName,
                wavelength=None,
                header=spectrumName,
                CRVAL1=getHeaderValue(spectrumName,'CRVAL1'),
                CRPIX1=getHeaderValue(spectrumName,'CRPIX1'),
                CDELT1=getHeaderValue(spectrumName,'CDELT1'),
               )


    fittingFunction = np.polynomial.legendre.legfit
    evalFunction = np.polynomial.legendre.legval
    order = 9
    nIterReject = 2
    nIterFit = 2
    lowReject = 3.
    highReject = 1.
    useMean = False
    continuum(sigmaSpecFileName,
              sigmaFitSpecFileName,
              fittingFunction,
              evalFunction,
              order,
              nIterReject,
              nIterFit,
              lowReject,
              highReject,
              type='fit',
              adjustSigLevels=False,
              useMean=useMean,
              display=show)

    if show:
        fit = getImageData(sigmaFitSpecFileName,0)
        plt.plot(wLen,spec)
        plt.plot(wLen,sigmas)
        plt.plot(wLen,fit)
        plt.show()
    return [wLen,fit]

def calculateErrors(spectrumFileName,idPNMain,csvLinesFileName):
    wLen, sigmaFit = getSigmaVsWLen(spectrumFileName,show=True)
    csvLines = csvFree.readCSVFile(csvLinesFileName,'\t',False)
    header = csvLines.header
    for i in range(csvLines.size()):
        if idPNMain == csvLines.getData('NAME'):
            for iLine in range(len(linesOfInterest)):
                keys = list(linesOfInterest.keys())
                area = csvLines.getData(keys[iLine],i)
                print('key = ',key,': area = ',area)


if __name__ == '__main__':
    spectrumFileName = '/Users/azuri/spectra/GTC/LDu1_sum.fits'
    csvLinesFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/observation.dat'
    hash_fitsFiles = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/HASH/fitsfiles.csv')
    idPNMain = None
    for i in range(hash_fitsfiles.size()):
        if spectrumFileName[spectrumFileName.rfind('/')+1:] == hash_fitsFiles.getData('fileName',i):
            idPNMain = hash_fitsFiles.getData('idPNMain',i)
    if idPNMain is None:
        print('ERROR: did not find ',spectrumFileName)
        STOP
    errors = calculateErrors(spectrumFileName,idPNMain,csvLinesFileName)
