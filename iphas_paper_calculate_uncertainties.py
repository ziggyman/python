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
from drUtils import applyVRadCorrection#(wavelength, vRad)
from fits_fit_2gauss import gauss#(x,a,x0,sigma,yBackground=0.)
from fits_fit_2gauss import getAreaGauss#(x,imageData,a1,x01,sigma1,addOnBothSidesOfX=0.,show=True,save=None)
#from myUtils import 
#from myUtils import 
#from myUtils import 
#from myUtils import 

imPath = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/iphas-gtc-images'


linesOfInterest = {'H1r_6563A':6562.801,
                   'H1r_4861A':4861.363,
                   'S2_6716A':6716.44,
                   'S2_6731A':6730.82,
                   'N2_5755A':5754.59,
                   'N2_6548A':6548.05,
                   'N2_6584A':6583.45,
                   'O3_4363A':4363.209,
                   'O3_5007A':5006.843,
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
    order = 4
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

    fit = getImageData(sigmaFitSpecFileName,0)
    if show:
        plt.plot(wLen,spec)
        plt.plot(wLen,sigmas)
        plt.plot(wLen,fit)
        plt.show()
    return [wLen,fit]

def calculateErrors(spectrumFileName,idPNMain,csvLinesFileName,show=False):
    print('spectrumFileName = <'+spectrumFileName+'>')
    wLen, sigmaFit = getSigmaVsWLen(spectrumFileName,show=False)
    print('len(wLen) = ',len(wLen),', len(sigmaFit) = ',len(sigmaFit))
    flux = getImageData(spectrumFileName,0)
    print('len(flux) = ',len(flux))
    csvLines = csvFree.readCSVFile(csvLinesFileName,'\t',False)
    csvVRad = csvFree.readCSVFile(os.path.join(imPath[:imPath.rfind('/')],'vrad.csv'))
    print('csvVRad.header = ',csvVRad.header)
    filenames = csvVRad.getData('fileName')
    print('filenames = ',filenames)
    vradPos = csvVRad.find('fileName',spectrumFileName[spectrumFileName.rfind('/')+1:],0)[0]
    if vradPos < 0:
        print('error: did not find spectrumFileName <'+spectrumFileName+'>')
        #STOP
    else:
        vrad = float(csvVRad.getData('vrad',vradPos))
        print('vrad = ',type(vrad),': ',vrad)
        wLen = applyVRadCorrection(wLen, vrad)
        print('vradPos = ',vradPos)
        header = csvLines.header
        keys = list(linesOfInterest.keys())
        for i in range(csvLines.size()):
            if idPNMain == csvLines.getData('NAME',i):
                for iLine in range(len(linesOfInterest)):
                    area = float(csvLines.getData(keys[iLine],i))
                    print('key = ',keys[iLine],': area = ',area)
                    if area > 0.:
                        x0 = linesOfInterest[keys[iLine]]
                        x = wLen[np.where(np.abs(wLen - x0) < 20.)[0]]
                        thisFlux = flux[np.where(np.abs(wLen - x0) < 3.)[0]]
                        maxFlux = np.max(thisFlux)
                        sigma = area / (maxFlux * 2.13 * np.sqrt(2. * np.log(2.)))
                        print('x = ',x0,', a = ',maxFlux,', sigma = ',sigma)
                        thisFlux = flux[np.where(np.abs(wLen - x0) < 20.)[0]]
                        thisSDev = sigmaFit[np.where(np.abs(wLen - x0) < 20.)[0]]
                        gaussFit = gauss(x,maxFlux,x0,sigma)
                        if show:
                            plt.plot(x,thisFlux,label='flux')
                            plt.plot(x,thisSDev,label='sigma')
                            plt.plot(x,gaussFit,label='fit')
                            plt.legend()
                            plt.show()
                        newArea = getAreaGauss(x,thisFlux,maxFlux,x0,sigma,addOnBothSidesOfX=0.,show=False,save=None)
                        print('old area = ',area,', newly fitted area = ',newArea)
                        if show:
                            plt.plot(x,gaussFit,label='fit')
                            plt.plot(x,thisFlux,label='flux')
                        newAreas = []
                        for iRun in range(100):
                            thisFluxWithErr = np.zeros(x.shape,dtype='float32')
                            for thisFluxPos in range(x.shape[0]):
                                thisFluxWithErr[thisFluxPos] = gaussFit[thisFluxPos] + np.random.normal(0.,np.abs(thisSDev[thisFluxPos]))
                            if show:
                                plt.plot(x,thisFluxWithErr,label='%d' % (iRun))
                            try:
                                newAreas.append(getAreaGauss(x,thisFluxWithErr,maxFlux,x0,sigma,addOnBothSidesOfX=0.,show=False,save=None)[0])
                            except Exception as e:
                                plt.plot(x,thisFlux,label='original')
                                plt.plot(x,thisFluxWithErr,label='with errors')
                                plt.show()
                                newAreas.append(area)
                                STOP
                        if show:
                            plt.legend()
                            plt.show()
                        newAreas = np.array(newAreas)
                        print('newAreas = ',len(newAreas),': ',newAreas)
                        if show:
                            plt.hist(newAreas)
                            plt.show()
                        sDev = np.std(newAreas)
                        print('sDev = ',sDev)
                        csvLines.setData(keys[iLine]+'e',i,'%.3E' % (sDev))
        csvFree.writeCSVFile(csvLines,csvLinesFileName,'\t')

if __name__ == '__main__':
    spectraDir = '/Users/azuri/spectra/GTC'
    (_, _, filenames) = next(os.walk(spectraDir))
    csvLinesFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/observation.dat'
    hash_fitsFiles = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/fitsfiles.csv')
#    spectrumFileName = '/Users/azuri/spectra/GTC/LDu1_sum.fits'
    for spectrumFileName in filenames:
        print('spectrumFileName[-5:] = <'+spectrumFileName[-5:]+'>')
        if ((spectrumFileName[-5:] == '.fits') 
            and (spectrumFileName != 'strange_blue_star_GT220816.fits') 
            and ('SNR' not in spectrumFileName)
            and (spectrumFileName != 'K1-6a_GT160516.fits')):
            print('starting')
            spectrumFileName = os.path.join(spectraDir,spectrumFileName)
            idPNMain = None
            for i in range(hash_fitsFiles.size()):
                if spectrumFileName[spectrumFileName.rfind('/')+1:] == hash_fitsFiles.getData('fileName',i):
                    idPNMain = hash_fitsFiles.getData('idPNMain',i)
            if idPNMain is None:
                print('ERROR: did not find ',spectrumFileName)
                #STOP
            else:
                errors = calculateErrors(spectrumFileName,idPNMain,csvLinesFileName,show=False)
