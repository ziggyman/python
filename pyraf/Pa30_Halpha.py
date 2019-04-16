import astropy.io.fits as pyfits
import os
import numpy as np
import matplotlib.pyplot as plt
from pyraf import iraf
from pyraf.iraf import noao
#from pyraf.iraf.noao import imred
#from pyraf.iraf.noao.imred import ccdred
#from myUtils import getWavelength, subtractSky, populateObjectArray, boxCarMedianSmooth, markAreas, findFirstIdxWithValGT, getDateTime

xLeft = [799,1045]
xRight = [25,970]
yBoth = [1464,1522]

#imcopy ../temp/0000956344-20160707-OSIRIS-OsirisLongSlitSpectroscopy1[1:59,1:1193] blueprint
template = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/Halpha/blueprint.fits'

objLeft = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/temp/0000956344-20160707-OSIRIS-OsirisLongSlitSpectroscopy1.fits'
objRight = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/temp/0000956344-20160707-OSIRIS-OsirisLongSlitSpectroscopy2.fits'

bias = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/bias/0000956%3d-20160707-OSIRIS-OsirisBias.fits'
biasNumbers = [190,208]

flat = '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/flat/00009562%2d-20160707-OSIRIS-OsirisSpectralFlat.fits'
flatNumbers = [31,35]

def createImage(inFileLeft, hduLeft, inFileRight, hduRight, template, outFile):
    hdulistTemp = pyfits.open(template)
    headerTemp = hdulistTemp[0].header
    dataTemp = hdulistTemp[0].data
    print('dataTemp.shape = ',dataTemp.shape)

    dataInLeft = pyfits.getdata(inFileLeft, hduLeft)
    print('dataInLeft.shape = ',dataInLeft.shape)

    dataInRight = pyfits.getdata(inFileRight, hduRight)
    print('dataInRight.shape = ',dataInRight.shape)

    row = 0
    for x in np.arange(xLeft[0], xLeft[1]+1, 1):
        col = 0
        for y in np.arange(yBoth[0]+2, yBoth[1]+3, 1):
            dataTemp[row,col] = dataInLeft[y,x]
            col += 1
        row += 1
#        print('left: dataTemp[',row-1,',:] = ',dataTemp[row-1,:])

    for x in np.arange(xRight[0], xRight[1]+1, 1):
        col = 0
        for y in np.arange(yBoth[0], yBoth[1]+1, 1):
            dataTemp[row, col] = dataInRight[y,x]
            col += 1
        row += 1
#        print('right: dataTemp[',row-1,',:] = ',dataTemp[row-1,:])

    hdulistTemp[0].data = dataTemp
    sec = '[1:%d,1:%d]' % (dataTemp.shape[1], dataTemp.shape[0])
    hdulistTemp[0].header['CCDSEC'] = sec
    hdulistTemp[0].header['DATASEC'] = sec
    hdulistTemp.writeto(outFile, overwrite=True)
    hdulistTemp.close()
    print(outFile+' written')

def myCombine(inFiles, outFile, fileType):
    noao.imred()
    noao.imred.ccdred()
    fileList = inFiles[0][0:inFiles[0].rfind('/')+1]+fileType+'.list'
    with open(fileList, 'w') as f:
        for inFile in inFiles:
            f.write(inFile+'\n')

    if fileType == 'zero':
        iraf.noao.imred.ccdred.zerocombine(input=','.join(inFiles),
                                           output=outFile,
                                           combine='average',
                                           reject='avsigclip',
                                           ccdtype=fileType,
                                           process='no',
                                           delete='no',
                                           clobber='yes',
                                           scale='none',
                                           statsec='',
                                           nkeep=1,
                                           mclip='yes',
                                           lsigma=3.0,
                                           hsigma=3.0,
                                           rdnoise=4.5,
                                           gain=0.95)

createImage(objLeft, 0, objRight, 0, template, '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/Halpha/object.fits')

biases = []
for iBias in np.arange(biasNumbers[0], biasNumbers[1]+1, 1):
    biases.append('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/Halpha/bias%3d.fits' % (iBias))
    createImage(bias % (iBias), 1, bias % (iBias), 2, template, biases[len(biases)-1])

flats = []
for iFlat in np.arange(flatNumbers[0], flatNumbers[1]+1, 1):
    flats.append('/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/Halpha/flat%2d.fits' % (iFlat))
    createImage(flat % (iFlat), 1, flat % (iFlat), 2, template, flats[len(flats)-1])

#myCombine(biases, '/Volumes/obiwan/azuri/spectra/IPHAS_GTC/Laurence/GTC4-16AMEX/OB0065/Halpha/combinedBias.fits', 'zero')
