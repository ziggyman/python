from astropy.coordinates import EarthLocation
import astropy.io.fits as pyfits
from drUtils import addSuffixToFileName, combine, separateFileList, silentRemove,extractSum
from drUtils import subtractOverscan, subtractBias, cleanCosmic, flatCorrect,interpolateTraceIm
from drUtils import makeSkyFlat, makeMasterFlat, imDivide, extractAndReidentifyARCs, dispCor
from drUtils import readFluxStandardsList,calcResponse,applySensFuncs,extractObjectAndSubtractSky
from drUtils import scombine,continuum,subtractMedianSky,removeFilesFromListWithAngleNotEqualTo
from drUtils import getWavelengthArr,getListOfFiles,getHeaderValue,fixDBSHeaders,writeFits1D,invertY
from drUtils import cleanSpec,setHeaderValue,traceMultiApertureImage,calcProfile
import numpy as np
import os
from drUtils import getImageData, getHeader, getHeaderValue

setting = 'b'

workPath = os.path.join('/Users/azuri/daten/uni/HKU/Pa30/CAHA_3.5m/231016_PMAS',setting)
fNameFileList = os.path.join(workPath,'allFits.list')

overscanSection = '[2073:2109,:]'
if setting == 'b':
    trimSection = '[26:2072,45:990]'
    nPeaksShould = 184
elif setting == 'c':
    trimSection = '[26:2072,104:1028]'
    nPeaksShould = 184
else:
    print('could not identify setting <'+setting+'>')

if False:
    command = 'ls '+os.path.join(workPath,'/pma-*.fits')+' > '+fNameFileList
    os.system(command)
    with open(fNameFileList,'r') as f:
        fileList = f.readlines()
    fileList = [f.strip() for f in fileList]

    for file in fileList:
        naxis1 = int(getHeaderValue(file,'NAXIS1'))
        object = getHeaderValue(file,'OBJECT')
        if (naxis1 > 2000) and (not 'FOCUS' in object):
            newFileName = file[:file.rfind('/')+1]+object+'_'+file[file.rfind('/')+1:]
            command = 'cp '+file.replace(' ','\ ')+' '+newFileName.replace(' ','\ ')
            print('command = ',command)
            os.system(command)

if False:
    command = 'ls '+workPath+'/*??_pma-2023*pmas-?.fits > '+fNameFileList
    os.system(command)
    inputList = getListOfFiles(fNameFileList)
    outputList = [f.replace('.fits','_ot.fits') for f in inputList]
    subtractOverscan(inputList,
                    overscanSection,
                    trimSection=trimSection,
                    fitsFilesOut=outputList,
                    overwrite=True)

# create master bias
masterBias = os.path.join(workPath,'masterBias.fits')

if False:
    bias_ot_list = os.path.join(workPath,'bias_ot.list')
    command = 'ls '+workPath+'/BIAS*ot.fits > '+bias_ot_list
    os.system(command)
    combinedBias = combine(getListOfFiles(os.path.join(workPath,'bias_ot.list')),
                        combinerMethod='median',
                        clippingMethod='sigma',
                        clippingParameters={'niter':0,
                                            'low_thresh':-3.,
                                            'high_thresh':3.,
                                            'func':np.ma.median},
                        scaling=False,
                        fitsOutName=masterBias)
    print('average sigma 0: mean(combinedBias) = ',np.mean(combinedBias))

if False:
    for inputFiles in ['arc*_ot.fits', '?ont*_ot.fits', 'flat*_ot.fits', 'obj*_ot.fits','std*_ot.fits']:
        inputList = os.path.join(workPath,inputFiles[:inputFiles.find('*')]+'.list' if inputFiles[0] != '?' else 'cont.list')
        command = 'ls '+os.path.join(workPath,inputFiles)+' > '+inputList
        os.system(command)
        print('created input list for subtractBias <'+inputList+'>')
        #if 'ont' in inputFiles:
        #    STOP
        subtractBias(getListOfFiles(inputList),
                    masterBias,
                    fitsFilesOut=[f.replace('ot.fits','otz.fits') for f in getListOfFiles(inputList)],
                    overwrite=True)

#trace flats
im = getListOfFiles(os.path.join(workPath,'cont.list'))[1].replace('_ot.fits','_otz.fits')
outFileName = os.path.join(workPath,'database')
if not os.path.exists(outFileName):
    command = 'mkdir '+outFileName
    os.system(command)
dbFileName = os.path.join(outFileName,'ap'+im[im.rfind('/')+1:].replace('.fits',''))
if False:
    traceMultiApertureImage(im,
                            dbFileName,
                            redoApertureNumber=0,
                            startAtApertureNumber=None)

calcProfile(im,
            dbFileName,
            profWidth = 2.5,
            swathWidth = 300)
