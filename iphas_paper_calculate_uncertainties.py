from fnmatch import fnmatch
import matplotlib.pyplot as plt
import numpy as np
import os
import pyneb as pn
#import subprocess

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
from myUtils import degToDMS,degToHMS,angularDistancePyAsl,dmsToDeg,hmsToDeg
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
                        csvLines.setData(keys[iLine]+'e',i,'%0.2E' % (sDev))
        csvFree.writeCSVFile(csvLines,csvLinesFileName,'\t')

def getPNGsfromHashIDs(csvHashPNMain,hashIDs):
    pngs = []
    for hashID in hashIDs:
        found = csvHashPNMain.find('idPNMain',hashID.strip(),0)[0]
        png = csvHashPNMain.getData('PNG',found)
        print('found = ',found,', png = ',png)
        pngs.append(png)
    #STOP
    return pngs

def createLineIntensityTable(csvLinesFileName,hashPNMainFileName,texFileName,pngsWithGoodSNR):
    csvLines = csvFree.readCSVFile(csvLinesFileName,'\t',False)
    hashIDs = csvLines.getData('NAME')
    print('hashIDs = ',hashIDs)
    csvHashPNMain = csvFree.readCSVFile(hashPNMainFileName)
    pngs = getPNGsfromHashIDs(csvHashPNMain,hashIDs)
    print('pngs = ',pngs)
    pngsSortedIndices = np.argsort([float(png[0:5]) for png in pngs])
    print('pngsSortedIndices = ',pngsSortedIndices)
    with open(texFileName,'w') as f:
        f.write('\\clearpage\n')
        f.write('%\\onecolumn\n')
        f.write('%\\begin{landscape}\n')
        f.write('\\begin{longtable}{ | *{10}{l|} }\n')
        f.write('\\caption{Line intensities for each PN from the 27 PNe spectra with sufficient S/N of which 23 are for the new PNe discovered. "n.d." stands for not detected.}\n')
        f.write('\\\\	\\hline\n')
        f.write('IAU PNG & $\\mathrm{H}_\\alpha$ & $\\mathrm{H}_\\beta$ & [NII] 5755\\AA & [NII] 6548\\AA & [NII] 6584\\AA & [OIII] 4363\\AA & [OIII] 5007\\AA & [SII] 6716\\AA & [SII] 6731\\AA \\\\\n')
        f.write('\\endhead  % header material\n')
        f.write('%\\hline\\endfoot  % footer material\n')
        f.write('\\hline\n')
        for i in range(csvLines.size()):
            if pngs[pngsSortedIndices[i]] in pngsWithGoodSNR:
                f.write(pngs[pngsSortedIndices[i]]+' & '
                    +(('$%0.1E'%(float(csvLines.getData('H1r_6563A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('H1r_6563Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ & ').replace('$0\\pm0$','n.d.')
                    +(('$%0.1E'%(float(csvLines.getData('H1r_4861A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('H1r_4861Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ & ').replace('$0\\pm0$','n.d.')
                    +(('$%0.1E'%(float(csvLines.getData('N2_5755A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('N2_5755Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ & ').replace('$0\\pm0$','n.d.')
                    +(('$%0.1E'%(float(csvLines.getData('N2_6548A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('N2_6548Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ & ').replace('$0\\pm0$','n.d.')
                    +(('$%0.1E'%(float(csvLines.getData('N2_6584A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('N2_6584Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ & ').replace('$0\\pm0$','n.d.')
                    +(('$%0.1E'%(float(csvLines.getData('O3_4363A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('O3_4363Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ & ').replace('$0\\pm0$','n.d.')
                    +(('$%0.1E'%(float(csvLines.getData('O3_5007A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('O3_5007Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ & ').replace('$0\\pm0$','n.d.')
                    +(('$%0.1E'%(float(csvLines.getData('S2_6716A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('S2_6716Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ & ').replace('$0\\pm0$','n.d.')
                    +(('$%0.1E'%(float(csvLines.getData('S2_6731A',pngsSortedIndices[i]).replace('E-','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+('\\pm%0.1E'%(float(csvLines.getData('S2_6731Ae',pngsSortedIndices[i]).replace('E','e')))).replace('0.0E+00','0').replace('E-','\\text{e-}')+'$ ').replace('$0\\pm0$','n.d.')
                    +'\\\\\n')
        f.write('\\hline\n')
        f.write('\\label{tab:line_intensities}\n')
        f.write('\\end{longtable}\n')
        f.write('%\\end{landscape}\n')

def reorderTable1(fNameIn,fNameOut):
    csvTable = csvFree.readCSVFile(fNameIn,'&',False)
    print('csvTable.header = ',csvTable.header)
    csvNewTable = csvData.CSVData()
    tempHeader = [' ' for key in csvTable.header]
    tempHeader[0] = csvTable.header[1]
    tempHeader[1] = csvTable.header[0]
    for i in np.arange(2,len(csvTable.header),1):
        tempHeader[i] = csvTable.header[i]
        print('csvTable.header[',i,'] = <'+csvTable.header[i]+'>')
        print('tempHeader[',i,'] set to <'+tempHeader[i]+'>')
    print('csvTable.header = ',csvTable.header)
    print('tempHeader = ',tempHeader)
    csvNewTable.header = tempHeader
    print('csvNewTable.header = ',csvNewTable.header)
    pngNames = csvTable.getData(' IAU PNG ')
    sortedIndices = np.argsort([float(png.strip()[0:5]) for png in pngNames])
    for i in range(csvTable.size()):
        csvNewTable.append(['' for i in csvTable.header])
        for j in range(len(csvTable.header)):
            csvNewTable.setData(csvNewTable.header[j],i,csvTable.getData(csvNewTable.header[j],sortedIndices[i]))
    csvFree.writeCSVFile(csvNewTable,fNameOut,'&')

def reorderTableElectrons(fNameIn,fNameOut,hashPNMainFileName):
    csvTable = csvFree.readCSVFile(fNameIn,'&',False)
    print('csvTable.header = ',csvTable.header)
    csvNewTable = csvData.CSVData()
    tempHeader = [' ' for key in csvTable.header[1:]]
    tempHeader[0] = 'IAU PNG'
    for i in np.arange(2,len(csvTable.header),1):
        tempHeader[i-1] = csvTable.header[i]
        print('csvTable.header[',i,'] = <'+csvTable.header[i]+'>')
        print('tempHeader[',i-1,'] set to <'+tempHeader[i-1]+'>')
    print('csvTable.header = ',csvTable.header)
    print('tempHeader = ',tempHeader)
    csvNewTable.header = tempHeader
    print('csvNewTable.header = ',csvNewTable.header)
    csvHashPNMain = csvFree.readCSVFile(hashPNMainFileName)
    hashIDs = csvTable.getData(' HASH ID ')
    pngNames = getPNGsfromHashIDs(csvHashPNMain,hashIDs)
    sortedIndices = np.argsort([float(png[0:5]) for png in pngNames])
    for i in range(csvTable.size()):
        csvNewTable.append(['' for i in csvNewTable.header])
        csvNewTable.setData(csvNewTable.header[0],i,pngNames[sortedIndices[i]])
        for j in np.arange(1,len(csvNewTable.header)):
            csvNewTable.setData(csvNewTable.header[j],i,csvTable.getData(csvNewTable.header[j],sortedIndices[i]))
    csvFree.writeCSVFile(csvNewTable,fNameOut,'&')

def redoImageTable(fNameIn,fNameImagesOut,fNameSpectraImagesOut,hashPNMainFileName):
    csvImages = csvFree.readCSVFile(fNameIn,'&',False)
    hashIDs = csvImages.getData(" HASH ID ")
    csvHashPNMain = csvFree.readCSVFile(hashPNMainFileName)
    pngs = getPNGsfromHashIDs(csvHashPNMain,hashIDs)
    print('pngs = ',pngs)
    pngsSortedIndices = np.argsort([float(png[0:5]) for png in pngs])
    print('pngsSortedIndices = ',pngsSortedIndices)

    csvImages.addColumn('PNG',pngs)
    csvImages.sort(pngsSortedIndices)
    spectraTable = csvData.CSVData()
    spectraTable.header = ['specLeft','specRight']
    for i in np.arange(0,csvImages.size(),2):
        spectraTable.append(['',''])
        spectraTable.setData('specLeft',spectraTable.size()-1,csvImages.getData(' spectrum\\\\',i).replace('\\\\','').replace('width=2.3cm,height=19.5mm','width=0.48\\textwidth'))
        if i+1 < csvImages.size():
            spectraTable.setData('specRight',spectraTable.size()-1,csvImages.getData(' spectrum\\\\',i+1).replace('width=2.3cm,height=19.5mm','width=0.48\\textwidth'))
    csvFree.writeCSVFile(spectraTable,fNameSpectraImagesOut,'&')
    print('removing column "'+' spectrum\\\\'+'"')
    print('csvImages.header = ',csvImages.header)
    csvImages.removeColumn(' spectrum\\\\')
    csvImages.removeColumn('        Target Name ')
    tempHead = [' ' for key in csvImages.header]
    tempHead[0] = 'PNG'
    for i in np.arange(1,len(tempHead),1):
        tempHead[i] = csvImages.header[i-1]
    print('tempHead = ',tempHead)
    csvOut = csvData.CSVData()
    csvOut.header = tempHead
    for i in range(csvImages.size()):
        line = [csvImages.getData(key,i) for key in tempHead]
        line[len(line)-1] = line[len(line)-1]+'\\\\\n'
        csvOut.append(line)
    csvFree.writeCSVFile(csvOut,fNameImagesOut,'&')

def getPAandExpTimes(spectraDir,csvHashPNMain,csvHashFitsFiles,tableName,tableNameNew):
    (_, _, filenames) = next(os.walk(spectraDir))
    print('filenames = ',filenames)
    csvTable = csvFree.readCSVFile(tableName,'&',False)
    for i in range(csvTable.size()):
        for j in range(len(csvTable.header)):
            csvTable.setData(csvTable.header[j],i,csvTable.getData(csvTable.header[j],i).strip().replace('\\\\','').replace('IPHASX','IX'))
    csvTable.addColumn('PA [$^\\circ$]')
    csvTable.addColumn('$\\mathrm{t_{exp}}$ [s]')
    for fName in filenames:
        if fName[fName.rfind('.')+1:] == 'fits':
            slitpa = 0.
            try:
                slitpa = float(getHeaderValue(os.path.join(spectraDir,fName),'SLITPA'))
            except:
                print('fName = ',fName,': header keyword SLITPA does not exist')
            print('fName = '+fName+': slitpa = ',slitpa)
            pa = 0.
            try:
                pa = float(getHeaderValue(os.path.join(spectraDir,fName),'PA'))
            except:
                print('fName = ',fName,': header keyword PA does not exist')
            print('fName = '+fName+': pa = ',pa)
            ipa = 0.
            try:
                ipa = float(getHeaderValue(os.path.join(spectraDir,fName),'IPA'))
            except:
                print('fName = ',fName,': header keyword IPA does not exist')
            print('fName = '+fName+': ipa = ',ipa)
            ipa0 = 150.53
            try:
                ipa0 = float(getHeaderValue(os.path.join(spectraDir,fName),'IPA0'))
            except:
                print('fName = ',fName,': header keyword IPA0 does not exist')
            print('fName = '+fName+': ipa0 = ',ipa0)
            print('ipa-ipa0 = ',ipa-ipa0)
            pa = round(ipa-ipa0)
            print('fName = ',fName,': PA set to ',pa)
            idPNMain = getIDPNMainFromFitsFileName(fName,csvHashFitsFiles)
            if idPNMain != 'notFound':
                png = getPNGsfromHashIDs(csvHashPNMain,[idPNMain])[0]
                print('fName = ',fName,': idPNMain = ',idPNMain,': PNG = ',png)
                findIdx = csvTable.find(' IAU PNG ',png,0)[0]
                if findIdx > -1:
                    csvTable.setData('PA [$^\\circ$]',findIdx,str(pa))
                    expTime = getHeaderValue(os.path.join(spectraDir,fName),'EXPTIME')
                    print('fName = ',fName,': EXPTIME = ',type(expTime),': ',expTime)
                    csvTable.setData('$\\mathrm{t_{exp}}$ [s]',findIdx,str(round(expTime))+'\\\\')
    csvFree.writeCSVFile(csvTable,tableNameNew,'&')

def getIDPNMainFromFitsFileName(fName,csvHashFitsFiles):
    if ((fName[-5:] == '.fits')
        and (fName != 'strange_blue_star_GT220816.fits')
        and ('SNR' not in fName)
        and (fName != 'K1-6a_GT160516.fits')
        and ('sub' not in fName)
        and ('_t.fits' not in fName)
        and ('Pa30_GT' not in fName)):
        print('starting')
        idPNMain = None
        for i in range(csvHashFitsFiles.size()):
            if fName == hash_fitsFiles.getData('fileName',i):
                idPNMain = hash_fitsFiles.getData('idPNMain',i)
        if idPNMain is None:
            print('ERROR: did not find ',fName)
            STOP
        return idPNMain
    elif 'Pa30_GT' in fName:
        return '15569'
    else:
        return 'notFound'

def readCSdata(fNameTablea1a,fNameTablea1c,fNameTablea2):
    with open(fNameTablea1a,'r') as f:
        lines = f.readlines()
    csvCSdata1 = csvData.CSVData()
    csvCSdata1.header = ['PNG',
                         'Name',
                         'GaiaEDR3',
                         'Group',
                         'RAdeg',
                         'DEdeg',
                         'Dang',
                         'Gmag',
                         'BP-RP',
                         'AV',
                         '(BP-RP)0']
    for line in lines:
        print('line = ',len(line),': ',line)
        csvCSdata1.append([line[0:15].strip()[4:],
                           line[16:45].strip(),
                           line[46:65].strip(),
                           line[66],
                           line[68:76].strip(),
                           line[77:85].strip(),
                           line[86:91].strip(),
                           line[92:97].strip(),
                           line[98:103].strip(),
                           line[104:110].strip(),
                           line[111:len(line)].strip(),
                          ])
        print('csvCSdata1[',csvCSdata1.size()-1,'] = ',csvCSdata1.getData(csvCSdata1.size()-1))

    with open(fNameTablea1c,'r') as f:
        lines = f.readlines()
    csvCSdata1c = csvData.CSVData()
    csvCSdata1c.header = ['PNG',
                          'Name',
                          'GaiaEDR3',
                          'Group',
                          'RAdeg',
                          'DEdeg',
                          'Dang',
                          'Gmag',
                          'BP-RP',
                          'AV',
                          '(BP-RP)0']
    for line in lines:
        print('line = ',len(line),': ',line)
        csvCSdata1c.append([line[0:15].strip()[4:],
                            line[16:45].strip(),
                            line[46:65].strip(),
                            line[66],
                            line[68:76].strip(),
                            line[77:85].strip(),
                            line[86:91].strip(),
                            line[92:97].strip(),
                            line[98:103].strip(),
                            line[104:110].strip(),
                            line[111:len(line)].strip(),
                           ])
        print('csvCSdata1c[',csvCSdata1c.size()-1,'] = ',csvCSdata1c.getData(csvCSdata1c.size()-1))
    csvFree.writeCSVFile(csvCSdata1c,fNameTablea1c[:fNameTablea1c.rfind('.')]+'.csv')
    csvCSdata1.append(csvCSdata1c.data)
    csvFree.writeCSVFile(csvCSdata1,fNameTablea1a[:fNameTablea1a.rfind('.')]+'.csv')

    with open(fNameTablea2,'r') as f:
        lines = f.readlines()
    csvCSdata2 = csvData.CSVData()
    csvCSdata2.header = ['PNG',
                        'Plx',
                        'ePlx',
                        'Dist',
                        'bDist',
                        'BDist',
                        'Height',
                        'radAng',
                        'radPhys',
                        'RV',
                        'Morph',
                        'SpType']
    for line in lines:
        print('line = ',len(line),': ',line)
        csvCSdata2.append([line[0:15].strip()[4:],
                           line[16:22].strip(),
                           line[23:28].strip(),
                           line[29:34].strip(),
                           line[35:39].strip(),
                           line[40:45].strip(),
                           line[46:50].strip(),
                           line[52:58].strip(),
                           line[59:65].strip(),
                           line[66:71].strip() if len(line)>71 else '',
                           line[72] if len(line) > 72 else '',
                           line[74:len(line)].strip() if len(line) > 75 else '',
                          ])
        print('csvCSdata[',csvCSdata2.size()-1,'] = ',csvCSdata2.getData(csvCSdata2.size()-1))
    csvFree.writeCSVFile(csvCSdata2,fNameTablea2[:fNameTablea2.rfind('.')]+'.csv')
    return [csvCSdata1,csvCSdata2]

def findGTCpngInCSdata(csvCSdata,csvGTCdata,csvCSCoordsHash,CScoordsOutFileName=None):
    nFound = 0
    nBoth = 0
    nCSinGTC = 0
    both = []
    onlyCS = []
    onlyGTC = []
    distances = []
    csCoords = []
    for iGTC in range(csvGTCdata.size()):
        if csvGTCdata.getData('CS',iGTC) == 'y':
            nCSinGTC += 1
        found = csvCSdata.find('PNG',csvGTCdata.getData('IAU PNG',iGTC),0)[0]
        if found >= 0:
            print('found ',csvGTCdata.getData('IAU PNG',iGTC),' at position ',found)
            nFound += 1
            pngFound = csvCSdata.find('PNG',csvGTCdata.getData('IAU PNG',iGTC),0)[0]
            print('found PNG ',csvGTCdata.getData('IAU PNG',iGTC),' in position ',pngFound)
            if csvGTCdata.getData('CS',iGTC) == 'y':
                nBoth += 1
                idPNMain = csvGTCdata.getData('HASH ID',iGTC)
                raGTC = csvCSCoordsHash.getData('CS_RAJ2000',csvCSCoordsHash.find('idPNMain',idPNMain)[0])
                decGTC = csvCSCoordsHash.getData('CS_DECJ2000',csvCSCoordsHash.find('idPNMain',idPNMain)[0])
                if 'Group' in csvCSdata.header:
                    raCS = degToHMS(float(csvCSdata.getData('RAdeg',pngFound)))
                    decCS = degToDMS(float(csvCSdata.getData('DEdeg',pngFound)))
                    distance = angularDistancePyAsl(hmsToDeg(raGTC),dmsToDeg(decGTC),hmsToDeg(raCS),dmsToDeg(decCS))
                both.append([csvGTCdata.getData('IAU PNG',iGTC),
                             raGTC,
                             raCS if 'Group' in csvCSdata.header else ' ',
                             decGTC,
                             decCS if 'Group' in csvCSdata.header else ' ',
                             distance*3600. if 'Group' in csvCSdata.header else ' ',
                            ])
            else:
                onlyCS.append([csvGTCdata.getData('IAU PNG',iGTC),
                            csvCSdata.getData('Group',pngFound) if 'Group' in csvCSdata.header else csvCSdata.getData('Dist',pngFound)])
                print('CS found in CS but not in HASH: ')
            if 'Dist' in csvCSdata.header:
                distances.append([csvGTCdata.getData('IAU PNG',iGTC),
                                    csvCSdata.getData('Dist',pngFound),
                                    ])

        else:
            if csvGTCdata.getData('CS',iGTC) == 'y':
                idPNMain = csvGTCdata.getData('HASH ID',iGTC)
                onlyGTC.append([csvGTCdata.getData('IAU PNG',iGTC),
                                csvCSCoordsHash.getData('CS_RAJ2000',csvCSCoordsHash.find('idPNMain',csvGTCdata.getData('HASH ID',iGTC))[0]),
                                csvCSCoordsHash.getData('CS_DECJ2000',csvCSCoordsHash.find('idPNMain',csvGTCdata.getData('HASH ID',iGTC))[0])])
        if csvGTCdata.getData('CS',iGTC) == 'y':
            idPNMain = csvGTCdata.getData('HASH ID',iGTC)
            csCoords.append(csvCSCoordsHash.getData('CS_RAJ2000',csvCSCoordsHash.find('idPNMain',csvGTCdata.getData('HASH ID',iGTC))[0])+' '+
                            csvCSCoordsHash.getData('CS_DECJ2000',csvCSCoordsHash.find('idPNMain',csvGTCdata.getData('HASH ID',iGTC))[0])+'\n')
    if CScoordsOutFileName is not None:
        with open(CScoordsOutFileName,'w') as f:
            for row in csCoords:
                f.write(row)
    print('found ',nFound,' GTC CS in csvCSdata, nBoth = ',nBoth,', nCSinGTC = ',nCSinGTC)
    print('onlyCS = ',onlyCS)
    print('onlyGTC = ',len(onlyGTC))
    for i in range(len(onlyGTC)):
        print(onlyGTC[i])
    print('distances = ',distances)
    print('both = ',len(both))
    for i in range(len(both)):
        print(both[i])

def getIntegrators(csvTable1,outFile):
    with open(outFile,'w') as f:
        f.write('mkdir IPHAS_GTC_Integrators\n')
        for i in range(csvTable1.size()):
            idPNMain = csvTable1.getData('HASH ID',i)
            print('i = ',i,': idPNMain = ',idPNMain)
            f.write('cp -pr /data/kegs/Integrators/'+idPNMain+' IPHAS_GTC_Integrators/\n')

def readGaiaEDR3VOTable(fNameVOT):
    from astropy.io.votable import parse
    votable = parse(fNameVOT)
    table = votable.get_first_table()#.to_table(use_names_over_ids=True)
    data = table.array
    print('data = ',data)
    print('dir(table) = ',dir(table))
    fields = table.fields
    print('fields = ',fields)
    print('dir(fields) = ',dir(fields))
    csvTable = csvData.CSVData()
    csvTable.header = [field.name for field in fields]
    print('csvTable.header = ',len(csvTable.header),': ',csvTable.header)
    for i in range(len(data)):
        row = [str(dat) for dat in data[i]]
        #print('type(row) = ',type(row))
        #print('row = ',len(row),': ',row)
        csvTable.append(row)
        print('csvTable.size() = ',csvTable.size())
        #STOP
    for i in range(len(csvTable.header)):
        print('csvTable.getData(',csvTable.header[i],',0) = ',csvTable.getData(csvTable.header[i],0))
    target = ''
    thisTargetIdx = []
    removeRows = []
    for i in range(csvTable.size()):
        if target != csvTable.getData('target_id',i):
            # find closest to previous target
            print('old targets = ',thisTargetIdx)
            if len(thisTargetIdx) > 1:
                dists = [float(csvTable.getData('target_separation (deg)', targetIdx)) * 3600. for targetIdx in thisTargetIdx]
                print('dists = ',dists)
                minDist = min(dists)
                for targetIdx in range(len(thisTargetIdx)):
                    if dists[targetIdx] != minDist:
                        removeRows.append(thisTargetIdx[targetIdx])
                        print('removing row ',thisTargetIdx[targetIdx],' with dist ',dists[targetIdx],', minDist = ',minDist)
            print('new target found at i = ',i,': target_id = <'+csvTable.getData('target_id',i)+'>, dist = ',float(csvTable.getData('target_separation (deg)',i))*3600.)
            target = csvTable.getData('target_id',i)
            thisTargetIdx = [i]
        else:
            thisTargetIdx.append(i)
    print('len(removeRows) = ',len(removeRows))
    for i in np.arange(len(removeRows)-1,-1,-1):
        print('i = ',i)
        csvTable.removeRow(removeRows[i])
        print('row ',removeRows[i],' removed')
    csvFree.writeCSVFile(csvTable,fNameVOT+'.csv',',')

def readBailerJonesTable(fNameIn):
    csvOut = csvData.CSVData()
    csvOut.header = ['source_id',
                     'RAdeg',
                     'DEdeg',
                     'rgeo',
                     'b_rgeo',
                     'B_rgeo',
                     'rpgeo',
                     'b_rpgeo',
                     'B_rpgeo',
                     'Flag']
    with open(fNameIn,'r') as f:
        lines = f.readlines()
    for line in lines:
        csvLine = []
        csvLine.append(line[0:19].strip())
        csvLine.append(line[20:35].strip())
        csvLine.append(line[36:51].strip())
        csvLine.append(line[51:66].strip())
        csvLine.append(line[67:81].strip())
        csvLine.append(line[82:96].strip())
        csvLine.append(line[96:111].strip())
        csvLine.append(line[112:126].strip())
        csvLine.append(line[127:141].strip())
        csvLine.append(line[142:].strip())
        print('line = <'+line+'>')
        print('csvLine = ',csvLine)
        csvOut.append(csvLine)
    csvFree.writeCSVFile(csvOut,fNameIn[:fNameIn.rfind('.')]+'.csv')

def plotEDR3Distances(csvFileNameIn):
    csvDists = csvFree.readCSVFile(csvFileNameIn)
    plt.hist([float(num) for num in csvDists.getData('r_med_geo')],bins=20)
    plt.show()

if __name__ == '__main__':
    spectraDir = '/Users/azuri/spectra/GTC'
    csvLinesFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/observation.dat'
    hash_fitsFiles = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/fitsfiles.csv')
    hashPNNamesFileName = ''
#    spectrumFileName = '/Users/azuri/spectra/GTC/LDu1_sum.fits'
    if False:
        (_, _, filenames) = next(os.walk(spectraDir))
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

    if False:
        reorderTable1('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_PNe.tex','/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_PNe_sorted.tex')
        reorderTable1('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_oldPNe.tex','/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_oldPNe_sorted.tex')
    #    reorderTableElectrons('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table_electrons.tex','/Users/azuri/daten/uni/HKU/IPHAS-GTC/table_electrons_sorted.tex','/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_PNMain.csv')
    #    csvData = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table_electrons_sorted.tex','&',False)
    #    pngs=csvData.getData('IAU PNG')
    #    createLineIntensityTable(csvLinesFileName,hashPNMainFileName='/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_PNMain.csv',texFileName='/Users/azuri/daten/uni/HKU/IPHAS-GTC/lineIntensities.tex',pngsWithGoodSNR=pngs)
        redoImageTable(fNameIn='/Users/azuri/daten/uni/HKU/IPHAS-GTC/imageTable1.tex',fNameImagesOut='/Users/azuri/daten/uni/HKU/IPHAS-GTC/imageTable1-spectra.tex',fNameSpectraImagesOut='/Users/azuri/daten/uni/HKU/IPHAS-GTC/spectraTable1.tex',hashPNMainFileName='/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_PNMain.csv')
        redoImageTable(fNameIn='/Users/azuri/daten/uni/HKU/IPHAS-GTC/imageTable2.tex',fNameImagesOut='/Users/azuri/daten/uni/HKU/IPHAS-GTC/imageTable2-spectra.tex',fNameSpectraImagesOut='/Users/azuri/daten/uni/HKU/IPHAS-GTC/spectraTable2.tex',hashPNMainFileName='/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_PNMain.csv')
        getPAandExpTimes(spectraDir,
                         csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_PNMain.csv'),
                         hash_fitsFiles,
                         '/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_oldPNe_sorted.tex',
                         '/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_oldPNe_sorted_withPAandExpTime.tex',)
        getPAandExpTimes(spectraDir,
                         csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_PNMain.csv'),
                         hash_fitsFiles,
                         '/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_PNe_sorted.tex',
                         '/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_PNe_sorted_withPAandExpTime.tex',)

    if True:
        CSdata1aFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51/tablea1.dat'
        CSdata1cFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51/tablea1c.dat'
        CSdata2FileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/CS_in_GAIA_DR3/J_A+A_656_A51/tablea2.dat'
        csvCSdata1,csvCSdata2 = readCSdata(CSdata1aFileName,CSdata1cFileName,CSdata2FileName)
        csvGTCdata = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_PNe_sorted_withPAandExpTime.tex','&',False)
        csvGTCdata2 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_oldPNe_sorted_withPAandExpTime.tex','&',False)
        print('csvGTCdata.size() = ',csvGTCdata.size())
        print('csvGTCdata2.size() = ',csvGTCdata2.size())
        csvGTCdata.append(csvGTCdata2.data)
        print('csvGTCdata.size() = ',csvGTCdata.size())
        headerNew = []
        for i in range(len(csvGTCdata.header)):
            print('setting <'+csvGTCdata.header[i]+'> to <'+csvGTCdata.header[i].strip()+'>')
            headerNew.append(csvGTCdata.header[i].strip())
        csvGTCdata.header = headerNew
        csvCSCoordsHash = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_CSCoords_250522.csv')
        print('csvGTCdata.header = ',csvGTCdata.header)
        if os.path.exists('/Users/azuri/daten/uni/HKU/IPHAS-GTC/IPHAS_GTC_CSCoords.dat'):
            os.remove('/Users/azuri/daten/uni/HKU/IPHAS-GTC/IPHAS_GTC_CSCoords.dat')
        findGTCpngInCSdata(csvCSdata1,csvGTCdata,csvCSCoordsHash,'/Users/azuri/daten/uni/HKU/IPHAS-GTC/IPHAS_GTC_CSCoords.dat')
        #STOP
        print('csvGTCdata.header = ',csvGTCdata.header)
        findGTCpngInCSdata(csvCSdata2,csvGTCdata,csvCSCoordsHash)
        #getIntegrators(csvGTCdata,'/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_getIntegrators')
    if False:
        readGaiaEDR3VOTable('/Users/azuri/daten/uni/HKU/IPHAS-GTC/IPHAS_GTC_CSCoords_resultsGaiaEDR3')
#        readBailerJonesTable('/Users/azuri/daten/uni/HKU/IPHAS-GTC/bailer-jones2021_distsGaiaEDR3/gedr3dis.sam')
    if True:
        plotEDR3Distances('/Users/azuri/daten/uni/HKU/IPHAS-GTC/bailer-jones2021_distsGaiaEDR3/IPHAS_GTC_CS_dists.csv')
