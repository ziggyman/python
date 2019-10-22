import matplotlib.pyplot as plt
import numpy as np
import os

spectraPath = '/Volumes/work/azuri/data/cmfgen/ritter/'
fname = '/Volumes/work/azuri/data/cmfgen/ritter/n19/obs/spectrum.dat'
#fname = '/Volumes/work/azuri/data/cmfgen/ritter/n21/obs/spectrum.dat'
#fname = os.path.join(spectraPath,'n19/obs_osix/spectrum.dat')
#fname = '/Volumes/work/azuri/data/cmfgen/ritter/a2_nd60/spectrum.dat'
#fname = '/Volumes/work/azuri/data/cmfgen/ritter/n5/spectrum.dat'
#fname = '/Volumes/work/azuri/data/cmfgen/ritter/n5/spectrum_unreddened.dat'
#fname = '/Volumes/work/azuri/data/cmfgen/ritter/n12b/spectrum_unreddened.dat'

spectraList = os.path.join(spectraPath,'spectra_Jy.list')

spectrumFileDered = os.path.join(spectraPath,'obs/Pa30_GT080716_cal_sum_cleaned_scaled_dereddened.data')
spectrumFileRed = os.path.join(spectraPath,'obs/Pa30_GT080716_cal_sum_cleaned_scaled.data')

ewDataFName = '/Volumes/work/azuri/data/cmfgen/ritter/n19/obs_ew/ewdata_fin'

xRange = [3000.,9000]
plotLog = True

# return number of columns and columns
def countCols(str):
    nCols = 1
    tempStr = str
    cols = []
    while tempStr.find(' ') >= 0:
        nCols += 1
        cols.append(tempStr[0:tempStr.find(' ')])
        tempStr = tempStr[tempStr.find(' '):].strip(' ')
    cols.append(tempStr)
    return nCols,cols

def readcmfFile(fname):
    dataArr = []
    nColsOld = 0
    with open(fname) as f:
        for lineStr in f:
            line = lineStr.strip(' \n')
            nCols, cols = countCols(line)
#            print('line = <'+line+'>: nCols = ',nCols,': cols = ',cols)
            if nCols != nColsOld:
                dataArr = []
                dataArrRow = []
                for dat in cols:
                    dataArrRow.append(float(dat))
                dataArr.append(dataArrRow)
                nColsOld = nCols
            else:
                dataArrRow = []
                for dat in cols:
                    dataArrRow.append(float(dat))
                dataArr.append(dataArrRow)
    dataArr = np.array(dataArr)
    print('dataArr.shape = ',dataArr.shape,': ',dataArr)
    return dataArr

def readEwData(fName):
    header = []
    dataLines = []
    with open(fName,'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.rstrip('\n')
            if line[0] != '!':
                keys = line.split(' ')
                if len(header) == 0:
                    print('keys = ',keys)
                    for key in keys:
                        if key != '':
                            if key == keys[len(keys)-2]:
                                header.append(key+keys[len(keys)-1])
                            elif key == keys[len(keys)-1]:
                                #do nothing
                                print('')
                            else:
                                header.append(key)
                    print('header = ',header)
                else:
                    dataLine = {}
                    iHeaderKey = 0
                    for key in keys:
                        if key != '':
                            dataLine.update({header[iHeaderKey]:key})
                            iHeaderKey += 1
                    print('dataLine = ',dataLine)
                    dataLines.append(dataLine)
    print('dataLines = ',dataLines)
    return dataLines

def getDataValueAt(x,y,xWhere):
    i = 0
    while True:
        if x[i] > xWhere:
            yWhere = y[i-1] + (xWhere - x[i-1]) * (y[i]-y[i-1]) / (x[i]-x[i-1])
            print('x,y[i-1] = [',x[i-1],',',y[i-1],'], x,y[i] = [',x[i],',',y[i],']: xWhere = ',xWhere,': yWhere = yWhere')
            return yWhere
        i += 1
        if i >= len(x):
            print('ERROR: xWhere=',xWhere,' not found in x = ',x)
            break

ewData = readEwData(ewDataFName)
print('ewData = ',ewData)

#with open(spectraList,'r') as f:
#    for fname in f:


if True:
    if True:
        if fname[0] != '#':
            dataArr = readcmfFile(os.path.join(spectraPath,fname.rstrip('\n')))#np.loadtxt(os.path.join(spectraPath,fname.rstrip('\n')))
            print('dataArr.shape = ',dataArr.shape)
            fig,ax = plt.subplots()
            print('ax = ',ax)
    #        print('dir(ax) = ',dir(ax))
            maxFlux = 0.
            minFlux = 0.
            maxFluxOld = 0.
            minFluxOld = 0.
            sym = 'g-'
            title=fname.rstrip('\n')
            for i in np.arange(0,dataArr.shape[1],2):
                print('i = ',i)
                #wavelength,flux = np.loadtxt(os.path.join(spectraPath,fname.rstrip('\n')),unpack=True)
                wavelength = dataArr[np.where(dataArr[:,i] != 0.0),i]
                flux = dataArr[np.where(dataArr[:,i] != 0.0),i+1]

                if fname[fname.rfind('/'):].find('ebv0.0') > 0:
                    spectrumFile = spectrumFileDered
                    yMin = 1.8e-3
                else:
                    spectrumFile = spectrumFileRed
                    yMin = 6.e-4
                if i > 1:
                    sym = 'b-'
                print('title = ',title)
                print('wavelength = ',wavelength.shape,': ',wavelength)
                print('flux = ',flux.shape,': ',flux)
    #            ax.plot(wavelength,flux,sym)
    #            flux = 2.99792458E-05 * flux / (wavelength * wavelength * wavelength)
                ax.plot(wavelength[0,:],flux[0,:],sym)
    #            plt.show()
                maxFlux = np.max(flux[np.where((wavelength > 3000.) & (wavelength < 9000.))])
                minFlux = np.min(flux[np.where((wavelength > 3000.) & (wavelength < 9000.) & (flux > 0.0))])
                print('plotted wavelength and flux: maxFlux = ',maxFlux,', minFlux = ',minFlux)
                maxFlux = np.max([1.1 * maxFlux, maxFluxOld])
                minFlux = np.min([0.9 * minFlux, minFluxOld])
                print('plotted wavelength and flux: maxFlux = ',maxFlux,', minFlux = ',minFlux)

                maxFluxOld = maxFlux
                minFluxOld = minFlux

            if dataArr.shape[1] < 3:
                wavePa,fluxPa = np.loadtxt(spectrumFile, unpack=True)
                print('read spectrumFile <'+spectrumFile+'>')

                fluxPa = 3.33564095E+04 * fluxPa * wavePa * wavePa
                maxFluxPa = 1.1*np.max(fluxPa)
                minFluxPa = 0.9*np.min(fluxPa[np.where(fluxPa > 0.0)])
                maxFlux = np.max([maxFlux,maxFluxPa])
                minFlux = np.min([minFlux,minFluxPa])
                print('plotted wavelPa and fluxPa: minFlux = ',minFlux,', maxFlux = ',maxFlux)
                scalingFac = 1.0#maxFlux / maxFluxPa
                print('scalingFac = ',scalingFac)
    #            if scalingFac < 10.0:
    #                scalingFac = 1.0
                #fluxPa =
                ax.plot(wavePa,fluxPa * scalingFac, 'c-')
                print('plotted wavePa and fluxPa: maxFluxPa = ',maxFluxPa,': ylim = ',ax.get_ylim())
                ax.set_title(fname+': scaling factor = %.2f' % (scalingFac))
            else:
                ax.set_title(fname)
            ax.set_xlim(xRange)
#            ax.set_ylim(minFlux,maxFlux)
            ax.set_ylim(yMin,maxFlux)
            print('minFlux = ',minFlux,', maxFlux = ',maxFlux,': ylim = ',ax.get_ylim())
            ax.set_xlabel('Wavelength [$\AA$]')
            ax.set_ylabel('$F_\lambda$ [Jy]')
            if plotLog:
                ax.set_yscale('log')
    #        ax.set_ylabel('$F_\lambda$ $ergs/s/cm^2/\AA$')

            """show ions"""
            for ew in ewData:
                if float(ew['AEW(Ang)']) > 10.:
                    x = float(ew['Lam(Ang)'])
                    y = getDataValueAt(wavelength[0,:],flux[0,:],x)
                    yRange = ax.get_ylim()
                    print('yRange = ',yRange)
                    ax.plot([x,x],[y,yRange[0]+((yRange[1]-yRange[0])/2.)])
                    ax.annotate(ew['Trans.Name'],xy=(x,yRange[0]+((yRange[1]-yRange[0])/2.)),rotation=90)

            plt.show()

