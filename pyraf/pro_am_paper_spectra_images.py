import astropy.io.fits as pyfits
import lineid_plot
import os
import numpy as np
import matplotlib.pyplot as plt
from myUtils import getWavelength,getImageData,getHeader
from drUtils import readMyPNLineList



def plotSpectraImage(fitsList,lineList,objectName,outFileName):
    spectra = []
    for fits in fitsList:
        header = getHeader(fits,0)
        wLen = getWavelength(header,1)
        spec = getImageData(fits,0)
        spec = 100. * spec/np.max(spec)
        color = 'r'
        if 'SAAO' in fits:
            color = 'g'
        spectra.append([wLen,spec,color])
    xMin = np.min([np.min(spec[0]) for spec in spectra])
    lineWave,lineLabel = lineList
    fig, ax = plt.subplots()
#    ax.set_xlim([4000,7000])
    ax.set_xlim(xMin,7500.)
    ax.set_ylim([-20,140])
    ak = lineid_plot.initial_annotate_kwargs()
#    ak['arrowprops']['arrowstyle'] = "->"
    pk = lineid_plot.initial_plot_kwargs()
    for spectrum in spectra:
        ax.plot(spectrum[0],spectrum[1],spectrum[2]+'-')
#    pk['color'] = "blue"# if spectra[0][2] == 'r' else 'green'
    print('spectra[0][0] = ',spectra[0][0])
    print('spectra[0][1] = ',spectra[0][1])
    print('lineWave = ',lineWave)
    print('lineLabel = ',lineLabel)
    lineid_plot.plot_line_ids(spectra[0][0],5.+spectra[0][1],lineWave,lineLabel,ax=ax,arrow_tip=110,box_loc=130)#, annotate_kwargs=ak, plot_kwargs=pk)
    for spectrum in spectra:
        ax.plot(spectrum[0],5.+spectrum[1],'w-')
    for spectrum in spectra:
        ax.plot(spectrum[0],spectrum[1],spectrum[2]+'-')
    ax.set_xlabel('wavelength [$\mathrm{\AA}$]')
    ax.set_ylabel('relative intensity')
    ax.set_title(objectName)
    plt.savefig(outFileName, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    lineList = readMyPNLineList('/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/pnLineList_plot.dat')
    print('lineList = ',lineList)
    fitsList = ['/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/FRA_Kn42_CO170620_id10887.fits',
                '/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/SAAO201808_Kn42_SA130818_id10887.fits']
    plotSpectraImage(fitsList,lineList,'Kn 42','/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/Kn42.pdf')
    fitsList = ['/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/FRA_MPAJ1827-1328_CO280520_id2387.fits',
                '/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/SAAO202005_MPA_J1827-1328_SA150520_id2387.fits']
    plotSpectraImage(fitsList,lineList,'MPA J1827-1328','/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/MPAJ1827-1328.pdf')
    fitsList = ['/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/FRA_PM1-322_CO210620_id4808.fits',
                '/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/SAAO201909_PM_1-322_SA090919_id4808.fits']
    plotSpectraImage(fitsList,lineList,'PM 1-322','/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/PM1-322.pdf')
    fitsList = ['/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/FRA_Kn134_CO220220_id31922.fits',
                '/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/SAAO201909_PNG228.6-09.3_SA090919_id31922.fits']
    plotSpectraImage(fitsList,lineList,'Kn 134','/Users/azuri/daten/uni/HKU/publications/pro-am-paper/spectra/Kn134.pdf')
