import astropy.io.fits as pyfits
import lineid_plot
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import numpy as np
from drUtils import getWavelengthArr,getImageData,readMyPNLineList
import csvFree,csvData

refArcFile = '/Users/azuri/spectra/saao/saao_may2007/REDUCED/night2/done/wave_calib/190STARWL.fits'
wLenArr = getWavelengthArr(refArcFile)
print('wLenArr = ',wLenArr)

#lineListFile = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_gr7_angle15_70_lines_identified_good_aug2018.dat'
#specFile = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_gr7_15_70_otzxfifEc_aug2018.fits'
lineListFile = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_gr7_angle0_0_lines_identified_good_aug2013.dat'
specFile = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_otzxfifEc_aug2013.fits'
#specFile = '/Users/azuri/spectra/saao/saao_may2007/RAW/night1/ARC_ARC_a2150061_otzxfifEc.fits'

with open(lineListFile,'r') as f:
    lines = f.readlines()
linePos = []
lineLabel = []
for line in lines:
    line = line.strip().split()
    print('line = ',line)
    linePos.append(float(line[0]))
    lineLabel.append(line[1])

wLen = getWavelengthArr(specFile,0)
spec = getImageData(specFile,0)

def plotRefSpec():
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    ax.plot(wLenArr[0:len(wLenArr)-1],spec)

    ax.set_xlabel('Wavelength [$\mathrm{\AA}$]')
    ax.set_ylabel('$\mathrm{F}_\lambda$ [$\mathrm{erg}$ $\mathrm{s^{-1}} \mathrm{cm^{-2}} \mathrm{\AA^{-1}}]$')
    #ax.set_xlim([0,2000])
    yLim = ax.get_ylim()[1]
    ax.set_ylim([0,yLim*1.3])
    lineid_plot.plot_line_ids(list(reversed(wLen)),list(reversed(spec)),[float(i) for i in list(reversed(lineLabel))],list(reversed(lineLabel)),ax=ax,arrow_tip=yLim*1.1,box_loc=yLim*1.2)
    plt.savefig(specFile[:specFile.rfind('.')]+'.pdf', bbox_inches='tight')
    plt.show()
    plt.plot(wLenArr,spec)


def plotSpecVsWLen():
    for i in range(len(wLen)):
        print('wLen = ',wLen[i],', spec = ',spec[i])

    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    ax.plot(wLen,spec)

    ax.set_xlabel('Wavelength [$\mathrm{\AA}$]')
    ax.set_ylabel('$\mathrm{F}_\lambda$ [$\mathrm{erg}$ $\mathrm{s^{-1}} \mathrm{cm^{-2}} \mathrm{\AA^{-1}}]$')
    #ax.set_xlim([0,2000])
    yLim = 1e6#ax.get_ylim()[1]
    ax.set_ylim([0,yLim*1.3])
    lineid_plot.plot_line_ids(wLen,spec,linePos,lineLabel,ax=ax,arrow_tip=yLim*1.1,box_loc=yLim*1.2)
    plt.savefig(specFile[:specFile.rfind('.')]+'.pdf', bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    plotRefSpec()
