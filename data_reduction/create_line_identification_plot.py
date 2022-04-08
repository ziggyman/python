import astropy.io.fits as pyfits
import lineid_plot
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import numpy as np
from drUtils import getWavelengthArr,getImageData,readMyPNLineList
import csvFree,csvData

refArcFile = '/Users/azuri/spectra/saao/saao_may2007/REDUCED/night2/done/wave_calib/190STARWL.fits'
refWLenArr = getWavelengthArr(refArcFile)
print('refWLenArr = ',refWLenArr)

#lineListFile = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_gr7_angle15_70_lines_identified_good_aug2018.dat'
#specFile = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_gr7_15_70_otzxfifEc_aug2018.fits'
#lineListFileRef = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_gr7_angle0_0_lines_identified_good_aug2013.dat'
#refSpecFile = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_otzxfifEc_aug2013.fits'
lineListFileRef = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_may2007.new'
refSpecFile = '/Users/azuri/stella/referenceFiles/spupnic/refArc_spupnic_may2007_d.fits'

lineListFile = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_may2007.temp'
specFile = '/Users/azuri/spectra/saao/saao_may2007/RAW/070511/refArc_spupnic_may2007_Ec.fits'

aug2013Arc = '/Users/azuri/spectra/saao/saao-aug2013/070513/ARC_ARC_a3160074_otzxfifEcd.fits'

with open(lineListFileRef,'r') as f:
    linesRef = f.readlines()
linePosRef = []
lineLabelRef = []
for line in linesRef:
    line = line.strip().split()
    print('line = ',line)
    linePosRef.append(float(line[0]))
    lineLabelRef.append(line[1])


with open(lineListFile,'r') as f:
    lines = f.readlines()
linePos = []
lineLabel = []
for line in lines:
    line = line.strip().split()
    print('line = ',line)
    linePos.append(float(line[0]))
    lineLabel.append(line[1])

refSpec = getImageData(refSpecFile,0)

def plotRefSpec(lam,flux,title='refArc_spupnic_may2007_d.fits'):
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    #ax.plot(wLenArr[0:len(wLenArr)-1],spec)
    ax.plot(lam,flux)

    ax.set_xlabel('Wavelength [$\mathrm{\AA}$]')
    ax.set_ylabel('$\mathrm{F}_\lambda$ [$\mathrm{erg}$ $\mathrm{s^{-1}} \mathrm{cm^{-2}} \mathrm{\AA^{-1}}]$')
    #ax.set_xlim([0,2000])
    yLim = ax.get_ylim()[1]
    ax.set_ylim([0,yLim*1.4])
    lineid_plot.plot_line_ids(lam,flux,[float(i) for i in lineLabelRef],lineLabelRef,ax=ax,arrow_tip=yLim*1.1,box_loc=yLim*1.27)
    plt.title=title
    plt.savefig(refSpecFile[:refSpecFile.rfind('.')]+'.pdf', bbox_inches='tight')
    plt.show()
#    plt.plot(wLenArr,spec)


def plotSpecVsWLen():
    wLen = getWavelengthArr(specFile,0)
    spec = getImageData(specFile,0)
#    wLen = range(len(spec))
    for i in range(len(wLen)):
        print('wLen = ',wLen[i],', spec = ',spec[i])

    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    ax.plot(wLen,spec)
    #plt.show()

    ax.set_xlabel('Wavelength [$\mathrm{\AA}$]')
    ax.set_ylabel('$\mathrm{F}_\lambda$ [$\mathrm{erg}$ $\mathrm{s^{-1}} \mathrm{cm^{-2}} \mathrm{\AA^{-1}}]$')
    #ax.set_xlim([0,2000])
    yLim = 1e6#ax.get_ylim()[1]
    ax.set_ylim([0,yLim*1.3])
    lineid_plot.plot_line_ids(wLen,spec,linePos,lineLabel,ax=ax,arrow_tip=yLim*1.1,box_loc=yLim*1.27)
    plt.title=specFile[specFile.rfind('/')+1:]
    plt.savefig(specFile[:specFile.rfind('.')]+'.pdf', bbox_inches='tight')
    plt.show()

def plotAug2013():
    spec = getImageData(aug2013Arc,0)
    wLen = getWavelengthArr(aug2013Arc,0)
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    ax.plot(wLen,spec)

    ax.set_xlabel('Wavelength [$\mathrm{\AA}$]')
    ax.set_ylabel('$\mathrm{F}_\lambda$ [$\mathrm{erg}$ $\mathrm{s^{-1}} \mathrm{cm^{-2}} \mathrm{\AA^{-1}}]$')
    #ax.set_xlim([0,2000])
    yLim = ax.get_ylim()[1]
    ax.set_ylim([0,yLim*1.4])
    lineid_plot.plot_line_ids(wLen,spec,[float(l) for l in lineLabel],lineLabel,ax=ax,arrow_tip=yLim*1.1,box_loc=yLim*1.27)
    plt.savefig(aug2013Arc[:aug2013Arc.rfind('.')]+'.pdf', bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    if False:
        plotAug2013()
    if True:
        refWLenArr = getWavelengthArr(refSpecFile)
        plotRefSpec(refWLenArr,refSpec,title=refSpecFile[refSpecFile.rfind('/')+1:])
    plotSpecVsWLen()
