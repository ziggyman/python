#import glob
import csv
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import math

from astropy import units as u
#from astropy.coordinates import SkyCoord
from astropy.stats import circstats
import pycircstat
import subprocess

from PIL import Image, ImageDraw

import csvData
import csvFree
#from galaxyMath import raDecToLB, parallaxToDistance,degToRad,radToDeg
from hammer import Pixel,XY,LonLat,Hammer
#from myUtils import getStarWithMinDist

reesZijlstra = False
mockSample = False
mockRandomSample = False

rotateBipolars = False

path = '/Users/azuri/daten/uni/HKU/PN alignment'
outPath = '/Volumes/discovery/azuri/data/pnAlignment'
if rotateBipolars:
    outPath = outPath + '/bipolars_rotated'
xRangeMax = [-2.83,2.83]
yRangeMax = [-1.4143,1.4143]

if reesZijlstra:
    #Rees+Zijlstra
    dataFileName = os.path.join(path, 'Rees_Zijlstra_table_with_HASH-ID.csv')
    hashFileName = os.path.join(path, 'HASH_true_PNe+004.2-05.9+005.9-02.6.csv')
    #bulge:
    areas = [[-0.172482,0.172482],[-0.174476,0.174476]]

    nBins = 40

    mainClasses = [['B','E','I','S']]
    latexFileName = os.path.join(path, 'PN-alignments_Rees+Zijlstra.tex')
    #gaiaFileNameRoot = '/Volumes/work/azuri/data/gaia/dr2/xy/GaiaSource_%.6f-%.6f_%.6f-%.6f.csv'
    imOutPath = os.path.join(outPath,'/Rees+Zijlstra/')

else:
    #my data:
    dataFileName = os.path.join(path, 'PN-alignments_finished.csv')
#    hashFileName = os.path.join(path, 'HASH_bipolar+elliptical_true_PNe.csv')
    if mockSample:
        hashFileName = os.path.join(path, 'mock/HASH_bipolar+elliptical_true_PNe_withPA_Bmean45_sdev10_Emean135_sdev20.csv')
#        hashFileName = os.path.join(path, 'mock/HASH_bipolar+elliptical_true_PNe_withPA_Bmean15_sdev10_Emean135_sdev20.csv')
    elif mockRandomSample:
        hashFileName = os.path.join(path, 'mock/HASH_bipolar+elliptical_true_PNe_withPA_random.csv')
    else:
        hashFileName = os.path.join(path, 'HASH_bipolar+elliptical_true+likely_PNe+004.2-05.9+005.9-02.6+001.4+05.3+357.5+03.2+357.6-03.3.csv')#HASH_bipolar+elliptical_true_PNe.csv')

    #xRange = [-0.263,0.263]
    #yRange = [-0.0873,0.0873]#yRangeMax#
    areas = [[xRangeMax,yRangeMax],                       # all sky
             #[[0.,xRangeMax[1]],[0.,yRangeMax[1]]],       # all sky 1st quadrant
             #[[xRangeMax[0],0],[0.,yRangeMax[1]]],        # all sky 2nd quadrant
             #[[xRangeMax[0],0],[yRangeMax[0],0.]],        # all sky 3rd quadrant
             #[[0.,xRangeMax[1]],[yRangeMax[0],0.]],       # all sky 4th quadrant
             #[[-0.172482,0.172482],[-0.174476,0.174476]], # bulge
             #[[0.,0.172482],[0.,0.174476]],               # bulge 1st quadrant
             #[[-0.172482,0.],[0.,0.174476]],              # bulge 2nd quadrant
             #[[-0.172482,0.],[-0.174476,0.]],             # bulge 3rd quadrant
             #[[0.,0.172482],[-0.174476,0.]],              # bulge 4th quadrant
            ]
    nBins = 36

    mainClasses = [['E'],['B'],['B','E']]
    latexFileName = os.path.join(outPath, 'PN-alignments.tex')
    #gaiaFileNameRoot = '/Volumes/work/azuri/data/gaia/dr2/xy/GaiaSource_%.6f-%.6f_%.6f-%.6f.csv'


#data = csvFree.readCSVFile(dataFileName)
#flags = data.getData('flag')
#flags = np.array([int(a) for a in flags])
#wh = np.where(flags==1,True,False)
#ones = flags[wh]
#print('found ',ones.shape[0],' flag 1')
#wh = np.where(flags==2,True,False)
#twos = flags[wh]
#print('found ',twos.shape[0],' flag 2')
#wh = np.where(flags==3,True,False)
#threes = flags[wh]
#print('found ',threes.shape[0],' flag 3')
#hashData = csvFree.readCSVFile(hashFileName)

show = False
# Turn interactive plotting off
#plt.ioff()

ham = Hammer()
pixels = ham.getPixels()

#def findInHash(hashData, hashID):
#    ids = hashData.getData('idPNMain')
#    for iLine in np.arange(0,hashData.size(),1):
#        if ids[iLine] == hashID:
#            return iLine
#    return -1

def findPNeWithPAinHASH(inFileNamePAs, inFileNameHASH, outFileName):

    print('csvPAs.inFileNamePAs = <'+inFileNamePAs+'>')
    csvPAs = csvFree.readCSVFile(inFileNamePAs)
    print('csvPAs.header = ',csvPAs.header)
    print('csvPAs.size() = ',csvPAs.size())

    print('csvHASH.inFileNameHASH = <'+inFileNameHASH+'>')
    csvHASH = csvFree.readCSVFile(inFileNameHASH)
    print('csvHASH.header = ',csvHASH.header)
    print('csvHASH.size() = ',csvHASH.size())

    idsWithPA = csvPAs.getData('HASH ID')
    allIDs = csvHASH.getData('idPNMain')

    csvOut = csvData.CSVData()
    csvOut.header = csvHASH.header
    print('csvOut.header = ',csvOut.header)
    csvOut.addColumn('flag')
    csvOut.addColumn('source')
    print('csvOut.header = ',csvOut.header)
    print('csvOut.size() = ',csvOut.size())

    nPA = 0
    for iPA in np.arange(0,csvPAs.size(),1):
        if csvPAs.getData('EPA',iPA) != '':
            print('searching for idsWithPA[iPA]=',idsWithPA[iPA])
            found = False
            previouslyFound = False
            for iHASH in np.arange(0,csvHASH.size(),1):
                if (not found) and (idsWithPA[iPA] == allIDs[iHASH]):
                    if previouslyFound:
                        print('ERROR: previouslyFound')
                        STOP
                    newLine = csvHASH.getData(iHASH)
                    if len(newLine) < len(csvOut.header):
                        newLine.append(csvPAs.getData('flag',iPA))
                        newLine.append(csvPAs.getData('source',iPA))
    #                print('newLine = ',len(newLine),': ',newLine)
    #                print('len(csvOut.header) = ',len(csvOut.header))
                    csvOut.append(newLine)
                    csvOut.setData('EPA',csvOut.size()-1,csvPAs.getData('EPA',iPA))
                    csvOut.setData('GPA',csvOut.size()-1,csvPAs.getData('GPA',iPA))
                    if len(newLine) == len(csvOut.header):
                        csvOut.setData('flag',csvOut.size()-1,csvPAs.getData('flag',iPA))
                        csvOut.setData('source',csvOut.size()-1,csvPAs.getData('source',iPA))

    #                print('iPA = ',iPA,', csvOut.size()-1 = ',csvOut.size()-1)
    #                print('csvPAs.getData(',iPA,') = ',csvPAs.getData(iPA))
    #                print('csvHASH.getData(',iHASH,') = ',csvHASH.getData(iHASH))
    #                print('csvOut.getData(',iPA,') = ',csvOut.getData(iPA))
                    if nPA != csvOut.size()-1:
                        print('nPA(=',nPA,') != csvOut.size()-1(=',csvOut.size()-1,')')
                        print('csvPAs.getData(',iPA,') = ',csvPAs.getData(iPA))
                        print('csvOut.getData(',csvOut.size()-1,') = ',csvOut.getData(csvOut.size()-1))
                        STOP
                    if csvPAs.getData('GPA',iPA) != csvOut.getData('GPA',nPA):
                        print('ERROR with GPAs at iPA = ',iPA,', nPA = ',nPA)
                        print("csvPAs.getData('GPA',",iPA,") = ",csvPAs.getData('GPA',iPA),", csvOut.getData('GPA',",nPA,") = ",csvOut.getData('GPA',nPA))
                        STOP
                    if csvPAs.getData('EPA',iPA) != csvOut.getData('EPA',nPA):
                        print('ERROR with EPAs at iPA = ',iPA,', nPA = ',nPA)
                        print("csvPAs.getData('EPA',",nPA,") = ",csvPAs.getData('EPA',nPA),", csvOut.getData('EPA',",nPA,") = ",csvOut.getData('EPA',nPA))
                        STOP
                    if csvPAs.getData('flag',iPA) != csvOut.getData('flag',nPA):
                        print('ERROR with flags at iPA = ',iPA,', nPA = ',nPA)
                        print("csvPAs.getData('flag',",iPA,") = ",csvPAs.getData('flag',iPA),", csvOut.getData('flag',",nPA,") = ",csvOut.getData('flag',nPA))
                        STOP
                    if csvPAs.getData('source',iPA) != csvOut.getData('source',nPA):
                        print('ERROR with sources at iPA = ',iPA,', nPA = ',nPA)
                        print("csvPAs.getData('source',",iPA,") = ",csvPAs.getData('source',iPA),", csvOut.getData('source',",nPA,") = ",csvOut.getData('source',nPA))
                        STOP

                    found = True
                    previouslyFound = True
                    nPA += 1
                    print(' ')
            if not found:
                print('ERROR: ID ',idsWithPA[iPA],' not found in HASH data')
                STOP

    csvFree.writeCSVFile(csvOut,outFileName)


# @brief Calculate first 4 circular moments of angles
# @param angles: np.array of angles in degrees
# @return: [[first moment direction (degrees), first moment length],[second moment direction, second moment length],...]
def calcMoments(angles):
    anglesDouble = np.array([math.radians(angle * 2.0) for angle in angles])
#    print('anglesDouble = ',anglesDouble)
#    print('first moment = ',circstats.circmoment(anglesDouble))
    moments = np.array([circstats.circmoment(anglesDouble,p=1),
                        circstats.circmoment(anglesDouble,p=2),
                        circstats.circmoment(anglesDouble,p=3),
                        circstats.circmoment(anglesDouble,p=4)]) / 2.
#    mean = circstats.mean(anglesDouble)
#    disp = circstats.disp(anglesDouble)
#    print('mean in degrees = ',math.degrees(mean))
#    print('disp in degrees = ',math.degrees(disp))
    return np.array([[math.degrees(moment[0]), moment[1]] for moment in moments])

# @param angles: np.array of angles of the half circle in degrees
def calcMean(angles):
    anglesDouble = np.array([math.radians(angle * 2.0) for angle in angles])
    r = np.zeros(2)
    for i in np.arange(0,anglesDouble.shape[0],1):
        r[0] += np.cos(anglesDouble[i])
        r[1] += np.sin(anglesDouble[i])
    print('r = ',r)
    lengthR = np.sqrt(r[0]**2 + r[1]**2)
    print('lengthR = ',lengthR)
    alpha = math.degrees(np.arccos(r[0]/lengthR))/2.
    beta = math.degrees(np.arcsin(r[1]/lengthR))/2.
    gamma = math.degrees(np.arctan2(r[1],r[0]))
    if gamma < 0.:
        gamma = 360. + gamma
    gamma = gamma / 2.
    print('alpha = ',alpha,', beta = ',beta,', gamma = ',gamma)
    return np.array([gamma,lengthR])

def plotHammerProjection(x,y,GPA,oneFlags,eFlags,bFlags,fNameOut=None):
    fig = plt.figure(figsize=(25,10))
    plt.axis('off')
    plt.scatter(x[oneFlags][eFlags],y[oneFlags][eFlags],c=GPA[oneFlags][eFlags],marker='o',s=20,cmap='viridis')
    plt.scatter(x[oneFlags][bFlags],y[oneFlags][bFlags],c=GPA[oneFlags][bFlags],marker='d',s=20,cmap='viridis')
    plotLBMarks(10)
    cbarTicks = np.arange(0.,1.0001,20./180.)
    cbar = plt.colorbar(ticks=cbarTicks)
    cbar.set_label('Galactic Position Angle')
    cbarTicks = cbarTicks*180.
    cbarTicks = np.round(cbarTicks)
    cbarTicks = [int(x) for x in cbarTicks]
    cbar.ax.set_yticklabels(cbarTicks)
    cbar.ax.tick_params(labelsize=26)
    ax = cbar.ax
    text = ax.yaxis.label
    font = matplotlib.font_manager.FontProperties(size=26)#family='times new roman', style='italic',
    text.set_font_properties(font)
    plt.tick_params(axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom=False,      # ticks along the bottom edge are off
                    top=False,         # ticks along the top edge are off
                    labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(axis='y',          # changes apply to the y-axis
                    which='both',      # both major and minor ticks are affected
                    left=False,      # ticks along the bottom edge are off
                    right=False,         # ticks along the top edge are off
                    labelleft=False) # labels along the bottom edge are off
    if (xRange[0] != xRangeMax[0]) or (xRange[1] != xRangeMax[1]) or (yRange[0] != yRangeMax[0]) or (yRange[1] != yRangeMax[1]):
        plt.plot([xRange[0], xRange[0]], [yRange[0], yRange[1]], 'b-')
        plt.plot([xRange[1], xRange[1]], [yRange[0], yRange[1]], 'b-')
        plt.plot([xRange[0], xRange[1]], [yRange[0], yRange[0]], 'b-')
        plt.plot([xRange[0], xRange[1]], [yRange[1], yRange[1]], 'b-')
    fig.tight_layout()
    if fNameOut is not None:
        print('plotHammerProjection: fNameOut = <'+fNameOut+'>')
        fNameOutTemp = fNameOut[:-3]+'.tmp.pdf'
        fig.savefig(fNameOutTemp, bbox_inches='tight')
        subprocess.run(["gs","-sDEVICE=pdfwrite","-dCompatibilityLevel=1.4","-dPDFSETTINGS=/ebook","-dNOPAUSE", "-dQUIET", "-dBATCH", "-sOutputFile="+fNameOut, fNameOutTemp])
        subprocess.run(["rm",fNameOutTemp])
    if show:
        plt.show()
    else:
        plt.close(fig)
    print('after all_pne: ',plt.get_fignums(),' figures open')


# FROM https://stackoverflow.com/questions/22562364/circular-histogram-for-python
#@param: angles: angles over full circle in degrees
#@param: mean: mean of angles over full circle in degrees
#@param: sigma: deviation of angles over full circle in degrees
def rose_plot(ax, angles, bins=16, density=None, offset=0, lab_unit="degrees",
              start_zero=False, fill=False, color='white', max_count=None,
              max_size=None, smooth=False, mean=None, sigma=None, **param_dict):
    """
    Plot polar histogram of angles on ax. ax must have been created using
    subplot_kw=dict(projection='polar').
    """
    # Wrap angles to [-pi, pi)
    radians = np.array([math.radians(angle) for angle in angles])
    radians = (radians + np.pi) % (2*np.pi) - np.pi

    # Set bins symetrically around zero
    if start_zero:
        # To have a bin edge at zero use an even number of bins
        if bins % 2:
            bins += 1
        bins = np.linspace(-np.pi, np.pi, num=bins+1)

    # Bin data and record counts
    count, bin = np.histogram(radians, bins=bins)
    if smooth:
        smoothedCount = np.zeros(len(count))
        for i in np.arange(0,len(count),1):
            if i == len(count)-1:
                smoothedCount[i] = (count[i-1]+count[i]+count[0]) / 3.
            else:
                smoothedCount[i] = (count[i-1]+count[i]+count[i+1]) / 3.
        count = smoothedCount
    maxCount = np.max(count)

    # Compute width of each bin
    widths = np.diff(bin)

    # By default plot density (frequency potentially misleading)
#    maxArea = None
    if density is None or density is True:
        # Area to assign each bin
        if max_size and max_count:
            area = count / max_size
        else:
            area = count / radians.size
        # Calculate corresponding bin radius
        radius = (area / np.pi)**.5
    else:
        radius = count
    ax.bar(bin[:-1], radius, zorder=1, align='edge', width=widths,
           edgecolor='C0', fill=fill, linewidth=1, color=color)
    if max_count and max_size:
        if density is None or density is True:
            max_area = max_count / max_size
            max_radius = (max_area / np.pi)**.5
        else:
            max_radius = max_count
        ax.bar(0,max_radius,width=0.001)

    if mean is not None:
        maxR = np.max(radius)*1.05
        ax.bar([math.radians(mean),math.radians(mean+180)],[maxR,maxR],width=0.01,linewidth=5,color='red')
        if sigma is not None:
            maxR = np.max(radius)*1.05
            ax.bar([math.radians(mean+sigma),math.radians(mean+sigma+180)],[maxR,maxR],width=0.01,linewidth=5,color='black')
            ax.bar([math.radians(mean-sigma),math.radians(mean-sigma+180)],[maxR,maxR],width=0.01,linewidth=5,color='black')

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels, they are mostly obstructive and not informative
    ax.set_yticks([])

    if lab_unit == "radians":
        label = ['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$',
                  r'$\pi$', r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$']
        ax.set_xticklabels(label)
    return maxCount,radians.size

#@brief construct a vector diagram
#@param angles: angles over full circle in degrees
def vectorDiagram(angles, fNameOut=None):
    fig = plt.figure(figsize=(5,4.7))
    plt.axis('off')
    a = np.zeros(2)
    b = np.zeros(2)
    radius = []
    for angle in angles:
        radian = math.radians(angle)
        x = [a[0],a[0]-np.sin(radian)]
        y = [a[1],a[1]+np.cos(radian)]
        plt.plot(x,y)
        a[0] = x[1]
        a[1] = y[1]
#        print('vectorDiagram: x[1] = ',x[1],', y[1] = ',y[1],': a = ',a)

#        m = [b[0],b[0]+np.cos(radian)]
#        n = [b[1],b[1]+np.sin(radian)]
#        plt.plot(m,n)
#        b[0] = m[1]
#        b[1] = n[1]
#        print('vectorDiagram: m[1] = ',m[1],', n[1] = ',n[1],': b = ',b)
#        print(' ')
        radius.append(np.sqrt(a[0]**2 + a[1]**2))
    maxRadius = np.amax(radius)
    print('vectorDiagram: a = ',a)
#    print('vectorDiagram: b = ',b)

    radius = np.sqrt(a[0]**2 + a[1]**2)
    alpha = np.degrees(np.arcsin(-a[0]/radius))
    beta = np.degrees(np.arccos(a[1]/radius))
    gamma = np.degrees(np.arctan2(a[1],a[0]))-90.
    if gamma < 0:
        gamma += 360.
    print('vectorDiagram: in degrees: radius = ',radius,', alpha = ',alpha,', beta = ',beta,', gamma = ',gamma)

#    radiusb = np.sqrt(b[0]**2 + b[1]**2)
#    alphab = np.degrees(np.arcsin(b[0]/radius))
#    betab = np.degrees(np.arccos(b[1]/radius))
#    gammab = np.degrees(np.arctan2(b[1],b[0]))
#    print('vectorDiagram: in degrees: radiusb = ',radiusb,', alphab = ',alphab,', betab = ',betab,', gamma = ',gammab)

    plt.plot([0,a[0]],[0,a[1]],'r-',linewidth=3)
    plt.plot([0,b[0]],[0,b[1]],'b-',linewidth=3)
    circle1=plt.Circle((0,0),maxRadius,color='b', fill=False)
    plt.plot([0,0],[-maxRadius,maxRadius],'grey',linewidth=0.5)
    plt.plot([-maxRadius,maxRadius],[0,0],'grey',linewidth=0.5)
    plt.plot([-maxRadius/np.sqrt(2),maxRadius/np.sqrt(2)],[-maxRadius/np.sqrt(2),maxRadius/np.sqrt(2)],'grey',linewidth=0.5)
    plt.plot([-maxRadius/np.sqrt(2),maxRadius/np.sqrt(2)],[maxRadius/np.sqrt(2),-maxRadius/np.sqrt(2)],'grey',linewidth=0.5)
    plt.gcf().gca().add_artist(circle1)
    plt.xlim(-maxRadius, maxRadius)
    plt.ylim(-maxRadius, maxRadius)
    plt.text(-maxRadius/50.,maxRadius+maxRadius/100.,'$0^\degree$')
    plt.text(-maxRadius/15.,-maxRadius-maxRadius/13.,'$180^\degree$')
    plt.text(-maxRadius-maxRadius/9.,-maxRadius/40.,'$90^\degree$')
    plt.text(maxRadius+maxRadius/100.,-maxRadius/40.,'$270^\degree$')
    plt.text(-maxRadius/np.sqrt(2.)-maxRadius/11,maxRadius/np.sqrt(2.),'$45^\degree$')
    plt.text(maxRadius/np.sqrt(2.)+maxRadius/100,maxRadius/np.sqrt(2.),'$315^\degree$')
    plt.text(-maxRadius/np.sqrt(2.)-maxRadius/5,-maxRadius/np.sqrt(2.)-maxRadius/20,'$135^\degree$')
    plt.text(maxRadius/np.sqrt(2.)+maxRadius/30,-maxRadius/np.sqrt(2.)-maxRadius/20,'$225^\degree$')

    fig.tight_layout()
    if fNameOut is not None:
        fig.savefig(fNameOut, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)
    return [gamma,radius]

#@param angles: angles over full circle in degrees
def linearOrderDiagram(angles, fNameOut=None):
    print('linearOrderDiagram: angles = ',angles)
    rad = np.array([math.radians(angle) for angle in angles])
    linearOrderX = np.sort(rad/(2.*np.pi))
    print('linearOrderDiagram: linearOrderX = ',linearOrderX)
#    STOP
    linearOrderY = np.zeros(linearOrderX.shape)
    for iX in np.arange(0,linearOrderX.shape[0],1):
        linearOrderY[iX] = iX / (linearOrderX.shape[0]+1)

    fig = plt.figure(figsize=(5,5))
    plt.scatter(linearOrderY,linearOrderX)
    plt.plot([0,1],[0,1])
    plt.xlabel('Uniform quantiles',fontsize=14)
    plt.ylabel('Sample quantiles',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    fig.tight_layout()
    if fNameOut is not None:
        fig.savefig(fNameOut, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)

# @brief return input array lbxyGPA with values inside [x0,x1], [y0,y1]
def selectXY(lbxyGPA_in, x0, x1, y0, y1):
    x = np.array([i[2] for i in lbxyGPA_in])
    y = np.array([i[3] for i in lbxyGPA_in])

    lbxyGPA_out = []
    for i in np.arange(0,len(lbxyGPA_in),1):
        if (x[i] >= x0) and (x[i] <= x1) and (y[i] >= y0) and (y[i] <= y1):
            lbxyGPA_out.append(lbxyGPA_in[i])
    return lbxyGPA_out

# @brief plot l and b every x degrees
def plotLBMarks(x):
    lArr = np.arange(0,360.1,0.1)
    bArr = np.arange(-90, 90.1, 0.1)
    xArr = []
    yArr = []
    for l in lArr:
        for b in np.arange(-90,91,x):
            xy = ham.lonLatToXY(l,b)
            xArr.append(xy.x)
            yArr.append(xy.y)
    for b in bArr:
        for l in np.arange(0,361,x):
            xy = ham.lonLatToXY(l,b)
            xArr.append(xy.x)
            yArr.append(xy.y)
        l=180.0001
        xy = ham.lonLatToXY(l,b)
        xArr.append(xy.x)
        yArr.append(xy.y)
    plt.scatter(xArr,yArr,s=0.1)
#    xy = ham.lonLatToXY(181.,0.)
#    plt.scatter([xy.x],[xy.y],s=100)

def makeLatexName(inputStr):
    tempStr = inputStr
    for c in ['_','=',' ']:
        if tempStr.find(c) == -1:
            return tempStr
        tmpStr = tempStr[:tempStr.find(c)]+'\\'+c
        tempStr = tempStr[tempStr.find(c)+1:]
        while tempStr.find(c) != -1:
            tmpStr += tempStr[:tempStr.find(c)]+'\\'+c
            tempStr = tempStr[tempStr.find(c)+1:]
        tmpStr += tempStr
        tempStr = tmpStr
    return tempStr

def plotEPA(imageName,imageNameEPA,epa):
    if os.path.isfile(imageName):
        print('imageName = ',imageName)
        im = Image.open(imageName)
        width, height = im.size
        d = ImageDraw.Draw(im)

        if width < height:
            r = width / 2.
        else:
            r = height / 2.
        print('r = ',r)
        center = [width/2,height/2]

        epaRad = math.radians(epa)

        topLeft = (center[0] - r*math.sin(epaRad), center[1] - r*math.cos(epaRad))
        bottomRight = (center[0] + r*math.sin(epaRad), center[1] + r*math.cos(epaRad))
        print('width = ',width,', height = ',height,': center = ',center,'; topLeft = ',topLeft,', bottomRight = ',bottomRight)

        line_color = (255, 255, 0)

    #    d.line([(center[0]-(width/10.),center[1]),(center[0]+(width/10.),center[1])],fill=line_color, width=2)
    #    d.line([(center[0],center[1]-(width/10.)),(center[0],center[1]+(width/10.))],fill=line_color, width=2)
        d.line([topLeft,bottomRight], fill=line_color, width=2)

        im.save(imageNameEPA)

if __name__ == "__main__":
    #print('data.getData("EPA") = ',data.getData('EPA'))
    hashData = None
    useHashFileName = None
    if mockSample or mockRandomSample:
        useHashFileName = hashFileName
        hashData = csvFree.readCSVFile(useHashFileName)
    else:
        newHashFileName = hashFileName[:-4]+'_withGPA.csv'
        print('dataFileName = <'+dataFileName+'>, hashFileName = <'+hashFileName+'>')
        findPNeWithPAinHASH(dataFileName, hashFileName, newHashFileName)
        hashData = csvFree.readCSVFile(newHashFileName)
        useHashFileName = newHashFileName

    iSec = 1

    with open(os.path.join(outPath,'numbers.txt'),'w') as numbers:
        for area in areas:
            xRange = area[0]
            yRange = area[1]
            tempOutPath = os.path.join(outPath,hashFileName[hashFileName.rfind('/')+1:-4])
            print('tempOutPath = <'+tempOutPath+'>')
            imOutPath = os.path.join(tempOutPath,'x=%.5f-%.5f_y=%.5f-%.5f/' % (xRange[0],xRange[1],yRange[0],yRange[1]))
            print('imOutPath = <'+imOutPath+'>')
            if not os.path.exists(imOutPath):
                os.makedirs(imOutPath)
            latexFileNameSummary = os.path.join(imOutPath,'latex_x=%.5f-%.5f_y=%.5f-%.5f.tex' % (xRange[0],xRange[1],yRange[0],yRange[1]))
            with open(latexFileNameSummary,'w') as texFileSummary:
                texFileSummary.write('\\documentclass{article}\n')
                texFileSummary.write('\\usepackage{graphicx}\n')
                texFileSummary.write('\\usepackage{xcolor}\n')
                texFileSummary.write('\\usepackage{float}\n')
                texFileSummary.write('\\begin{document}\n')
                for maxFlag in [3]:#np.arange(1,4,1):
                    for mainClass in mainClasses:
                        for rosePlotArea in [True,]:#False,True]:#,False]:
                            for rosePlotSmooth in [False,]:#True,False]:
                                for truePNeOnly in [True]:#,False]:
                                    if False:#(maxFlag == 2) and (mainClass == ['B']) and (not rosePlotSmooth):
                #                        print('show = ',show)
                #                        STOP
                                        show = True
                #                        print('turned on "show"')
                #                        STOP
                                    else:
                                        show = False
                                    numbers.write('%s %d: ' % (mainClass,maxFlag))
                                    fNameSuffix = '_flag_le_'+str(maxFlag)+'_'
                                    for c in mainClass:
                                        fNameSuffix += c
                                    fNameSuffix += '_true'
                                    if not truePNeOnly:
                                        fNameSuffix += '+likely'
                                    fNameSuffix += '_x=%.3f-%.3f_y=%.3f-%.3f' % (xRange[0], xRange[1], yRange[0], yRange[1])
                                    fNameSuffix += '_p=%.4f'
                                    if rosePlotArea:
                                        fNameSuffix += '_area'
                                    if rosePlotSmooth:
                                        fNameSuffix += '_smoothed'
                                    fNameSuffix += '_mean%.1f'
                                    fNameSuffix += '_dev%.1f'
                                    fNameSuffix += '.pdf'

                                    lbxyGPA = [] #Glon, Glat, x, y, GPA, flag, csGlon, csGlat, [dist,...]
            #                        hashIDs = []
                                    latexFileName = os.path.join(imOutPath, 'PN-alignments.tex')

                                    with open(latexFileName,'w') as texFile:
                                        texFile.write('\\documentclass{article}\n')
                                        texFile.write('\\usepackage{longtable}\n')
                                        texFile.write('\\usepackage{xcolor}\n')
                                        texFile.write('\\usepackage{graphicx}\n')
                                        texFile.write('\\usepackage[space]{grffile}')
                                        texFile.write('\\usepackage[legalpaper, landscape, margin=0.5in]{geometry}\n')
                                        texFile.write('\\begin{document}\n')
                                        texFile.write('\\begin{longtable}{| p{.03\\textwidth} | p{.05\\textwidth} | p{.05\\textwidth} | p{.05\\textwidth} | p{.05\\textwidth} | p{.03\\textwidth} | p{.03\\textwidth} | p{.03\\textwidth} | p{.02\\textwidth} | p{.06\\textwidth} | p{.04\\textwidth} | p{.03\\textwidth} | p{.06\\textwidth} | p{.30\\textwidth} |}\n')
                                        texFile.write('\\hline\n')
                                        texFile.write('HASH ID & DRAJ2000 & DDECJ2000 & l & b & EPA & QAP & diff & flag & GPA & class & majDiam & source & image\\\\ \n')
                                        texFile.write('\\hline\n')
                                        prev = -1
                                        prevPrev = -1
                                        for iLine in np.arange(0,hashData.size(),1):
                                            hashID = hashData.getData('idPNMain', iLine)
            #                                hashIDs.append(hashID)
            #                                hashLine = hashData.getData('idPNMain', iLine)
            #                                if hashLine == -1:
            #                                    STOP
                                            imageName = os.path.join(path, 'thumbnails/'+hashID+'.png')
                                #            print('imageName = <'+imageName+'>')
                                            if not os.path.isfile(imageName):
                                                print('file <'+imageName+'> not found')
                                #                files = glob.glob(os.path.join(path,'thumbnails/'+str(hashID)+'[a-z].png'))
                                #                print('files = ',files)
                                                if prevPrev == hashID:
                                                    imageName = os.path.join(path, 'thumbnails/'+hashID+'c.png')
                                                elif prev == hashID:
                                                    imageName = os.path.join(path, 'thumbnails/'+hashID+'b.png')
                                                else:
                                                    imageName = os.path.join(path, 'thumbnails/'+hashID+'a.png')
                                                if not os.path.isfile(imageName):
                                                    print('file <'+imageName+'> not found')
                                            imageNameEPA = imageName[:imageName.rfind('.')]+'_EPA.png'
                                            plotEPA(imageName,imageNameEPA,float(hashData.getData('EPA', iLine)))
                                            texFile.write(hashID+' & ')
                                            texFile.write(hashData.getData('DRAJ2000', iLine)+' & ')
                                            texFile.write(hashData.getData('DDECJ2000', iLine)+' & ')
                                            texFile.write(hashData.getData('Glon', iLine)+' & ')
                                            texFile.write(hashData.getData('Glat', iLine) + ' & ')
                                            texFile.write(hashData.getData('EPA', iLine)+' & ')
                                            qapTable = csv.DictReader(open('/Users/azuri/daten/uni/HKU/PN_orientation_angles/HASH-True-PN-EPA-measures.csv'))
                                            color = 'black'
                                            foundQAP = False
                                            for row in qapTable:
                                                if (row['HASH ID'] == hashID) and (row['EPA Ziggy'] == hashData.getData('EPA', iLine)):
                                                    texFile.write(row['QAP EPA']+' & ')
                                                    difference = row['ZIG-QAP(EPA)']
                                                    if (difference == '#VALUE!') or (difference == '') or (difference == ' '):
                                                        difference = '0'
                                                    difference = float(difference)
                                                    if (abs(difference) > 20) and (abs(180.-abs(difference)) > 20.):
                                                        color='red'
                                                    texFile.write(row['ZIG-QAP(EPA)'].strip('#')+' & ')
                                                    foundQAP = True
                                                    break
                                            if not foundQAP:
                                                texFile.write('N/A & ')
                                                texFile.write('N/A & ')

                                            texFile.write('\\textcolor{'+color+'}{'+hashData.getData('flag', iLine)+'} & ')
                                            texFile.write('\\textcolor{'+color+'}{'+hashData.getData('GPA', iLine)+'} & ')
                                            texFile.write('\\textcolor{'+color+'}{'+hashData.getData('mainClass', iLine) + hashData.getData('subClass', iLine) + '} & ')
                                            texFile.write('\\textcolor{'+color+'}{'+hashData.getData('MajDiam', iLine) + '} & ')
                                            texFile.write('\\textcolor{'+color+'}{'+hashData.getData('source', iLine)+'} & ')
                                            if os.path.isfile(imageName):
                                                texFile.write('\\centerline{\\includegraphics[height=0.3\\textheight]{'+imageNameEPA+'}}\n')
                                            texFile.write('\\\\ \\hline\n')
                                            prevPrev = prev
                                            prev = hashID

                                            if (truePNeOnly and (hashData.getData('PNstat',iLine) == 'T')) or (not truePNeOnly):
                                                if hashData.getData('mainClass', iLine) in mainClass:
                                                    l = float(hashData.getData('Glon', iLine))
                                                    b = float(hashData.getData('Glat', iLine))
                                                    xy = ham.lonLatToXY(l,b)
                #                                    print('xy = ',xy)
                                                    x = xy.x
                                                    y = xy.y
                #                                    print('x = ',x,', y = ',y)

                                                    gpa = float(hashData.getData('GPA', iLine))
                                                    if rotateBipolars and (hashData.getData('mainClass', iLine) == 'B'):
                                                        gpa += 90.
                                                    if gpa < 0.:
                                                        gpa += 180.
                                                    if gpa > 180.:
                                                        gpa -= 180.

                                                    csGlon = hashData.getData('CS_Glon', iLine)
                                                    csGlat = hashData.getData('CS_Glat', iLine)

                                                    dist = None
                #                                    if (csGlon != '') and (csGlat != ''):
                                        #                getStarWithMinDist(gaiaData, ra, dec, iStar=0)
                                        #                c = SkyCoord(ra=float(csRa)*u.degree, dec=float(csDec)*u.degree, frame='icrs')
                                        #                lb = c.galactic
                                        #                print('dir(lb) = ',dir(lb))
                #                                        print('l = ',l,', b = ',b)
                                                    if (x >= xRange[0]) and (x <= xRange[1]) and (y >= yRange[0]) and (y <= yRange[1]):
                                                        #print('hashData.getData(',iLine,') = ',hashData.getData(iLine))
                                                        lbxyGPA.append({'l':l,# 0
                                                                        'b':b,# 1
                                                                        'x':x,# 2
                                                                        'y':y,# 3
                                                                        'gpa':gpa,# 4
                                                                        'flag':int(hashData.getData('flag', iLine)),# 5
                                                                        'csGlon':csGlon,# 6
                                                                        'csGlat':csGlat,# 7
                                                                        'dist':dist,# 8
                                                                        'mainClass':hashData.getData('mainClass', iLine),# 9
                                                                        'ID':hashID,# 10
                                                        })
                                                    if hashID == '113':
                                                        print(lbxyGPA)
                #                                            STOP

                                        texFile.write('\\caption{Your caption here} % needs to go inside longtable environment\n')
                                        texFile.write('\\label{tab:myfirstlongtable}\n')
                                        texFile.write('\\end{longtable}\n')
                                        texFile.write('%\\end{table}\n')

                                        texFile.write('Table \\ref{tab:myfirstlongtable} shows my first longtable.\n')
                                        texFile.write('\\end{document}\n')

                #                        print('len(lbxyGPA) = ',len(lbxyGPA))
                #                        print('lbxyGPA[:][0] = ',lbxyGPA[:][0])
        #                                print('lbxyGPA[0] = ',lbxyGPA[0])
        #                                print('lbxyGPA[0]["l"] = ',lbxyGPA[0]['l'])
                                        l = np.array([i['l'] for i in lbxyGPA])
                                        b = np.array([i['b'] for i in lbxyGPA])
                                        x = np.array([i['x'] for i in lbxyGPA])
                                        y = np.array([i['y'] for i in lbxyGPA])
                                        GPA = np.array([i['gpa'] for i in lbxyGPA])
        #                                print('lbxyGPA[0] = ',lbxyGPA[0])
                                        flag = np.array([i['flag'] for i in lbxyGPA])
                                        print('flag = ',flag)
                                        mClass = np.array([i['mainClass'] for i in lbxyGPA])
                                        oneFlags = np.array([f <= maxFlag for f in flag])
                                        print('oneFlags = ',oneFlags)
        #                                STOP
                                        ids = np.array([i['ID'] for i in lbxyGPA])

                                        # calculate moments
                                        moments = calcMoments(GPA[oneFlags])
                                        if moments[0][0] < 0.:
                                            moments[0][0] += 180.
                                        print('moments = ',moments)
                                        mean = calcMean(GPA[oneFlags])
                                        print('mean = ',mean)
                                        dev = (GPA[oneFlags].shape[0]-mean[1]) / 2.
                                        print('oneFlags.shape[0] = ',oneFlags.shape[0])
                                        print('GPA[oneFlags].shape[0] = ',GPA[oneFlags].shape[0])

                                        # calculate Rayleigh probability that the distribution is uniform
                                        gpaRad = np.array([np.radians(2.*g) for g in GPA[oneFlags]])
                                        pRad = circstats.rayleightest(gpaRad)
                                        print('pRad = ',pRad)
                                        print('length of R = ',mean[1])
                                        print('dev = ',dev)
                                        print('fNameSuffix = <'+fNameSuffix+'>')

                                        pltArr = []
                                        pltMClass = []
                                        nEllipticals = 0
                                        nBipolars = 0
                                        for i in np.arange(0,len(oneFlags),1):
                                            if oneFlags[i]:
                                                pltArr.append(GPA[i])
                                                pltArr.append(GPA[i]+180.)
                                                if mClass[i] == 'E':
                                                    pltMClass.append(True)
                                                    pltMClass.append(True)
                                                    nEllipticals += 1
                                                else:
                                                    pltMClass.append(False)
                                                    pltMClass.append(False)
                                                    nBipolars += 1
                                        print('found ',nEllipticals,' ellipticals and ',nBipolars,' bipolars')
                                        pltArr = np.array(pltArr)
                                        fNameSuffix = fNameSuffix % (pRad, mean[0], dev)

                                        texFileSummary.write('\\section{Main morphological class = ')
                                        for cl in mainClass:
                                            texFileSummary.write(cl)
                                        texFileSummary.write(', maximum flag = '+str(maxFlag)+', length of rose plot fingers proportional to '+('area' if rosePlotArea else 'length')+', rose plot is '+('smoothed' if rosePlotSmooth else 'not smoothed}\n'))
                                        tempHashFileName = useHashFileName[useHashFileName.rfind('/')+1:]
                                        texFileSummary.write('inputFile = '+makeLatexName(tempHashFileName)+'\\\\\n')
                                        texFileSummary.write(str(nEllipticals)+' elliptical PNe, '+str(nBipolars)+' bipolar PNe\\\\\n')
                                        texFileSummary.write('mean = '+str(mean)+', deviation = '+str(dev)+'\\\\~\\\\\n')
                                        texFileSummary.write('\\begin{tabular}{l | l}\n\\hline\n')
                                        texFileSummary.write('Test & results\\\\\n\\hline\n')
                                        statsFile = os.path.join(imOutPath,fNameSuffix[1:-3]+'dat')
                                        print('statsFile = <'+statsFile+'>')
                                        with open(statsFile,'w') as sFile:
                                            sFile.write('found %d ellipticals and %d bipolars\n' % (nEllipticals,nBipolars))
                                            sFile.write('pRayleigh = %.10f\n' % pRad)
                                            texFileSummary.write('Rayleigh & '+'%.10f\\\\\n' % pRad)
                                            pPyCircRad = pycircstat.rayleigh(gpaRad)
                                            print('pPyCircRad = ',pPyCircRad)

                                            pOmnibus = pycircstat.omnibus(gpaRad)
                                            print('Omnibus = ',pOmnibus)
                                            print('Omnibus[0] = ',pOmnibus[0])
                                            print('Omnibus[1] = ',pOmnibus[1])
                                            sFile.write('Omnibus = (%.10f, %.10f)\n' % (pOmnibus[0], pOmnibus[1][0]))
                                            texFileSummary.write('Omnibus & '+'%.10f, %.10f\\\\\n' % (pOmnibus[0], pOmnibus[1][0]))

                                            try:
                                                pRaoSpacing = pycircstat.raospacing(gpaRad)
                                                print('pRaoSpacing = ',pRaoSpacing)
                                                print('pRaoSpacing[0] = ',pRaoSpacing[0])
                                                print('pRaoSpacing[1] = ',pRaoSpacing[1])
                                                print('pRaoSpacing[2] = ',pRaoSpacing[2])
                                                sFile.write('pRaoSpacing = (%.10f, %.10f, %.10f)\n' % (pRaoSpacing[0], pRaoSpacing[1], pRaoSpacing[2]))
                                                texFileSummary.write('pRaoSpacing & '+'%.10f, %.10f, %.10f\\\\\n' % (pRaoSpacing[0], pRaoSpacing[1], pRaoSpacing[2]))
                                            except:
                                                sFile.write('pRaoSpacingTest failed\n')
                                                texFileSummary.write('pRaoSpacing & failed\\\\\n')
                                                pass

                                            pVTest = pycircstat.vtest(gpaRad, np.radians(mean))
                                            print('pVTest = ',pVTest)
                                            sFile.write('pVTest = [%.10f, %.10f], [%.10f, %.10f]\n' % (pVTest[0][0], pVTest[0][1], pVTest[1][0], pVTest[1][1]))
                                            texFileSummary.write('pVTest & '+'[%.10f, %.10f], [%.10f, %.10f]\\\\\n' % (pVTest[0][0], pVTest[0][1], pVTest[1][0], pVTest[1][1]))

                                            pSymTest = pycircstat.symtest(gpaRad, axis=None)
                                            print('pSymTest = ',pSymTest)
                                            sFile.write('pSymTest = (%.10f, %.10f)\n' % (pSymTest[0], pSymTest[1]))
                                            texFileSummary.write('pSymTest & '+'%.10f, %.10f\\\\\n' % (pSymTest[0], pSymTest[1]))

                                            try:
                                                pMTest = pycircstat.mtest(gpaRad, np.radians(mean))
                                                print('pMTest = ',pMTest)
                                                sFile.write('pMTest = [%r, %r], %.10f, confidence interval(lower=%.10f, upper=%.10f)\n' % (pMTest[0][0], pMTest[0][1], pMTest[1], pMTest[2].lower, pMTest[2].upper))
                                                texFileSummary.write('pMTest & '+'[%r, %r], %.10f, confidence interval(lower=%.10f, upper=%.10f)\\\\\n' % (pMTest[0][0], pMTest[0][1], pMTest[1], pMTest[2].lower, pMTest[2].upper))

                                            except:
                                                sFile.write('pMTest failed\n')
                                                texFileSummary.write('pMTest & failed\\\\\n')
                                                pass
                                            p = pRad
                                            texFileSummary.write('\\hline\\end{tabular}\n')

                                        uni = np.unique(ids)
                                        print('unique ids = ',uni)
                                        uni = np.unique(ids[oneFlags])
                                        print('unique ids with flag = ',uni)
                                        outf = os.path.join(imOutPath,'uniqueIDs_'+fNameSuffix[1:-3]+'dat')
                                        for c in mainClass:
                                            outf += c
                                        outf += str(maxFlag)+'.csv'
                                        with open(outf,'w') as tf:
                                            for u in uni:
                                                tf.write('%d\n' % int(u))
                                        numbers.write('%d EPAs from %d PNe\n' % (GPA[oneFlags].shape[0],uni.shape[0]))


                                        # plot number of PNe per quadrant
                                        nStarsPerQuadrant = np.zeros(4)
                                        for i in np.arange(0,l[oneFlags].shape[0],1):
                                            if x[oneFlags][i] > 0:
                                                if y[oneFlags][i] > 0:
                                                    nStarsPerQuadrant[0] += 1
                                                else:
                                                    nStarsPerQuadrant[3] += 1
                                            else:
                                                if y[oneFlags][i] > 0:
                                                    nStarsPerQuadrant[1] += 1
                                                else:
                                                    nStarsPerQuadrant[2] += 1
                                        ang = 45.
                                        qAngles = []
                                        for i in np.arange(0,4,1):
                                            for j in np.arange(0,nStarsPerQuadrant[i],1):
                                                qAngles.append(ang)
                                            ang += 90.
                                        fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
                                        rose_plot(ax,
                                                  qAngles,
                                                  offset=0,
                                                  bins=4,
                                                  start_zero=True,
                                                  color='blue',
                                                  fill=True,
                                                  density=rosePlotArea,
                                                 )
                                        fig.tight_layout()
                                        fig.savefig(os.path.join(imOutPath,'rose_plot_nStars_per_quadrant'+fNameSuffix), bbox_inches='tight')
                                        if show:
                                            plt.show()
                                        else:
                                            plt.close(fig)

                                        # plot Hammer projection
                                        eFlags = []
                                        bFlags = []
                                        for i in np.arange(0,len(mClass[oneFlags]),1):
                                            if mClass[oneFlags][i] == 'B':
                                                eFlags.append(False)
                                                bFlags.append(True)
                                            else:
                                                bFlags.append(False)
                                                eFlags.append(True)

                                        plotHammerProjection(x,y,GPA,oneFlags,eFlags,bFlags,os.path.join(imOutPath,'all_pne'+fNameSuffix))
                                        texFileSummary.write('\\begin{figure}[H]\n')
                                        texFileSummary.write('\\centering\n')
                                        tempName = 'all_pne'+fNameSuffix
                                        texFileSummary.write('\\includegraphics[width=\\textwidth]{{{'+tempName[:tempName.rfind('.')]+'}}}\n')
                                        texFileSummary.write('\\caption{All measured PNe in our sample. Circles are elliptical PNe, diamonds are biploar PNe. The colour code is their Galactic Position Angle (GPA).\n')
                                        texFileSummary.write('}\n')
                                        texFileSummary.write('\\label{fig1}\n')
                                        texFileSummary.write('\\end{figure}\n')

                                        fig = plt.figure(figsize=(5,4.7))
                                        plt.hist(GPA[oneFlags], bins=int(nBins/2))
                                        plt.xlabel('GPA [degrees]',fontsize=14)
                                        plt.ylabel('Number of PNe',fontsize=14)
                                        plt.xticks(fontsize=14)
                                        plt.yticks(fontsize=14)
                                        fig.tight_layout()
                                        fig.savefig(os.path.join(imOutPath,'histogram'+fNameSuffix), bbox_inches='tight')
                                        if show:
                                            plt.show()
                                        else:
                                            plt.close(fig)
                                        print('after histogram: ',plt.get_fignums(),' figures open')

                #                        print('np.max(GPA[oneFlags]) = ',np.max(GPA[oneFlags]))

                                        fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
                                        max_count, max_size = rose_plot(ax,
                                                                        pltArr,
                                                                        offset=np.pi/2.,
                                                                        bins=nBins,
                                                                        start_zero=True,
                                                                        color='blue',
                                                                        fill=True,
                                                                        density=rosePlotArea,
                                                                        smooth=rosePlotSmooth,
                                                                        mean=mean[0],
                                                                        #sigma=dev,
                                                                        )
                                        max_count, max_size = rose_plot(ax,
                                                                        pltArr[pltMClass],
                                                                        offset=np.pi/2.,
                                                                        bins=nBins,
                                                                        start_zero=True,
                                                                        color='green',
                                                                        fill=True,
                                                                        max_count=max_count,
                                                                        max_size=max_size,
                                                                        density=rosePlotArea,
                                                                        smooth=rosePlotSmooth,
                                                                        mean=mean[0],
                                                                        #sigma=dev,
                                                                        )
                #                        if rosePlotArea:
                #                            plt.plot([-np.pi,np.pi],[0,max_size],co)
                #                        else:
                #                            plt.plot([-np.pi,np.pi],[0,max_count])
                                        fig.tight_layout()
                                        fig.savefig(os.path.join(imOutPath,'rose_plot'+fNameSuffix), bbox_inches='tight')
                                        if show:
                                            plt.show()
                                        else:
                                            plt.close(fig)
                                        print('after rose_plot: ',plt.get_fignums(),' figures open')

                                        alpha, radius = vectorDiagram(GPA[oneFlags]*2.,os.path.join(imOutPath,'vectorDiagram'+fNameSuffix))
                                        print('alpha = ',alpha/2.,', radius = ',radius)
                                        print('after vectorDiagram: ',plt.get_fignums(),' figures open')
                                        alpha, radius = vectorDiagram(GPA[oneFlags],os.path.join(imOutPath,'vectorDiagramHalf'+fNameSuffix))
                #                        alpha, radius = vectorDiagram(np.arange(0,90,1))
                #                        alpha, radius = vectorDiagram(np.arange(0,180,1))
                #                        alpha, radius = vectorDiagram(np.arange(0,360,1))
                                        linearOrderDiagram(GPA[oneFlags]*2., os.path.join(imOutPath,'linearOrderDiagram'+fNameSuffix))
                                        print('after linearOrderDiagram: ',plt.get_fignums(),' figures open')
                #                        if show:
                #                            STOP

                #                        xyBulge = ham.lonLatToXY(-15., -5.)
                #                        print('xyBulge = ',xyBulge.x,xyBulge.y)
                                        texFileSummary.write('\\begin{figure}[H]\n')
                                        texFileSummary.write('\\centering\n')
                                        tempName = 'rose_plot_nStars_per_quadrant'+fNameSuffix
                                        texFileSummary.write('\\hfill\\includegraphics[width=0.48\\textwidth]{{{'+tempName[:tempName.rfind('.')]+'}}}\\hfill\n')
                                        tempName = 'histogram'+fNameSuffix
                                        texFileSummary.write('\\includegraphics[width=0.48\\textwidth]{{{'+tempName[:tempName.rfind('.')]+'}}}\\hfill~\\\\\n')
                                        tempName = 'linearOrderDiagram'+fNameSuffix
                                        texFileSummary.write('\\hfill\\includegraphics[width=0.48\\textwidth]{{{'+tempName[:tempName.rfind('.')]+'}}}\\hfill\n')
                                        tempName = 'rose_plot'+fNameSuffix
                                        texFileSummary.write('\\includegraphics[width=0.48\\textwidth]{{{'+tempName[:tempName.rfind('.')]+'}}}\\hfill~\\\\\n')
                                        tempName = 'vectorDiagram'+fNameSuffix
                                        texFileSummary.write('\\hfill\includegraphics[width=0.48\\textwidth]{{{'+tempName[:tempName.rfind('.')]+'}}}\\hfill\n')
                                        tempName = 'vectorDiagramHalf'+fNameSuffix
                                        texFileSummary.write('\\includegraphics[width=0.48\\textwidth]{{{'+tempName[:tempName.rfind('.')]+'}}}\\hfill\\ \n')
                                        texFileSummary.write('\\caption{Rose plot for the number of PNe per quadrant (top left), Histogram (top right), Uniformity plot (center left), '+('smoothed ' if rosePlotSmooth else '')+('area ' if rosePlotArea else 'number ')+'Rose plot (center right), Vector plot for doubled GPAs [0,360) (bottom left), Vector plot for actual GPAs [0,180) (bottom right) of GPAs for '+str(nEllipticals)+' elliptical and '+str(nBipolars)+' bipolar PNe with the confidence flag $\\leq'+str(maxFlag)+'$. Deviations from the diagonal in the Uniformity plot indicate non-uniformity. Only the left side of the Rose plot has been measured and then rotated by 180 degrees for better visibility. The mean GPA has been marked with a red line. \n')
                                        texFileSummary.write('}\n')
                                        texFileSummary.write('\\label{fig2}\n')
                                        texFileSummary.write('\\end{figure}\n')


                        iSec += 1
                texFileSummary.write('\\end{document}\n')
#    subprocess.run(["cd",imOutPath])
#    subprocess.run(["pdflatex", latexFileNameSummary[latexFileNameSummary.rfind('/')+1:]])
#    pdfFileName = latexFileNameSummary[latexFileNameSummary.rfind('/')+1:]
#    subprocess.run(["open",pdfFileName])
#    subprocess.run(["cd","/Users/azuri/entwicklung/python"])
