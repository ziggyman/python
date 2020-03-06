import glob
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import math

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import circstats

import csvData
import csvFree
from galaxyMath import raDecToLB, parallaxToDistance,degToRad,radToDeg
from hammer import Pixel,XY,LonLat,Hammer
from myUtils import getStarWithMinDist

path = '/Users/azuri/daten/uni/HKU/PN alignment'
dataFileName = os.path.join(path, 'PN-alignments.csv')
hashFileName = os.path.join(path, 'HASH_bipolar+elliptical_true_PNe.csv')
latexFileName = os.path.join(path, 'PN-alignments.tex')
gaiaFileNameRoot = '/Volumes/work/azuri/data/gaia/dr2/xy/GaiaSource_%.6f-%.6f_%.6f-%.6f.csv'
imOutPath = '/Users/azuri/entwicklung/tex/Poster/workplans2/images/'

data = csvFree.readCSVFile(dataFileName)
hashData = csvFree.readCSVFile(hashFileName)

show = False
# Turn interactive plotting off
#plt.ioff()

ham = Hammer()
pixels = ham.getPixels()
xRangeMax = [-2.83,2.83]
yRangeMax = [-1.4143,1.4143]
xRange = xRangeMax#[-0.236,0.261]
yRange = yRangeMax#[-0.0873,0.0873]

def findInHash(hashData, hashID):
    ids = hashData.getData('idPNMain')
    for iLine in np.arange(0,hashData.size(),1):
        if ids[iLine] == hashID:
            return iLine

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

# FROM https://stackoverflow.com/questions/22562364/circular-histogram-for-python
#@param: angles: angles over full circle in degrees
#@param: mean: mean of angles over full circle in degrees
#@param: sigma: deviation of angles over full circle in degrees
def rose_plot(ax, angles, bins=16, density=None, offset=0, lab_unit="degrees",
              start_zero=False, fill=False, color='white', max_count=None,
              max_size=None, smooth=None, mean=None, sigma=None, **param_dict):
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
    if smooth is not None:
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
        ax.bar([math.radians(mean),math.radians(mean+180)],[maxR,maxR],width=0.01,linewidth=5,color='black')
        if sigma is not None:
            maxR = np.max(radius)*1.05
            ax.bar([math.radians(mean+sigma),math.radians(mean+sigma+180)],[maxR,maxR],width=0.01,linewidth=5,color='red')
            ax.bar([math.radians(mean-sigma),math.radians(mean-sigma+180)],[maxR,maxR],width=0.01,linewidth=5,color='red')

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
def vectorDiagram(angles, fNameOut):
    fig = plt.figure(figsize=(5,4.7))
    plt.axis('off')
    a = np.zeros(2)
    radius = []
    for angle in angles:
        radian = math.radians(angle)
        x = [a[0],a[0]-np.sin(radian)]
        y = [a[1],a[1]+np.cos(radian)]
        plt.plot(x,y)
        a[0] = a[0]-np.sin(radian)
        a[1] = a[1]+np.cos(radian)
        radius.append(np.sqrt(a[0]**2 + a[1]**2))
    maxRadius = np.amax(radius)
    plt.plot([0,a[0]],[0,a[1]],'r-',linewidth=3)
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
    fig.savefig(fNameOut, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)
    radius = np.sqrt(a[0]**2 + a[1]**2)
    alpha = np.degrees(np.arcsin(-a[0]/radius))
    beta = np.degrees(np.arccos(a[1]/radius))
    gamma = np.degrees(np.arctan2(a[1],a[0]))
    print('vectorDiagram: alpha = ',alpha,', beta = ',beta,', gamma = ',gamma)
    return [gamma,radius]

#@param angles: angles over full circle in degrees
def linearOrderDiagram(angles, fNameOut):
    rad = np.array([math.radians(angle) for angle in angles])
    linearOrderX = np.sort(rad/(2.*np.pi))
    linearOrderY = np.zeros(linearOrderX.shape)
    for iX in np.arange(0,linearOrderX.shape[0],1):
        linearOrderY[iX] = iX / (linearOrderX.shape[0]+1)

    fig = plt.figure(figsize=(5,5))
    plt.scatter(linearOrderY,linearOrderX)
    plt.plot([0,1],[0,1])
    fig.tight_layout()
    fig.savefig(fNameOut, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)

# @brief return input array lbxyGPA with values inside [x0,x1], [y0,y1]
def selectXY(lbxyGPA_in, x0, x1, y0, y1):
    l = np.array([i[0] for i in lbxyGPA_in])
    b = np.array([i[1] for i in lbxyGPA_in])
    x = np.array([i[2] for i in lbxyGPA_in])
    y = np.array([i[3] for i in lbxyGPA_in])
    GPA = np.array([i[4] for i in lbxyGPA_in])
    flag = np.array([i[5] for i in lbxyGPA_in])

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

if __name__ == "__main__":
    #print('data.getData("EPA") = ',data.getData('EPA'))
    for maxFlag in np.arange(1,4,1):
        for mainClass in [['E'],['B'],['B','E']]:
            for rosePlotArea in [True,False]:
                for rosePlotSmooth in [True,False]:
                    if (maxFlag == 2) and (mainClass == ['B']) and (not rosePlotSmooth):
#                        print('show = ',show)
#                        STOP
                        show = True
#                        print('turned on "show"')
#                        STOP
                    else:
                        show = False
                    fNameSuffix = '_flag_le_'+str(maxFlag)+'_'
                    for c in mainClass:
                        fNameSuffix += c
                    fNameSuffix += '_x=%.3f-%.3f_y=%.3f-%.3f' % (xRange[0], xRange[1], yRange[0], yRange[1])
                    fNameSuffix += '_p=%.4f'
                    if rosePlotArea:
                        fNameSuffix += '_area'
                    if rosePlotSmooth:
                        fNameSuffix += '_smoothed'
                    fNameSuffix += '.eps'

                    lbxyGPA = [] #Glon, Glat, x, y, GPA, flag, csGlon, csGlat, [dist,...]
                    hashIDs = []

                    with open(latexFileName,'w') as texFile:
                        texFile.write('\\documentclass{article}\n')
                        texFile.write('\\usepackage{longtable}\n')
                        texFile.write('\\usepackage{graphicx}\n')
                        texFile.write('\\usepackage[space]{grffile}')
                        texFile.write('\\usepackage[legalpaper, landscape, margin=0.5in]{geometry}\n')
                        texFile.write('\\begin{document}\n')
                        texFile.write('\\begin{longtable}{| p{.03\\textwidth} | p{.05\\textwidth} | p{.05\\textwidth} | p{.05\\textwidth} | p{.05\\textwidth} | p{.03\\textwidth} | p{.02\\textwidth} | p{.06\\textwidth} | p{.04\\textwidth} | p{.03\\textwidth} | p{.12\\textwidth} | p{.30\\textwidth} |}\n')
                        texFile.write('\\hline\n')
                        texFile.write('HASH ID & DRAJ2000 & DDECJ2000 & l & b & EPA & flag & GPA & class & majDiam & source & image\\\\ \n')
                        texFile.write('\\hline\n')
                        prev = -1
                        prevPrev = -1
                        for iLine in np.arange(0,data.size(),1):
                            if data.getData('EPA',iLine) != '':
                                hashID = data.getData('HASH ID', iLine)
                                hashIDs.append(hashID)
                                hashLine = findInHash(hashData, hashID)
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
                                texFile.write(hashID+' & ')
                                texFile.write(data.getData('DRAJ2000', iLine)+' & ')
                                texFile.write(data.getData('DDECJ2000', iLine)+' & ')
                                texFile.write(data.getData('l', iLine)+' & ')
                                texFile.write(data.getData('b', iLine) + ' & ')
                                texFile.write(data.getData('EPA', iLine)+' & ')
                                texFile.write(data.getData('flag', iLine)+' & ')
                                texFile.write(data.getData('GPA', iLine)+' & ')
                                texFile.write(hashData.getData('mainClass', hashLine) + hashData.getData('subClass', hashLine) + ' & ')
                                texFile.write(hashData.getData('MajDiam', hashLine) + ' & ')
                                texFile.write(data.getData('source', iLine)+' & ')
                                if os.path.isfile(imageName):
                                    texFile.write('\\centerline{\\includegraphics[height=0.3\\textheight]{'+imageName+'}}\n')
                                texFile.write('\\\\ \\hline\n')
                                prevPrev = prev
                                prev = hashID

                                if hashData.getData('mainClass', iLine) in mainClass:
                                    l = float(hashData.getData('Glon', hashLine))
                                    b = float(hashData.getData('Glat', hashLine))
                                    xy = ham.lonLatToXY(l,b)
#                                    print('xy = ',xy)
                                    x = xy.x
                                    y = xy.y
#                                    print('x = ',x,', y = ',y)

                                    gpa = float(data.getData('GPA', iLine))
                                    if gpa < 0.:
                                        gpa += 180
                                    if gpa > 180:
                                        gpa -= 180.

                                    csGlon = hashData.getData('CS_Glon', hashLine)
                                    csGlat = hashData.getData('CS_Glat', hashLine)

                                    dist = None
#                                    if (csGlon != '') and (csGlat != ''):
                        #                getStarWithMinDist(gaiaData, ra, dec, iStar=0)
                        #                c = SkyCoord(ra=float(csRa)*u.degree, dec=float(csDec)*u.degree, frame='icrs')
                        #                lb = c.galactic
                        #                print('dir(lb) = ',dir(lb))
#                                        print('l = ',l,', b = ',b)
                                    if (x >= xRange[0]) and (x <= xRange[1]) and (y >= yRange[0]) and (y <= yRange[1]):
                                        lbxyGPA.append([l, b, x, y, gpa, int(data.getData('flag', iLine)), csGlon, csGlat, dist, hashData.getData('mainClass', iLine)])

                        texFile.write('\\caption{Your caption here} % needs to go inside longtable environment\n')
                        texFile.write('\\label{tab:myfirstlongtable}\n')
                        texFile.write('\\end{longtable}\n')
                        texFile.write('%\\end{table}\n')

                        texFile.write('Table \\ref{tab:myfirstlongtable} shows my first longtable.\n')
                        texFile.write('\\end{document}\n')

#                        print('len(lbxyGPA) = ',len(lbxyGPA))
#                        print('lbxyGPA[:][0] = ',lbxyGPA[:][0])
                        l = np.array([i[0] for i in lbxyGPA])
                        b = np.array([i[1] for i in lbxyGPA])
                        x = np.array([i[2] for i in lbxyGPA])
                        y = np.array([i[3] for i in lbxyGPA])
                        GPA = np.array([i[4] for i in lbxyGPA])
                        flag = np.array([i[5] for i in lbxyGPA])
                        mClass = np.array([i[9] for i in lbxyGPA])
                        oneFlags = np.array([f <= maxFlag for f in flag])
#                        print('oneFlags = ',type(oneFlags),': ',oneFlags)

                        p = circstats.rayleightest(GPA[oneFlags] * 2.)
                        print('p(GPA*2) = ',p)

                        # plot Hammer projection
                        fig = plt.figure(figsize=(25,10))
                        plt.scatter(x[oneFlags],y[oneFlags],c=GPA[oneFlags],s=20,cmap='viridis')
                        plotLBMarks(10)
                        cbarTicks = np.arange(0.,1.0001,20./180.)
                        cbar = plt.colorbar(ticks=cbarTicks)
                        cbar.set_label('Galactic Position Angle')
                        cbarTicks = cbarTicks*180.
                        cbarTicks = np.round(cbarTicks)
                        cbarTicks = [int(x) for x in cbarTicks]
                        cbar.ax.set_yticklabels(cbarTicks)
                        cbar.ax.tick_params(labelsize=18)
                        ax = cbar.ax
                        text = ax.yaxis.label
                        font = matplotlib.font_manager.FontProperties(size=18)#family='times new roman', style='italic',
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
                        fig.savefig(os.path.join(imOutPath,'all_pne'+fNameSuffix % (p)), bbox_inches='tight')
                        if show:
                            plt.show()
                        else:
                            plt.close(fig)
                        print('after all_pne: ',plt.get_fignums(),' figures open')

                        fig = plt.figure()
                        plt.hist(GPA[oneFlags], bins=18)
                        plt.xlabel('GPA [degrees]')
                        plt.ylabel('Number of PNe')
                        fig.savefig(os.path.join(imOutPath,'histogram'+fNameSuffix % (p)), bbox_inches='tight')
                        if show:
                            plt.show()
                        else:
                            plt.close(fig)
                        print('after histogram: ',plt.get_fignums(),' figures open')

                        moments = calcMoments(GPA[oneFlags])
                        print('moments = ',moments)
                        mean = calcMean(GPA[oneFlags])
                        print('mean = ',mean)

#                        print('np.max(GPA[oneFlags]) = ',np.max(GPA[oneFlags]))

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

                        fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
                        max_count, max_size = rose_plot(ax,
                                                        pltArr,
                                                        offset=np.pi/2.,
                                                        bins=36,
                                                        start_zero=True,
                                                        color='blue',
                                                        fill=True,
                                                        density=rosePlotArea,
                                                        smooth=rosePlotSmooth,
                                                        mean=moments[0][0]*2.,
                                                        sigma=moments[1][0]*2.)
                        max_count, max_size = rose_plot(ax,
                                                        pltArr[pltMClass],
                                                        offset=np.pi/2.,
                                                        bins=36,
                                                        start_zero=True,
                                                        color='green',
                                                        fill=True,
                                                        max_count=max_count,
                                                        max_size=max_size,
                                                        density=rosePlotArea,
                                                        smooth=rosePlotSmooth,
                                                        mean=moments[0][0]*2.,
                                                        sigma=moments[1][0]*2.)
#                        if rosePlotArea:
#                            plt.plot([-np.pi,np.pi],[0,max_size],co)
#                        else:
#                            plt.plot([-np.pi,np.pi],[0,max_count])
                        fig.tight_layout()
                        fig.savefig(os.path.join(imOutPath,'rose_plot'+fNameSuffix % (p)), bbox_inches='tight')
                        if show:
                            plt.show()
                        else:
                            plt.close(fig)
                        print('after rose_plot: ',plt.get_fignums(),' figures open')

                        alpha, radius = vectorDiagram(GPA*2.,os.path.join(imOutPath,'vectorDiagram'+fNameSuffix % (p)))
                        print('alpha = ',alpha/2.,', radius = ',radius)
                        print('after vectorDiagram: ',plt.get_fignums(),' figures open')
                        linearOrderDiagram(GPA*2., os.path.join(imOutPath,'linearOrderDiagram'+fNameSuffix % (p)))
                        print('after linearOrderDiagram: ',plt.get_fignums(),' figures open')
                        if show:
                            STOP

#                        xyBulge = ham.lonLatToXY(-15., -5.)
#                        print('xyBulge = ',xyBulge.x,xyBulge.y)
