import glob
import numpy as np
import os
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

ham = Hammer()
pixels = ham.getPixels()
lbxyGPA = [] #Glon, Glat, x, y, GPA, flag, csGlon, csGlat, [dist,...]
hashIDs = []
maxFlag = 2
xRangeMax = [-2.83,2.83]
yRangeMax = [-1.4143,1.4143]
xRange = xRangeMax#[-0.236,0.261]
yRange = yRangeMax#[-0.0873,0.0873]
fNameSuffix = '_flag_le_'+str(maxFlag)+'_'
mainClass = ['E','B']#
for c in mainClass:
    fNameSuffix += c
fNameSuffix += '_x=%.3f-%.3f_y=%.3f-%.3f' % (xRange[0], xRange[1], yRange[0], yRange[1])
fNameSuffix += '_p=%.4f'
fNameSuffix += '.eps'

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
    return np.array([[math.degrees(moment[0]), moment[1]] for moment in moments])

# FROM https://stackoverflow.com/questions/22562364/circular-histogram-for-python
def rose_plot(ax, angles, bins=16, density=None, offset=0, lab_unit="degrees",
              start_zero=False, fill=False, color='white', max_count=None,
              max_size=None, **param_dict):
    """
    Plot polar histogram of angles on ax. ax must have been created using
    subplot_kw=dict(projection='polar'). Angles are expected in radians.
    """
    # Wrap angles to [-pi, pi)
    angles = (angles + np.pi) % (2*np.pi) - np.pi

    # Set bins symetrically around zero
    if start_zero:
        # To have a bin edge at zero use an even number of bins
        if bins % 2:
            bins += 1
        bins = np.linspace(-np.pi, np.pi, num=bins+1)

    # Bin data and record counts
    count, bin = np.histogram(angles, bins=bins)
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
            area = count / angles.size
        # Calculate corresponding bin radius
        radius = (area / np.pi)**.5
    else:
        radius = count
#    print('radius = ',radius)
#    print('max(radius) = ',np.max(radius))

#    print('maxArea = ',maxArea)

    # Plot data on ax
    ax.bar(bin[:-1], radius, zorder=1, align='edge', width=widths,
           edgecolor='C0', fill=fill, linewidth=1, color=color)
    if max_count and max_size:
        if density is None or density is True:
            max_area = max_count / max_size
            max_radius = (max_area / np.pi)**.5
        else:
            max_radius = max_count
        print('plotting max_radius = ',max_radius)
        ax.bar(0,max_radius,width=0.001)

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels, they are mostly obstructive and not informative
    ax.set_yticks([])

    if lab_unit == "radians":
        label = ['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$',
                  r'$\pi$', r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$']
        ax.set_xticklabels(label)
    return maxCount,angles.size

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

#print('data.getData("EPA") = ',data.getData('EPA'))
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
                print('xy = ',xy)
                x = xy.x
                y = xy.y
                print('x = ',x,', y = ',y)

                gpa = float(data.getData('GPA', iLine))
                if gpa < 0.:
                    gpa += 180
                if gpa > 180:
                    gpa -= 180.

                csGlon = hashData.getData('CS_Glon', hashLine)
                csGlat = hashData.getData('CS_Glat', hashLine)

                dist = None
                if (csGlon != '') and (csGlat != ''):
    #                getStarWithMinDist(gaiaData, ra, dec, iStar=0)
    #                c = SkyCoord(ra=float(csRa)*u.degree, dec=float(csDec)*u.degree, frame='icrs')
    #                lb = c.galactic
    #                print('dir(lb) = ',dir(lb))
                    print('l = ',l,', b = ',b)
                if (x >= xRange[0]) and (x <= xRange[1]) and (y >= yRange[0]) and (y <= yRange[1]):
                    lbxyGPA.append([l, b, x, y, gpa, int(data.getData('flag', iLine)), csGlon, csGlat, dist, hashData.getData('mainClass', iLine)])

    texFile.write('\\caption{Your caption here} % needs to go inside longtable environment\n')
    texFile.write('\\label{tab:myfirstlongtable}\n')
    texFile.write('\\end{longtable}\n')
    texFile.write('%\\end{table}\n')

    texFile.write('Table \\ref{tab:myfirstlongtable} shows my first longtable.\n')
    texFile.write('\\end{document}\n')

    print('len(lbxyGPA) = ',len(lbxyGPA))
    print('lbxyGPA[:][0] = ',lbxyGPA[:][0])
    l = np.array([i[0] for i in lbxyGPA])
    b = np.array([i[1] for i in lbxyGPA])
    x = np.array([i[2] for i in lbxyGPA])
    y = np.array([i[3] for i in lbxyGPA])
    GPA = np.array([i[4] for i in lbxyGPA])
    flag = np.array([i[5] for i in lbxyGPA])
    mClass = np.array([i[9] for i in lbxyGPA])
    oneFlags = np.array([f <= maxFlag for f in flag])
    print('oneFlags = ',type(oneFlags),': ',oneFlags)

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
    plt.show()
    fig.savefig(os.path.join(imOutPath,'all_pne'+fNameSuffix % (p)), bbox_inches='tight')

    f = plt.figure()
    plt.hist(GPA[oneFlags], bins=36)
    plt.show()
    f.savefig(os.path.join(imOutPath,'histogram'+fNameSuffix % (p)), bbox_inches='tight')

    moments = calcMoments(GPA[oneFlags])
    print('moments = ',moments)

    print('np.max(GPA[oneFlags]) = ',np.max(GPA[oneFlags]))

    pltArr = []
    pltMClass = []
    nEllipticals = 0
    nBipolars = 0
    for i in np.arange(0,len(oneFlags),1):
        if oneFlags[i]:
            pltArr.append(math.radians(GPA[i]))
            pltArr.append(math.radians(GPA[i]+180.))
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
    max_count, max_size = rose_plot(ax, pltArr, offset=np.pi/2., bins=36, start_zero=True, color='blue', fill=True)
    rose_plot(ax, pltArr[pltMClass], offset=np.pi/2., bins=36, start_zero=True, color='green', fill=True, max_count=max_count, max_size=max_size)
    fig.tight_layout()
    plt.show()
    fig.savefig(os.path.join(imOutPath,'rose_plot'+fNameSuffix % (p)), bbox_inches='tight')

    xyBulge = ham.lonLatToXY(-15., -5.)
    print('xyBulge = ',xyBulge.x,xyBulge.y)