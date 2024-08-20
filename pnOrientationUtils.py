import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import math

from astropy import units as u
#from astropy.coordinates import SkyCoord
from astropy.stats import circstats
import pycircstat
import seaborn as sns
import subprocess

from PIL import Image, ImageDraw

from myUtils import plotLBMarks

show = True


def plotHammerProjection(x,y,GPA,oneFlags,eFlags,bFlags,xRange,yRange,fNameOut=None):
    xRangeMax = [-2.83,2.83]
    yRangeMax = [-1.4143,1.4143]
    fig = plt.figure(figsize=(25,10))
    plt.axis('off')
    colormap = plt.cm.get_cmap('hsv')
    if len(GPA[oneFlags][eFlags] > 0):
        colorsE = colormap((GPA[oneFlags][eFlags] - np.min(np.array(GPA[oneFlags][eFlags]))) / np.max(np.array(GPA[oneFlags][eFlags])))
        print('GPA[oneFlags][eFlags] = ',GPA[oneFlags][eFlags])
        print('colorsE = ',colorsE)
        plt.scatter(x[oneFlags][eFlags],y[oneFlags][eFlags],c=colorsE,marker='o',s=20)
    if len(GPA[oneFlags][bFlags]) > 0:
        print('GPA[oneFlags][bFlags] = ',GPA[oneFlags][bFlags])
        colorsB = colormap((GPA[oneFlags][bFlags] - np.min(np.array(GPA[oneFlags][bFlags]))) / np.max(np.array(GPA[oneFlags][bFlags])))
        print('colorsB = ',colorsB)
        plt.scatter(x[oneFlags][bFlags],y[oneFlags][bFlags],c=colorsB,marker='d',s=20)
#    plt.scatter(x[oneFlags][eFlags],y[oneFlags][eFlags],c=GPA[oneFlags][eFlags],marker='o',s=20,cmap='hsv')
#    plt.scatter(x[oneFlags][bFlags],y[oneFlags][bFlags],c=GPA[oneFlags][bFlags],marker='d',s=20,cmap='hsv')
    sm = plt.cm.ScalarMappable(cmap=colormap)
    plotLBMarks(10)
    cbarTicks = np.arange(0.,1.0001,20./180.)
    cbar = plt.colorbar(sm,ticks=cbarTicks)
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
