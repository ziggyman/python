import glob
import numpy as np
import os
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

import csvData
import csvFree
from galaxyMath import raDecToLB, parallaxToDistance
from hammer import Pixel,XY,LonLat,Hammer
from myUtils import getStarWithMinDist

path = '/Users/azuri/daten/uni/HKU/PN alignment'
dataFileName = os.path.join(path, 'PN-alignments.csv')
hashFileName = os.path.join(path, 'HASH_bipolar+elliptical_true_PNe.csv')
latexFileName = os.path.join(path, 'PN-alignments.tex')

data = csvFree.readCSVFile(dataFileName)
hashData = csvFree.readCSVFile(hashFileName)

ham = Hammer()
pixels = ham.getPixels()
lbxyGPA = [] #Glon, Glat, x, y, GPA, flag, csGlon, csGlat, [dist,...]

def findInHash(hashData, hashID):
    ids = hashData.getData('idPNMain')
    for iLine in np.arange(0,hashData.size(),1):
        if ids[iLine] == hashID:
            return iLine

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
                getStarWithMinDist(gaiaData, ra, dec, iStar=0)
#                c = SkyCoord(ra=float(csRa)*u.degree, dec=float(csDec)*u.degree, frame='icrs')
#                lb = c.galactic
#                print('dir(lb) = ',dir(lb))
                print('l = ',l,', b = ',b)
            lbxyGPA.append([l, b, x, y, gpa, int(data.getData('flag', iLine)), csGlon, csGlat, dist])

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
    oneFlags = [f < 4 for f in flag]
    print('oneFlags = ',oneFlags)
    plt.scatter(x[oneFlags],y[oneFlags],c=GPA[oneFlags],s=20,cmap='viridis')
    cbar = plt.colorbar()
    cbar.set_label('Galactic Position Angle')
    plt.show()
