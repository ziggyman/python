print('starting to import things')
from astroquery.vo_conesearch.conesearch import conesearch
import csv
import numpy as np
import os
import shutil
import math
import matplotlib.pyplot as plt
from PIL import Image
import pyneb as pn
from pyneb.utils.FortranRecordReader import T
from scipy.optimize import curve_fit
import subprocess
from iphas_paper_calculate_uncertainties import getPNGsfromHashIDs

import csvData
import csvFree
from hammer import Pixel,XY,LonLat,Hammer

#from drUtils import applyVRadCorrection
from myUtils import hmsToDeg, dmsToDeg, getImageData, getWavelength, getHeader,plotLBMarks,applyVRadCorrection,angularDistancePyAsl
from fits_fit_2gauss import gauss,gauss2,gauss3,getAreaGauss,getAreas2Gauss,getAreas3Gauss
from iphas_paper_calculate_uncertainties import calculateErrors#(spectrumFileName,idPNMain,csvLinesFileName)

print('imports done')

imPath = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/iphas-gtc-images'
latexPath = '/Users/azuri/entwicklung/tex/IPHAS-GTC'
publicationImagesPath = os.path.join(latexPath,'images')
spectraPath = '/Users/azuri/spectra/GTC'

surveys = [
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_shs.csv','/data/kegs/pngPNImages/%s/SHS/%s_thre*.png',['idPNMain','idPNMain'],'SHS'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_quotHaSr.csv','/data/kegs/pngPNImages/%s/SHS/%s_quot*.png',['idPNMain','idPNMain'],'quot'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_iphas.csv','/data/kegs/pngPNImages/%s/IPHAS/%s_r%s_iphas3colour*.png',['idPNMain','idPNMain','run_id'],'IPHAS'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_iquotHaSr.csv','/data/kegs/pngPNImages/%s/IPHAS/%s_r%s_iquot*.png',['idPNMain','idPNMain','run_id'],'iquot'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_vphas.csv','/data/kegs/pngPNImages/%s/VPHASplus/%s_vp*.png',['idPNMain','idPNMain'],'VPHAS'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_vquotHaSr.csv','/data/fermenter/PNImages/%s/VPHASplus/quot_vphas*',['idPNMain'],'vquot'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_sss.csv','/data/kegs/pngPNImages/%s/SSS_irb/%s_sss_irb*.png',['idPNMain','idPNMain'],'SSS_irb'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_sss.csv','/data/kegs/pngPNImages/%s/SSS_irb2/%s_sss_irb2*.png',['idPNMain','idPNMain'],'SSS_irb2'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_2mass.csv','/data/kegs/pngPNImages/%s/TWOMASS/%s_2mass*.png',['idPNMain','idPNMain'],'2MASS'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_nvss.csv','/data/kegs/pngPNImages/%s/NVSS/%s_nvss*.png',['idPNMain','idPNMain'],'NVSS'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_galex_rgb.csv','/data/kegs/pngPNImages/%s/Galex/%s_galex_rgb*.png',['idPNMain','idPNMain'],'GALEX_rgb'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_galex_nd.csv','/data/kegs/pngPNImages/%s/GALEX%s_galex_nd*.png',['idPNMain','idPNMain'],'GALEX_nd'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_wise.csv','/data/kegs/pngPNImages/%s/WISE/%s_wise321*.png',['idPNMain','idPNMain'],'WISE321'],
           ['/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_wise.csv','/data/kegs/pngPNImages/%s/WISE/%s_wise432*.png',['idPNMain','idPNMain'],'WISE432'],
           ]

diags = pn.Diagnostics()
# include in diags the relevant line ratios
diags.addDiag([
            '[NII] 5755/6584',
            '[NII] 5755/6548',
            '[NII] 5755/6584+',
            '[OIII] 4363/5007',
            '[SII] 6731/6716',
            ])
diags.addClabel('[SII] 6731/6716', '[SII]a')

# Tell PyNeb tu use parallelisation
pn.config.use_multiprocs()

### General settings
# Setting verbosity level. Enter pn.my_logging? for details
pn.log_.level = 1 # set this to 3 to have more details

c0 = 299792.458 #km/s
o3 = pn.Atom('O', 3)
s2 = pn.Atom('S', 2)
n2 = pn.Atom('N', 2)

linesOfInterest = {'SIIa':6716.44,
                    'SIIb':6730.82,
                    'Halpha':6562.801,
                    'NII':5754.59,
                    'NIIa':6548.05,
                    'NIIb':6583.45,
                    'OIIIa':4363.209,
                    'OIIIb':5006.843	,
                    'Hbeta':4861.363,
                    }

def coneSearch(csv):
    searches = []
    notFound = []
    for i in np.arange(0,csv.size(),1):
        try:
            search = conesearch(center=(hmsToDeg(csv.getData(' RA ',i).strip()), dmsToDeg(csv.getData(' DEC ',i).strip())),
                                radius=0.001,#in degrees
                                verb=3,
                                catalog_db="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&")
            print('dir(search) = ',dir(search))
            print('search.nrows() = ',search.nrows)
            print('search.array() = ',search.array)
            print('search.fields() = ',search.fields)
            print('search.params() = ',search.params)
            print('search.to_table() = ',search.to_table())
#            search.to_table().write('/Volumes/work/azuri/spectra/IPHAS_GTC/iphas-data.fits', format='fits')
            searches.append([csv.getData('Name ',i),search])
        except Exception as e:
            print('Exception caught: ',e)
            print(csv.getData('Name ',i)," not found")
            notFound.append(csv.getData('Name ',i))
            pass
    print('searches = ',len(searches),': ',searches)
    print('not found: ',len(notFound),': ',notFound)

def findHASHid(csvTablePaper, csvTableTargets):
    print('plotLBMarks starting')
    ham = Hammer()
    fig = plt.figure(figsize=(25,10))
    plotLBMarks(10.)
    print('plotLBMarks finished')
    pnMain = csvFree.readCSVFile(os.path.join(imPath[:imPath.rfind('/')],'PNMain.csv'))
    tbCNames = csvFree.readCSVFile(os.path.join(imPath[:imPath.rfind('/')],'tbCNames.csv'))

    pnStat = pnMain.getData('PNstat')
    idx = []#np.where(np.array(pnStat).all() in ['T','L','P'])[0]
    for i in range(pnMain.size()):
        if (pnMain.getData('PNstat',i) in ['T','L','P']) and (pnMain.getData('domain',i) == 'Galaxy'):
            idx.append(i)
    print('idx = ',len(idx),': ',idx)
    xsT = []
    xsL = []
    xsP = []
    ysT = []
    ysL = []
    ysP = []
    for i in idx:
        lon = float(pnMain.getData('Glon',i))
        lat = float(pnMain.getData('Glat',i))
        xy = ham.lonLatToXY(lon,lat)
        x = xy.x
        y = xy.y
        if pnStat[i] == 'L':
            xsL.append(x)
            ysL.append(y)
        elif pnStat[i] == 'P':
            xsP.append(x)
            ysP.append(y)
        else:
            xsT.append(x)
            ysT.append(y)
    print('found ',len(xsT),' True, ',len(xsL),' Likely, and ',len(xsP),' Possible PNe')
    plt.scatter(xsL,ysL,c='g',s=5,marker='v')
    plt.scatter(xsP,ysP,c='r',s=5,marker='s')
    plt.scatter(xsT,ysT,c='b',s=5,marker='o')
    plt.axis('off')
    ids = []
    for iObs in range(csvTablePaper.size()):
        name = csvTablePaper.getData("Name ",iObs).strip(' ')
        if name == 'IPHASX J014238+600947':
            name = 'We 2-5'
        elif name == 'IPHASX J031058.8+624755':
            name = 'Sh 2-200'
        elif name == 'Pa 21':
            name = 'DSH J192315+270734'
        elif name == 'Pa 27':
            name = 'DSH J204858+321815'
        elif name == 'Pa 22':
            name = 'DSH J195813+395440'
        elif name == 'Pa 29':
            name = 'DSH J205943+345423'
        elif name == 'Kn 62':
            name = 'DSH J062355+381515'
        elif name == 'Pa 15':
            name = 'DSH J202907+231108'
        elif name == 'IPHASX J191104.8+060845':
            name = 'IRAS 19086+0603'
#        elif name == ''
#        print('Name =<'+name+'>')
        found = False
        for iTarget in range(csvTableTargets.size()):
            targetName = csvTableTargets.getData('Name',iTarget)
#            print('checking targetName <'+targetName+'>')
            if name == targetName:
#                print('found it!')
                for iName in range(tbCNames.size()):
                    if (tbCNames.getData('idPNMain',iName) == csvTableTargets.getData('idPNMain', iTarget)) & (tbCNames.getData('InUse',iName) == '1'):
                        name = tbCNames.getData('Name',iName)
                        csvTablePaper.setData("Name ",iObs,name)
                ids.append([name,csvTableTargets.getData('idPNMain', iTarget)])
                found = True
        if not found:
            print('PROBLEM: did not find object named <'+name+'>')
            print('object data: ',csvTablePaper.getData(iObs))
        lon = float(csvTablePaper.getData(' l ',iObs))
        lat = float(csvTablePaper.getData(' b ',iObs))
        xy = ham.lonLatToXY(lon,lat)
        x = xy.x
        y = xy.y
        stat = csvTablePaper.getData(' Status ',iObs)
        print('id = ',ids[len(ids)-1],': stat = ',stat,', lon = ',lon,', lat = ',lat,', x = ',x,', y = ',y)
        if stat == ' T ':
            color = 'b'
        elif stat == ' L ':
            color = 'g'
        else:
            color = 'r'
        plt.scatter(x,y,c=color,s=60,marker='D', edgecolors='k',)
    if len(ids) == csvTablePaper.size():
        print('FOUND ALL TARGETS!')
    else:
        print('HMMM, NOT ALL TARGETS FOUND... :(')
    pdfNameTemp = os.path.join(latexPath,'images/hammer_tmp.pdf')
    plt.savefig(pdfNameTemp,bbox_inches='tight')
    #plt.show()
    plt.close()
#    subprocess.run(["gs","-sDEVICE=pdfwrite","-dCompatibilityLevel=1.4","-dPDFSETTINGS=/ebook","-dNOPAUSE", "-dQUIET", "-dBATCH", "-sOutputFile="+pdfNameTemp[:pdfNameTemp.rfind('_tmp')]+'.pdf', pdfNameTemp])
#    subprocess.run(["rm",pdfNameTemp])

#    """Set name from tbCNames"""
#    tbCNames = csvFree.readCSVFile(os.path.join(imPath[:imPath.rfind('/')],'tbCNames.csv'))
#    for id in ids:
#        idPNMain = id[1]
#        found = False
#        for i in range(tbCNames.size()):
#            if (tbCNames.getData('idPNMain',i) == idPNMain) & (tbCNames.getData('InUse',i) == '1'):
#                id[0] = tbCNames.getData('Name',i)
#                found = True
#        if not found:
#            print('Hmmm, did not find idPNMain ',idPNMain,' in tbCNames with an InUse name')
#            STOP

    sortedIndices = np.argsort(np.array([nam.lower() for nam in csvTablePaper.getData('Name ')]))
    print('sortedIndices = ',type(sortedIndices),': ',sortedIndices)
    print('sortedIndices[0] = ',type(sortedIndices[0]),': ',sortedIndices[0])

    csvTablePaper.sort(sortedIndices)
    csvFree.writeCSVFile(csvTablePaper,os.path.join(latexPath,'table_true_names_sorted.tex'),'&')

    return [ids[i] for i in sortedIndices]

def getImages(ids):
    print('getImages: ids = ',len(ids),': ',ids)
    hashPNMain = '/Users/azuri/daten/uni/HKU/HASH/PNMain_full_Aug-07-2020.csv'

    outFile = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/getIPHASImages.bash'
    with open(outFile,'w') as f:
        f.write('mkdir iphas-gtc-images\n')
        for survey in surveys:
            csv = csvFree.readCSVFile(survey[0])
            for id in ids:
                hashID = id[1]
                lineIDs = csv.find('idPNMain',hashID)
                print('getImages: hashID = ',hashID,': lineIDs = ',lineIDs)
                if lineIDs[0] == -1:
                    print('getImages: PROBLEM: did not find hashID ',hashID,' in ',survey[0])
                    #STOP
                else:
                    for lineID in lineIDs:
                        if csv.getData('inuse',lineID) == '1':
                            print('survey[1] = ',survey[1],' survey[2] = ',survey[2])
                            tmp = [csv.getData(num,lineID) for num in survey[2]]
                            print('tmp = ',tmp)
#                            print('dir(tmp) = ',dir(tmp))
                            print('(num for num in survey[2]) = ',(num for num in tmp))
                            print('survey = ',survey)
                            print('survey[1] = ',survey[1],' % (',tuple([csv.getData(num,lineID) for num in survey[2]]),')')
                            print('[csv.getData(num,lineID) for num in survey[2]] = ',[csv.getData(num,lineID) for num in survey[2]])
                            f.write('cp '+survey[1] % tuple([csv.getData(num,lineID) for num in survey[2]])+' iphas-gtc-images/\n')

def filterIphasImagesInBashFile():
    bashFile = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/getIPHASImages.bash'
    with open(bashFile,'r') as f:
        lines = f.readlines()
    with open(bashFile[:bashFile.rfind('.')]+'_iphas.bash','w') as f:
        for line in lines:
            if 'IPHAS' in line:
                f.write(line)

def filterSSSImagesInBashFile():
    bashFile = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/getIPHASImages.bash'
    with open(bashFile,'r') as f:
        lines = f.readlines()
    with open(bashFile[:bashFile.rfind('.')]+'_sss.bash','w') as f:
        for line in lines:
            if 'SSS' in line:
                f.write(line)

def filterWISE432ImagesInBashFile():
    bashFile = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/getIPHASImages.bash'
    with open(bashFile,'r') as f:
        lines = f.readlines()
    with open(bashFile[:bashFile.rfind('.')]+'_wise432.bash','w') as f:
        for line in lines:
            if 'wise432' in line:
                f.write(line)

def filterNVSSImagesInBashFile():
    bashFile = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/getIPHASImages.bash'
    with open(bashFile,'r') as f:
        lines = f.readlines()
    with open(bashFile[:bashFile.rfind('.')]+'_nvss.bash','w') as f:
        for line in lines:
            if 'nvss' in line:
                f.write(line)

def filterGalexImagesInBashFile():
    bashFile = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/getIPHASImages.bash'
    with open(bashFile,'r') as f:
        lines = f.readlines()
    with open(bashFile[:bashFile.rfind('.')]+'_galex.bash','w') as f:
        for line in lines:
            if 'galex' in line:
                f.write(line)

def combineImages(idPNMain):
    (_, _, filenames) = next(os.walk(imPath))
    for survey in surveys:
        print("survey[1][survey[1].rfind('/')+1:] = ",survey[1][survey[1].rfind('/')+1:])
        print("survey[2][1:] = ",survey[2][1:])
        print("tuple([idPNMain[1] for i in survey[2][1:]]) = ",tuple(survey[2][1:]))
        filename = survey[1][survey[1].rfind('/')+1:]
        filenameSearch = filename[filename.rfind('%s')+2:filename.find('*')]
        print('idPNMain = ',idPNMain,': looking for filenameSearch <'+filenameSearch+'> in filesnames')
        imagenames = []
        for filename in filenames:
            filename = filename[filename.rfind('/')+1:]
            if (filename[:filename.find('_')] == str(idPNMain)) and (filenameSearch in filename):
                print('found filename <'+filename+'>')
                imagenames.append(filename)
        if len(imagenames) > 1:
            background = ''
            foreground = []
            lens = [len(imagename) for imagename in imagenames]
            minLen = min(lens)
            for i in range(len(imagenames)):
                if lens[i] == minLen:
                    background = os.path.join(imPath,imagenames[i])
                else:
                    if ('overlay' not in imagenames[i]) and ('combined' not in imagenames[i]):
                        foreground.append(os.path.join(imPath,imagenames[i]))
            backgroundIm = Image.open(background)
            print('combining ',background,' and ',foreground)
            for f in foreground:
                f = Image.open(f)
                backgroundIm.paste(f, (0, 0), f)
            backgroundIm.save(os.path.join(publicationImagesPath,background[background.rfind('/')+1:background.rfind('.')]+'_combined.png'))
    #if idPNMain == '8214':
    #    STOP

def createImageTable(ids):
    (_, _, filenames) = next(os.walk(publicationImagesPath))
#    print('filenames = ',filenames)
    with open(os.path.join(latexPath,'image_table.tex'),'w') as f:
        f.write('\\documentclass[12pt]{article}\n')
        f.write('\\usepackage{graphicx}\n')
        f.write('\\usepackage{natbib}\n')
        f.write('\\usepackage{longtable}\n')
        f.write('\\usepackage{pdflscape}\n')
        f.write('\\setcitestyle{numbers}\n')
        f.write('\\setcitestyle{square}\n')
        f.write('\\begin{document}\n')
        f.write('\\begin{landscape}\n')
        f.write('\\clearpage\n')
        f.write('\\onecolumn\n')
        f.write('\\begin{longtable}{ | *{7}{l|} }\n')
        f.write('\\hline\n')
        f.write('Target Name & optical & $\\mathrm{H_\\alpha/Sr}$ & WISE321 & NVSS & GALEX & spectrum\\\\\n')
        f.write('\\endhead  % header material\n')
        f.write('\\hline\\endfoot  % footer material\n')
        f.write('\\hline\n')
        for id in ids:
            print('id = ',id)
            imNames = [None,None,None,None,None,None]
            idPNMain = id[1]
            for survey in surveys:
                for filename in filenames:
                    if filename[:filename.find('_')] == idPNMain:
                        fNameSur = survey[1][survey[1].rfind('%s')+2:survey[1].rfind('*')]
                        if '/' in fNameSur:
                            fNameSur = fNameSur[fNameSur.rfind('/')+1:]
                        print('fNameSur = <'+fNameSur+'>')
                        if fNameSur in filename:
                            if survey[3] == 'SHS':
                                imNames[0] = filename
                            elif survey[3] == 'quot':
                                imNames[1] = filename
                            elif survey[3] == 'IPHAS':
                                imNames[0] = filename
                            elif survey[3] == 'iquot':
                                imNames[1] = filename
                            elif survey[3] == 'VPHAS':
                                imNames[0] = filename
                            elif survey[3] == 'vquot':
                                imNames[1] = filename
                            elif survey[3] == 'WISE321':
                                imNames[2] = filename
                            elif survey[3] == 'NVSS':
                                imNames[3] = filename
                            elif survey[3] == 'GALEX_nd':
                                imNames[4] = filename
                            elif survey[3] == 'GALEX_rgb':
                                imNames[4] = filename
                            elif survey[3] == 'SSS_irb2':
                                if imNames[0] is None:
                                    imNames[0] = filename
                            elif survey[3] == 'SSS_irb':
                                if imNames[0] is None:
                                    imNames[0] = filename
                        imNames[5] = id[2]

#                            else:
#                                if len(imNames) < 4:
#                                    if filename not in imNames:
#                                        imNames.append(filename)
            print('id = ',id,': imNames = ',imNames)
            f.write('%s & ' % (id[0]))
            for imName in imNames[:5]:
                if imName is None:
                    f.write('& ')
                else:
                    f.write('\\includegraphics[width=2.2cm,height=2cm]{%s} &' % (os.path.join('images',imName)))
            if imNames[5] is None:
                f.write('\\\\\n')
            else:
                f.write('\\includegraphics[width=2.3cm,height=2cm]{%s}\\\\\n' % (os.path.join('images',imNames[5])))
            f.write('\\hline\n')
        f.write('\\end{longtable}\n')
        f.write('\\end{landscape}\n')
        f.write('\\clearpage\n')
        f.write('\\twocolumn\n')
        f.write('\\end{document}\n')

def fixInUseInIquote():
    iphasCSV = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_iphas.csv')
    iquotCSV = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/PNImages_iquotHaSr.csv')
    idPNe = iphasCSV.getData('idPNMain')
    print('len(idPNe) = ',len(idPNe))
    idPNeUnique = np.unique(idPNe)
    print('len(idPNeUnique) = ',len(idPNeUnique))
    problematic = []
    linesOut = []
    idsOut = []
    for idPN in idPNeUnique:
        lineIDsIphas = iphasCSV.find('idPNMain',idPN,0)
        print('idPNMain = ',idPN,': lineIDsIphas = ',lineIDsIphas)
        lineIDsIquot = iquotCSV.find('idPNMain',idPN,0)
        print('idPNMain = ',idPN,': lineIDsIquot = ',lineIDsIquot)
        for lineIDIphas in lineIDsIphas:
            if (iphasCSV.getData('inuse',lineIDIphas) == '1') and (iquotCSV.find('run_id',iphasCSV.getData('run_id',lineIDIphas))[0] == -1):
                print('ERROR: could not find run_id ',iphasCSV.getData('run_id',lineIDIphas),' in iquot')
                thisProblem = [iphasCSV.getData('idPNMain',lineIDIphas),iphasCSV.getData('run_id',lineIDIphas)]
                if thisProblem not in problematic:
                    problematic.append(thisProblem)
                #STOP
            if lineIDsIquot[0] != -1:
                for lineIDIquot in lineIDsIquot:
                    if iphasCSV.getData('run_id',lineIDIphas) == iquotCSV.getData('run_id',lineIDIquot):
                        if iphasCSV.getData('inuse',lineIDIphas) != iquotCSV.getData('inuse',lineIDIquot):
                            line = 'UPDATE `iquot_HaSr` SET `InUse` = %d WHERE `idiquot_HaSr` = %d;\n' % (int(iphasCSV.getData('inuse',lineIDIphas)),
                                                                                                           int(iquotCSV.getData('idiquot_HaSr',lineIDIquot)))
                            if line not in linesOut:
                                if iphasCSV.getData('idPNMain',lineIDIphas) not in idsOut:
                                    idsOut.append(iphasCSV.getData('idPNMain',lineIDIphas))
                                linesOut.append(line)
    with open('/Users/azuri/daten/uni/HKU/HASH/fixInUseInIquote.sql','w') as f:
        f.write('USE PNImages;\n')
        for line in linesOut:
            f.write(line)
    print('idsOut = ',idsOut)
    print('problematic = ',len(problematic),': ',problematic)

#        for i in range(iquoteCSV.size()):
#
#        for i in range(iphasCSV.size()):

def splitObservationFile(obsFileIn):
    with open(obsFileIn,'r') as f:
        lines = f.readlines()
    for line in lines[1:]:
        with open(obsFileIn[:obsFileIn.rfind('.')]+'_'+line[:line.find('\t')]+obsFileIn[obsFileIn.rfind('.'):],'w') as f:
            f.write(lines[0])
            f.write(line)

def makeSpectraTable(ids, hashPNMainFileName, calculateLineIntensities = False):
    (_, _, filenames) = next(os.walk(spectraPath))
    print('filenames = ',filenames)
    #STOP
    csvHashPNMain = csvFree.readCSVFile(hashPNMainFileName)

    areasSII6716 = []
    areasSII6731 = []
    areasHalpha = []
    areasHbeta = []
    areasNII5755 = []
    areasNII6548 = []
    areasNII6583 = []
    areasOIII4363 = []
    areasOIII5007 = []

    csvOut = csvData.CSVData()
    csvOut.header = ['idPNMain',
                     '$\mathrm{H_\\alpha}$',
                     '$\mathrm{H_\\beta}$',
                     '$\mathrm{[SII]_{6716}}$',
                     '$\mathrm{[SII]_{6731}}$',
                     '$\mathrm{[NII]_{5755}}$',
                     '$\mathrm{[NII]_{6548}}$',
                     '$\mathrm{[NII]_{6583}}$',
                     '$\mathrm{[OIII]_{4363}}$',
                     '$\mathrm{[OIII]_{5007}}$',
#                     '$\mathrm{T_{e^-}}$',
#                     '$\mathrm{\\rho_{e^-}}$',
                     ]

    has = csvFree.readCSVFile(os.path.join(imPath[:imPath.rfind('/')],'h_a.csv'))

    with open(os.path.join(imPath[:imPath.rfind('/')],'lines.csv'),'w') as fl:
        fl.write('fileName,line,strength\n')
        with open(os.path.join(imPath[:imPath.rfind('/')],'vrad.csv'),'w') as fv:
            fv.write('fileName,vrad\n')
            with open(os.path.join(latexPath,'spectra_table.tex'),'w') as f:
                f.write('\\documentclass[12pt]{article}\n')
                f.write('\\usepackage{graphicx}\n')
                f.write('\\usepackage{natbib}\n')
                f.write('\\usepackage{longtable}\n')
                f.write('\\setcitestyle{numbers}\n')
                f.write('\\setcitestyle{square}\n')

                f.write('\\begin{document}\n')

                f.write('\\begin{longtable}{ | *{3}{l|} }\n')
                f.write('\\hline\n')
        #        f.write('Target Name & image 1 & image 2 & image 3 & image 4\\\\\n')
                f.write('\\endhead  % header material\n')
                f.write('\\hline\\endfoot  % footer material\n')
                f.write('\\hline\n')

                idPNMains = [id[1] for id in ids]
                nSpec = 0
                for filename in filenames:
                    if filename[filename.rfind('.'):] == '.fits':
                        idPNMain = None
                        if filename == 'Pa30_GT080716.fits':
                            idPNMain = '15569'
                        elif filename == 'IPHASXJ230323_GT200816_t.fits':
                            idPNMain = '8248'
                        else:
                            hashFitsFiles = csv.DictReader(open('/Users/azuri/daten/uni/HKU/IPHAS-GTC/fitsfiles.csv'))
                            for hashFitsFile in hashFitsFiles:
                                if hashFitsFile['fileName'] == filename:
                                    idPNMain = hashFitsFile['idPNMain']
                        if idPNMain is None:
                            print('could not find filename ',filename,' in hashFitsFiles')
                            STOP
                        print('filename = ',filename,': idPNMain = ',idPNMain)
                        if idPNMain in idPNMains:
                            pdfFileName = os.path.join(latexPath,'images',filename[:filename.rfind('.')].replace('.','_')+'.pdf')
                            print('filename = ',filename,': idPNMain = ',idPNMain,' is a good one')
                            doIt = False
                            for id in ids:
                                if id[1] == idPNMain:
                                    objectName = id[0]
                                    png = getPNGsfromHashIDs(csvHashPNMain,[idPNMain])[0]
                                    print('idPNMain = ',idPNMain,': png = ',png)

                                    if filename in ['LDu1_sum.fits',
                                                    'We2-260_GT210816.fits',
                                                    'Kn24_GT230616.fits',
                                                    'IPHASXJ194301_GT120417.fits',
                                                    'We2-5_GT260816.fits',
                                                    'Ou2_GT140816.fits',
                                                    ]:
                                        doIt = True
#                                        if len(id) > 2:
#                                            id[2] = pdfFileName[pdfFileName.rfind('/')+1:]
#                                        else:
                                        id.append(pdfFileName[pdfFileName.rfind('/')+1:])
                                    else:
                                        if filename[:filename.find('_')] not in ['LDu1',
                                                                                 'We2-260',
                                                                                 'Kn24',
                                                                                 'IPHASXJ194301',
                                                                                 'We2-5',
                                                                                 'Ou2',
                                                                                ]:
                                            doIt = True
                                            id.append(pdfFileName[pdfFileName.rfind('/')+1:])
                            if doIt:
                                f.write('\\includegraphics[width=3.5cm,height=3cm]{%s}\n ' % (os.path.join('images',pdfFileName[pdfFileName.rfind('/')+1:])))
                                if nSpec < 3:
                                    f.write(' & ')
                                if nSpec == 2:
                                    f.write('\\\\')
                                    f.write('    \\hline\n')
                                    nSpec = 0
                                else:
                                    nSpec += 1
                                if calculateLineIntensities:
                                    lam = getWavelength(getHeader(os.path.join(spectraPath,filename),0),axis=1)
                                    flux = getImageData(os.path.join(spectraPath,filename),0)
                                    idx = np.where(lam < 7350.)[0]
                                    idx = np.where(lam[idx] > 3950.)[0]
                                    plt.plot(lam[idx], flux[idx], 'k-')
                                    plt.xlabel('wavelength [$\mathrm{\AA}$]')
                                    plt.ylabel('$\mathrm{F_\lambda}$ [$\mathrm{ergs s^{-1} cm^{-2} \AA^{-1}}$]')
                                    plt.title(png)#objectName)
                                    pdfFileName = pdfFileName[:-4].replace('.','_')+'.pdf'
                                    plt.savefig(pdfFileName,bbox_inches='tight')
                                    plt.close()
                                    print('pdf saved to <'+pdfFileName+'>')
#                                    if 'Ou3' in pdfFileName:
#                                        STOP
                                    hasPos = has.find('fileName',filename,0)
                                    print('hasPos = ',hasPos)
                                    show = False
                                    addx = 10.
                                    sigma = 3.
                                    if hasPos[0] != -1:
                                        HalphaPos = float(has.getData('H_alpha',hasPos[0]))
                                        vrad = c0 * (HalphaPos - linesOfInterest['Halpha'])/ linesOfInterest['Halpha']
                                        fv.write(filename+','+str(vrad)+'\n')
                                        print('HalphaPos = ',HalphaPos,': vrad = ',vrad)
                                        lam = applyVRadCorrection(lam,vrad)
                                        #plt.xlim(linesOfInterest['Halpha'] - 30.,linesOfInterest['Halpha'] + 30.)
                                        if True:
                                            indices = np.where(abs(lam-linesOfInterest['Halpha'])<9.)[0]
                                            print('indices = ',indices)
                                            try:
                                                plotName = os.path.join(imPath[:imPath.rfind('/')],'lines/'+idPNMain+'_Ha.png')
                                                areaHalpha,popt = getAreaGauss(lam[indices],
                                                                               flux[indices],
                                                                               np.max(flux[indices]),
                                                                               linesOfInterest['Halpha'],
                                                                               sigma,
                                                                               addOnBothSidesOfX=addx,
                                                                               show=show,
                                                                               save=plotName)
                                                if areaHalpha < 0. or abs(popt[1] - linesOfInterest['Halpha']) > 3. or abs(abs(popt[2]) - sigma) > 1.:
                                                    areaHalpha = 0.
                                                    os.remove(plotName)
                                            except:
                                                areaHalpha = 0.
                                            fl.write(filename+',Halpha,%s\n' % (str(areaHalpha)))
                                            indices = np.where(abs(lam-linesOfInterest['Hbeta'])<10.)[0]
                                            print('indices = ',indices)
                                            try:
                                                plotName = os.path.join(imPath[:imPath.rfind('/')],'lines/'+idPNMain+'_Hb.png')
                                                areaHbeta,popt = getAreaGauss(lam[indices],
                                                                              flux[indices],np.max(flux[indices]),
                                                                              linesOfInterest['Hbeta'],
                                                                              sigma,
                                                                              addOnBothSidesOfX=addx,
                                                                              show=show,
                                                                              save=plotName)
                                                if idPNMain == '4386':
                                                    print('idPNMain = 4386: areaHbeta = ',areaHbeta,', popt = ',popt)
                                                    STOP
                                                if areaHbeta < 0. or abs(popt[1]-linesOfInterest['Hbeta']) > 3. or abs(abs(popt[2]) - sigma) > 1.:
                                                    areaHbeta = 0.
                                                    os.remove(plotName)
                                            except:
                                                areaHbeta = 0.
                                            fl.write(filename+',Hbeta,%s\n' % (str(areaHbeta)))

                                            indices = np.where(abs(lam-linesOfInterest['OIIIa'])<10.)[0]
                                            try:
                                                plotName = os.path.join(imPath[:imPath.rfind('/')],'lines/'+idPNMain+'_OIII4363.png')
                                                areaOIII4363,popt = getAreaGauss(lam[indices],
                                                                                flux[indices],
                                                                                np.max(flux[indices]),
                                                                                linesOfInterest['OIIIa'],
                                                                                sigma,
                                                                                addOnBothSidesOfX=addx,
                                                                                show=show,
                                                                                save=plotName)
                                                if areaOIII4363 < 0. or abs(popt[1]-linesOfInterest['OIIIa']) > 3. or abs(abs(popt[2]) - sigma) > 1.:
                                                    areaOIII4363 = 0.
                                                    os.remove(plotName)
                                            except:
                                                areaOIII4363 = 0.
                                            fl.write(filename+',OIII4363,%s\n' % (str(areaOIII4363)))

                                            indices = np.where(abs(lam-linesOfInterest['OIIIb'])<10.)[0]
                                            try:
                                                plotName = os.path.join(imPath[:imPath.rfind('/')],'lines/'+idPNMain+'_OIII5007.png')
                                                areaOIII5007,popt = getAreaGauss(lam[indices],
                                                                                flux[indices],
                                                                                np.max(flux[indices]),
                                                                                linesOfInterest['OIIIb'],
                                                                                sigma,
                                                                                addOnBothSidesOfX=addx,
                                                                                show=show,
                                                                                save=plotName)
                                                if areaOIII5007 < 0. or abs(popt[1]-linesOfInterest['OIIIb']) > 3. or abs(abs(popt[2]) - sigma) > 1.:
                                                    areaOIII5007 = 0.
                                                    os.remove(plotName)
                                            except:
                                                areaOIII5007 = 0.
                                            fl.write(filename+',OIII5007,%s\n' % (str(areaOIII5007)))

                                            indices = np.where(abs(lam-linesOfInterest['NII'])<10.)[0]
                                            try:
                                                plotName = os.path.join(imPath[:imPath.rfind('/')],'lines/'+idPNMain+'_NII5755.png')
                                                areaNII5755,popt = getAreaGauss(lam[indices],
                                                                                flux[indices],
                                                                                np.max(flux[indices]),
                                                                                linesOfInterest['NII'],
                                                                                sigma,
                                                                                addOnBothSidesOfX=addx,
                                                                                show=show,
                                                                                save=plotName)
                                                if areaNII5755 < 0. or abs(popt[1]-linesOfInterest['NII']) > 3. or abs(abs(popt[2]) - sigma) > 1.:
                                                    areaNII5755 = 0.
                                                    os.remove(plotName)
                                            except:
                                                areaNII5755 = 0.
                                            fl.write(filename+',NII5755,%s\n' % (str(areaNII5755)))

                                            indices = np.where(abs(lam-linesOfInterest['NIIa'])<7.)[0]
                                            try:
                                                plotName = os.path.join(imPath[:imPath.rfind('/')],'lines/'+idPNMain+'_NII6548.png')
                                                areaNII6548,popt = getAreaGauss(lam[indices],
                                                                                flux[indices],
                                                                                np.max(flux[indices]),
                                                                                linesOfInterest['NIIa'],
                                                                                sigma,
                                                                                addOnBothSidesOfX=addx,
                                                                                show=show,
                                                                                save=plotName)
                                                if areaNII6548 < 0. or abs(popt[1]-linesOfInterest['NIIa']) > 3. or abs(abs(popt[2]) - sigma) > 1.:
                                                    areaNII6548 = 0.
                                                    os.remove(plotName)
                                            except:
                                                areaNII6548 = 0.
                                            fl.write(filename+',NII6548,%s\n' % (str(areaNII6548)))

                                            indices = np.where(abs(lam-linesOfInterest['NIIb'])<7.)[0]
                                            try:
                                                plotName = os.path.join(imPath[:imPath.rfind('/')],'lines/'+idPNMain+'_NII6583.png')
                                                areaNII6583,popt = getAreaGauss(lam[indices],
                                                                                flux[indices],
                                                                                np.max(flux[indices]),
                                                                                linesOfInterest['NIIb'],
                                                                                sigma,
                                                                                addOnBothSidesOfX=addx,
                                                                                show=show,
                                                                                save=plotName)
                                                if areaNII6583 < 0. or abs(popt[1]-linesOfInterest['NIIb']) > 3. or abs(abs(popt[2]) - sigma) > 1.:
                                                    areaNII6583 = 0.
                                                    os.remove(plotName)
                                            except:
                                                areaNII6583 = 0.
                                            fl.write(filename+',NII6583,%s\n' % (str(areaNII6583)))

                #                            areaNII6548,areaHalpha,areaNII6583,popt = getAreas3Gauss(lam,flux,np.max(flux[np.where(abs(lam-linesOfInterest['NIIa'])<5.)[0]]),
                #                                                                                              np.max(flux[np.where(abs(lam-linesOfInterest['Halpha'])<5.)[0]]),
                #                                                                                              np.max(flux[np.where(abs(lam-linesOfInterest['NIIb'])<5.)[0]]),
                #                                                                                              linesOfInterest['NIIa'],
                #                                                                                              linesOfInterest['Halpha'],
                #                                                                                              linesOfInterest['NIIb'],
                #                                                                                              10.,
                #                                                                                              10.,
                #                                                                                              10.,
                #                                                                                              show=True)
                                            indices = np.where(lam < linesOfInterest['SIIb'] + 10.)[0]
                                            print('indicesL = ',indices)
                                            indices = np.where(lam[indices] > linesOfInterest['SIIa'] - 10.)[0]
                                            print('indicesG = ',indices)
                #                            plt.plot(lam[indices],flux[indices])
                #                            plt.show()
                                            try:
                                                plotName = os.path.join(imPath[:imPath.rfind('/')],'lines/'+idPNMain+'_SII6716+6731.png')
                                                areaSII6716,areaSII6731,popt = getAreas2Gauss(lam[indices],
                                                                                              flux[indices],
                                                                                              np.max(flux[np.where(abs(lam-linesOfInterest['SIIa'])<5.)[0]]),
                                                                                              np.max(flux[np.where(abs(lam-linesOfInterest['SIIb'])<5.)[0]]),
                                                                                              linesOfInterest['SIIa'],
                                                                                              linesOfInterest['SIIb'],
                                                                                              sigma,
                                                                                              sigma,
                                                                                              addOnBothSidesOfX=addx,
                                                                                              show=show,
                                                                                              save=plotName
                                                                                            )
                                                print('areaSII6716 = ',areaSII6716,', areaSII6731 = ',areaSII6731)
                                                if areaSII6716 < 0. or abs(popt[2]-linesOfInterest['SIIa']) > 3. or abs(abs(popt[4]) - sigma) > 1.:
                                                    areaSII6716 = 0.
                                                if areaSII6731 < 0. or abs(popt[3]-linesOfInterest['SIIb']) > 3. or abs(abs(popt[5]) - sigma) > 1.:
                                                    areaSII6731 = 0.
                                                if (areaSII6716 < 0. or abs(popt[2]-linesOfInterest['SIIa']) > 3.) and (areaSII6731 < 0. or abs(popt[3]-linesOfInterest['SIIb']) > 3.):
                                                    os.remove(plotName)
                                            except:
                                                areaSII6716,areaSII6731 = [0.,0.]
                                            print('areaSII6716 = ',areaSII6716,', areaSII6731 = ',areaSII6731)
                                            fl.write(filename+',SII6716,%s\n' % (str(areaSII6716)))
                                            fl.write(filename+',SII6731,%s\n' % (str(areaSII6731)))
                                            print('areaHalpha = ',areaHalpha)
                                            areasHalpha.append(areaHalpha)
                                            areasHbeta.append(areaHbeta)
                                            areasOIII4363.append(areaOIII4363)
                                            areasOIII5007.append(areaOIII5007)
                                            areasNII5755.append(areaNII5755)
                                            areasNII6548.append(areaNII6548)
                                            areasNII6583.append(areaNII6583)
                                            areasSII6716.append(areaSII6716)
                                            areasSII6731.append(areaSII6731)
        #                                    tem, den = diags.getCrossTemDen(diag_tem='[NII] 5755/6584',
        #                                                                    diag_den='[SII] 6731/6716',
        #                                                                    value_tem=areaNII5755 / (areaNII6548 if areaNII6548 > 0. else 0.00001),
        #                                                                    value_den=areaSII6731 / (areaSII6716 if areaSII6716 > 0. else 0.00001),
        #                                                                   )
        #                                    print('areaNII5755 = ',areaNII5755,', areaNII6548 = ',areaNII6548,', areaSII6731 = ',areaSII6731,', areaSII6716 = ',areaSII6716,': tem = ',tem,', den = ',den)
        #                                    fl.write(filename+',tem,%s\n' % (str(tem)))
        #                                    fl.write(filename+',den,%s\n' % (str(den)))
                                            #csvOut.header = ['idPNMain','Halpha','Hbeta','SII6716','SII6731','NII5755','NII6548','NII6583','OIII4363','OIII5007','$T_{e^-}$','$\rho_{e^-}$']
                                            rowStr = '%s,'
                                            rowStr += ('%0.2E' if areaHalpha > 0. else '%.0f')+','
                                            rowStr += ('%0.2E' if areaHbeta > 0. else '%.0f')+','
                                            rowStr += ('%0.2E' if areaSII6716 > 0. else '%.0f')+','
                                            rowStr += ('%0.2E' if areaSII6731 > 0. else '%.0f')+','
                                            rowStr += ('%0.2E' if areaNII5755 > 0. else '%.0f')+','
                                            rowStr += ('%0.2E' if areaNII6548 > 0. else '%.0f')+','
                                            rowStr += ('%0.2E' if areaNII6583 > 0. else '%.0f')+','
                                            rowStr += ('%0.2E' if areaOIII4363 > 0. else '%.0f')+','
                                            rowStr += ('%0.2E' if areaOIII5007 > 0. else '%.0f')
                                            print('rowStr = <'+rowStr+'>')
                                            rowStr = rowStr % (idPNMain,
                                                            areaHalpha,
                                                            areaHbeta,
                                                            areaSII6716,
                                                            areaSII6731,
                                                            areaNII5755,
                                                            areaNII6548,
                                                            areaNII6583,
                                                            areaOIII4363,
                                                            areaOIII5007,
                                                            )
                                            #if idPNMain == '15569':
                                            #    print('rowStr = ',rowStr)
                                            #    STOP
                                            csvOut.append(rowStr.split(','))
                f.write('\\end{longtable}\n')
                f.write('\\end{document}\n')

    if calculateLineIntensities:
        outFileName = os.path.join(imPath[:imPath.rfind('/')],'lines_id.csv')
        csvFree.writeCSVFile(csvOut,outFileName,'&')
#    STOP
    return csvOut

def getMeanAndStd(arrayIn):
    arr = removeBadValues(arrayIn)
    return [np.mean(arr),np.std(arr)]

def removeBadValues(arrayIn):
    arrayOut = arrayIn[np.where(arrayIn > 0.)[0]]
#    print('removeBadValues: arrayOut = ',arrayOut.shape,', ',arrayOut)
    idx = [not math.isnan(elem) for elem in arrayOut]
#    print('removeBadValues: idx = ',len(idx),', ',idx)
    arrayOut = arrayOut[idx]
#    print('removeBadValues: arrayOut = ',arrayOut.shape,', ',arrayOut)
    idx = [not math.isinf(elem) for elem in arrayOut]
#    print('removeBadValues: idx = ',len(idx),', ',idx)
    arrayOut = arrayOut[idx]
#    print('removeBadValues: arrayOut = ',arrayOut.shape,', ',arrayOut)
    return arrayOut

#def removeBadValues(arrayA,arrayB)

def addLinesToTable(csvPaper,csvLines=None,obsFileName=None):
    idsPaper = csvPaper.getData(' HASH ID ')
    idsPaper = [id.strip() for id in idsPaper]
    print('idsPaper = ',len(idsPaper),': ',idsPaper)
#    STOP

    if csvLines is None:
        csvLines = csvFree.readCSVFile(os.path.join(imPath[:imPath.rfind('/')],'lines_id_corrected.csv'),'&',False)
    print('csvLines.header = ',csvLines.header)
    idsLines = csvLines.getData('idPNMain')
    print('idsLines = ',len(idsLines),': ',idsLines)
#    idsLines.header = ['idPNMain',
#                     '$\mathrm{H_\\alpha}$',
#                     '$\mathrm{H_\\beta}$',
#                     '$\mathrm{[SII]_{6716}}$',
#                     '$\mathrm{[SII]_{6731}}$',
#                     '$\mathrm{[NII]_{5755}}$',
#                     '$\mathrm{[NII]_{6548}}$',
#                     '$\mathrm{[NII]_{6583}}$',
#                     '$\mathrm{[OIII]_{4363}}$',
#                     '$\mathrm{[OIII]_{5007}}$',
    for i in range(csvLines.size()):
        if (float(csvLines.getData('$\mathrm{[SII]_{6716}}$',i)) + float(csvLines.getData('$\mathrm{[SII]_{6731}}$',i)) > 0.) and (float(csvLines.getData('$\mathrm{[NII]_{6548}}$',i)) + float(csvLines.getData('$\mathrm{[NII]_{6583}}$',i)) > 0.):
            plt.scatter(np.log10(float(csvLines.getData('$\mathrm{H_\\alpha}$',i)) / (float(csvLines.getData('$\mathrm{[SII]_{6716}}$',i))+float(csvLines.getData('$\mathrm{[SII]_{6731}}$',i)))),
                        np.log10(float(csvLines.getData('$\mathrm{H_\\alpha}$',i)) / (float(csvLines.getData('$\mathrm{[NII]_{6548}}$',i)) + float(csvLines.getData('$\mathrm{[NII]_{6583}}$',i)))),c='k')
    plt.xlabel(r'log($\mathrm{H_\alpha} / \mathrm{[SII]}$)')
    plt.ylabel(r'log($\mathrm{H_\alpha} / \mathrm{[NII]}$)')
    plt.savefig(os.path.join(latexPath,'images/Halpha_NII_vs_Halpha_SII.pdf'),bbox_inches='tight')
    plt.show()

    csvPaper.removeColumn(' comments')

    csvOut = csvData.CSVData()
    csvOut.header = csvLines.header
    csvOut.renameColumn('idPNMain','Target Name')
    csvOut.data = csvLines.data
    csvOut.addColumn('$\mathrm{(H_\\alpha)/(H_\\beta)}$')
    csvOut.addColumn('$\mathrm{(6717,31)/(H_\\alpha)}$')
    csvOut.addColumn('$\mathrm{(5007)/(H_\\beta)}$')
    csvOut.addColumn('$\mathrm{(6584)/(H_\\alpha)}$')
    csvOut.addColumn('E(B-V)')
    csvOut.addColumn('$\mathrm{T_{e^-}([NII],[SII])}$')
#    csvOut.addColumn('$\mathrm{T_{e^-}([NII],[SII])a}$')
#    csvOut.addColumn('$\mathrm{T_{e^-}([NII],[SII])b}$')
    csvOut.addColumn('$\mathrm{T_{e^-}([OIII],[SII])}$')
    csvOut.addColumn('$\mathrm{\\rho_{e^-}([NII],[SII])}$')
#    csvOut.addColumn('$\mathrm{\\rho_{e^-}([NII],[SII])a}$')
#    csvOut.addColumn('$\mathrm{\\rho_{e^-}([NII],[SII])b}$')
    csvOut.addColumn('$\mathrm{\\rho_{e^-}([OIII],[SII])}$')
#    csvOut.addColumn('$\mathrm{\\rho_{e^-}(SII,T=10.000)}$')
#    csvOut.addColumn('$\mathrm{\\rho_{e^-}(OIII,T=10.000)}$')
#    csvOut.addColumn('$\mathrm{\\rho_{e^-}(NII,T=10.000)}$')

    if obsFileName is None:
        obsFileName = os.path.join(imPath[:imPath.rfind('/')],'observation.dat')
        with open(obsFileName,'w') as f:
            f.write('NAME\tH1r_4861A\tH1r_4861Ae\tH1r_6563A\tH1r_6563Ae\tN2_5755A\tN2_5755Ae\tN2_6548A\tN2_6548Ae\tN2_6584A\tN2_6584Ae\tO3_4363A\tO3_4363Ae\tO3_5007A\tO3_5007Ae\tS2_6716A\tS2_6716Ae\tS2_6731A\tS2_6731Ae\n')
            for id in idsPaper:
                print('id = ',type(id),': <'+id+'>')
                idx = -1
                for i in range(len(idsLines)):
                    if idsLines[i] == id:
                        idx = i
                print('idx = ',idx)
        #        print('np.where(',idsLines,' == <'+id+'>) = ',np.where(idsLines == id))
        #        idx = np.where(idsLines == id)[0]
                if idx == -1:
                    STOP
        #        print("csvLines.getData('$\mathrm{[NII]_{5755}}$') = ",csvLines.getData('$\mathrm{[NII]_{5755}}$'))
                NII5755 = float(csvLines.getData('$\mathrm{[NII]_{5755}}$',idx))
                NII6548 = float(csvLines.getData('$\mathrm{[NII]_{6548}}$',idx))
                NII6583 = float(csvLines.getData('$\mathrm{[NII]_{6583}}$',idx))
                print('NII5755 = ',NII5755)
                print('NII6548 = ',NII6548)
                print('NII6583 = ',NII6583)

                SII6716 = float(csvLines.getData('$\mathrm{[SII]_{6716}}$',idx))
                SII6731 = float(csvLines.getData('$\mathrm{[SII]_{6731}}$',idx))
                print('SII6716 = ',SII6716)
                print('SII6731 = ',SII6731)

                OIII4363 = float(csvLines.getData('$\mathrm{[OIII]_{4363}}$',idx))
                OIII5007 = float(csvLines.getData('$\mathrm{[OIII]_{5007}}$',idx))
                print('OIII4363 = ',OIII4363)
                print('OIII5007 = ',OIII5007)

                Halpha = float(csvLines.getData('$\mathrm{H_\\alpha}$',idx))
                Hbeta = float(csvLines.getData('$\mathrm{H_\\beta}$',idx))
                if Hbeta == 0.:
                    Hbeta = Halpha / 1000.
                print('Halpha = ',Halpha)
                print('Hbeta = ',Hbeta)
                if Hbeta > 0.:
                    print('Halpha / Hbeta = ',Halpha / Hbeta)

    #        f.write('NAME\tH1r_4861A\tH1r_4861Ae\tH1r_6563A\tH1r_6563Ae\tN2_5755A\tN2_5755Ae\tN2_6548A\tN2_6548Ae\tN2_6584A\tN2_6584Ae\tO3_4363A\tO3_4363Ae\tO3_5007A\tO3_5007Ae\tS2_6716A\tS2_6716Ae\tS2_6731A\tS2_6731Ae\n')
                f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (id,
                                                                    str(Hbeta),
                                                                    str(Hbeta*0.1),
                                                                    str(Halpha),
                                                                    str(Halpha*0.1),
                                                                    str(NII5755),
                                                                    str(NII5755*0.1),
                                                                    str(NII6548),
                                                                    str(NII6548*0.1),
                                                                    str(NII6583),
                                                                    str(NII6583*0.1),
                                                                    str(OIII4363),
                                                                    str(OIII4363*0.1),
                                                                    str(OIII5007),
                                                                    str(OIII5007*0.1),
                                                                    str(SII6716),
                                                                    str(SII6716*0.1),
                                                                    str(SII6731),
                                                                    str(SII6731*0.1),
                                                                    ))
    (_, _, filenames) = next(os.walk(spectraPath))
    for filename in filenames:
        if filename[filename.rfind('.'):] == '.fits':
            idPNMain = None
            if filename == 'Pa30_GT080716.fits':
                idPNMain = '15569'
            else:
                hashFitsFiles = csv.DictReader(open('/Users/azuri/daten/uni/HKU/IPHAS-GTC/fitsfiles.csv'))
                for hashFitsFile in hashFitsFiles:
                    if hashFitsFile['fileName'] == filename:
                        idPNMain = hashFitsFile['idPNMain']
            if idPNMain is None:
                print('could not find filename ',filename,' in hashFitsFiles')
                STOP
            print('filename = ',filename,': idPNMain = ',idPNMain)
            calculateErrors(os.path.join(spectraPath,filename),idPNMain,obsFileName)

#    for i in range(len(idsPaper)):
#        intensitiesObs = obs.getIntens(obsName=idsPaper[i],returnObs=True)
#        intensitiesCor = obs.getIntens(obsName=idsPaper[i],returnObs=False)
#        print('id = ',idsPaper[i],': E(B-V) = ',E_BV[i])
#        print('id = ',idsPaper[i],': H_alpha = ',intensitiesObs['H1r_6563A'],' => ',intensitiesCor['H1r_6563A'])
#        print('id = ',idsPaper[i],', H_beta = ',intensitiesObs['H1r_4861A'],' => ',intensitiesCor['H1r_4861A'])
#        print('id = ',idsPaper[i],', N2_5755A = ',intensitiesObs['N2_5755A'],' => ',intensitiesCor['N2_5755A'])
#        print('id = ',idsPaper[i],', N2_6548A = ',intensitiesObs['N2_6548A'],' => ',intensitiesCor['N2_6548A'])
#        print('id = ',idsPaper[i],', N2_6584A = ',intensitiesObs['N2_6584A'],' => ',intensitiesCor['N2_6584A'])
#        print('id = ',idsPaper[i],', O3_4363A = ',intensitiesObs['O3_4363A'],' => ',intensitiesCor['O3_4363A'])
#        print('id = ',idsPaper[i],', O3_5007A = ',intensitiesObs['O3_5007A'],' => ',intensitiesCor['O3_5007A'])
#        print('id = ',idsPaper[i],', S2_6716A = ',intensitiesObs['S2_6716A'],' => ',intensitiesCor['S2_6716A'])
#        print('id = ',idsPaper[i],', S2_6731A = ',intensitiesObs['S2_6731A'],' => ',intensitiesCor['S2_6731A'])
    #STOP

    if len(idsPaper) != len(idsLines):
        sortedIndicesPaper = np.argsort(idsPaper)
        sortedIndicesLines = np.argsort(idsLines)
        for i in range(np.min([len(idsPaper),len(idsLines)])):
            print('idsPaper[',sortedIndicesPaper[i],'] = ',idsPaper[sortedIndicesPaper[i]],', idsLines[',sortedIndicesLines[i],'] = ',idsLines[sortedIndicesLines[i]])
        if len(idsPaper) < len(idsLines):
            for i in np.arange(len(idsPaper),len(idsLines),1):
                print('idsLines[',sortedIndicesLines[i],'] = ',idsLines[sortedIndicesLines[i]])
        else:
            for i in np.arange(len(idsLines),len(idsPaper),1):
                print('idsPaper[',sortedIndicesPaper[i],'] = ',idsPaper[sortedIndicesPaper[i]])
        STOP
    print('len(idsPaper) = ',len(idsPaper))
    print('len(idsLines) = ',len(idsLines))
    #STOP
    for i in range(len(idsPaper)):
        id = idsPaper[i]
        print('id = ',type(id),': <'+id+'>')
        idx = -1
        for j in range(len(idsLines)):
            if idsLines[j] == id:
                idx = j
        print('idx = ',idx)
#        print('np.where(',idsLines,' == <'+id+'>) = ',np.where(idsLines == id))
#        idx = np.where(idsLines == id)[0]
        if idx == -1:
            STOP
#        print("csvLines.getData('$\mathrm{[NII]_{5755}}$') = ",csvLines.getData('$\mathrm{[NII]_{5755}}$'))

#        if math.isinf(E_BV[i]):
#            NII5755 = float(csvLines.getData('$\mathrm{[NII]_{5755}}$',idx))
#            NII6548 = float(csvLines.getData('$\mathrm{[NII]_{6548}}$',idx))
#            NII6583 = float(csvLines.getData('$\mathrm{[NII]_{6583}}$',idx))
#
#            SII6716 = float(csvLines.getData('$\mathrm{[SII]_{6716}}$',idx))
#            SII6731 = float(csvLines.getData('$\mathrm{[SII]_{6731}}$',idx))
#
#            OIII4363 = float(csvLines.getData('$\mathrm{[OIII]_{4363}}$',idx))
#            OIII5007 = float(csvLines.getData('$\mathrm{[OIII]_{5007}}$',idx))
#            print('NII5755 = ',NII5755)
#            print('NII6548 = ',NII6548)
#            print('NII6583 = ',NII6583)
#            print('SII6716 = ',SII6716)
#            print('SII6731 = ',SII6731)
#            print('OIII4363 = ',OIII4363)
#            print('OIII5007 = ',OIII5007)
#            print('obs.getIntens(obsName = ',id,', returnObs = True) = ',obs.getIntens(obsName = id, returnObs = True))
#            STOP
#            Halpha = float(csvLines.getData('$\mathrm{H_\\alpha}$',idx))
#            Hbeta = float(csvLines.getData('$\mathrm{H_\\beta}$',idx))
#        else:

        thisObsFileName = obsFileName[:obsFileName.rfind('.')]+'_'+id+obsFileName[obsFileName.rfind('.'):]
        print('thisObsFileName = ',thisObsFileName)
        obs = pn.Observation(thisObsFileName, fileFormat='lines_in_cols', errIsRelative=False, delimiter='\t')
        obs.addMonteCarloObs(N = 500)
        obs.extinction.law = 'F99'
        obs.def_EBV(label1='H1r_6563A',label2='H1r_4861A',r_theo=2.86)
        print('E(B-V) = ',obs.extinction.E_BV)
        E_BV = obs.extinction.E_BV
#        for i in range(len(idsPaper)):
#            intensitiesObs = obs.getIntens(obsName=idsPaper[i],returnObs=True)
#            print('id = ',idsPaper[i],': intensities = ',type(intensitiesObs),': ',intensitiesObs)
#            print('id = ',idsPaper[i],': H_alpha = ',intensitiesObs['H1r_6563A'],
#                                    ', H_beta = ',intensitiesObs['H1r_4861A'],
#                                    ': H_a/H_b = ',intensitiesObs['H1r_6563A']/intensitiesObs['H1r_4861A'] if intensitiesObs['H1r_4861A'] > 0. else 0.,
#                                    ', E_BV = ',E_BV[i])
        E_BV[np.where(E_BV < 0.)[0]] = 0.
        print('len(idsPaper) = ',len(idsPaper),', len(E_BV) = ',len(E_BV),', len(idsLines) = ',len(idsLines))
        print('E(B-V) = ',obs.extinction.E_BV)
        obs.correctData()

        intensities = obs.getIntens(returnObs=False)#,obsName = id,
        """H1_4861A\tH1_6563A\tN2_5755A\tN2_6548A\tN2_6584A\tO3_4363A\tO3_5007A\tS2_6716A\tS2_6731A"""
        Halpha = intensities['H1r_6563A']
        Hbeta = intensities['H1r_4861A']
        NII5755 = intensities['N2_5755A']
        NII6548 = intensities['N2_6548A']
        NII6583 = intensities['N2_6584A']
        SII6716 = intensities['S2_6716A']
        SII6731 = intensities['S2_6731A']
        OIII4363 = intensities['O3_4363A']
        OIII5007 = intensities['O3_5007A']
        print('obs.getIntens(obsName = ',id,', returnObs = True) = ',obs.getIntens(returnObs = True))#,obsName = id,
        print('NII5755 = ',NII5755.shape,': ',NII5755)
        print('NII6548 = ',NII6548.shape,': ',NII6548)
        print('NII6583 = ',NII6583.shape,': ',NII6583)
        print('SII6716 = ',SII6716.shape,': ',SII6716)
        print('SII6731 = ',SII6731.shape,': ',SII6731)
        print('OIII4363 = ',OIII4363.shape,': ',OIII4363)
        print('OIII5007 = ',OIII5007.shape,': ',OIII5007)
        #STOP
#        rc = pn.RedCorr(E_BV = 1.2, R_V = 3.2, law = 'F99')
#        rc.setCorr(Halpha / Hbeta, 6563., 4861.)
#        wave = [linesOfInterest['SIIa'],
#                linesOfInterest['SIIb'],
#                linesOfInterest['NII'],
#                linesOfInterest['NIIa'],
#                linesOfInterest['NIIb'],
#                linesOfInterest['OIIIa'],
#                linesOfInterest['OIIIb'],
#                linesOfInterest['Halpha'],
#                linesOfInterest['Hbeta'],]
#        correc = rc.getCorr(wave)
#        print('wave = ',wave,': correc = ',correc)
#        STOP

        temNS,denNS,eTemNS,eDenNS = [None,None,None,None]
        temNSa,denNSa,eTemNSa,eDenNSa = [None,None,None,None]
        temNSb,denNSb,eTemNSb,eDenNSb = [None,None,None,None]
        temOS,denOS,eTemOS,eDenOS = [None,None,None,None]
#        denS = None
#        denO = None
#        denN = None
        NIIHa,eNIIHa = [None,None]
        OIIIHb,eOIIIHb = [None,None]
        SIIHa,eSIIHa = [None,None]
        HaHb,eHaHb = [None,None]
        if (NII5755[0] > 0.) and (NII6548[0] > 0.) and (NII6583[0] > 0.):
            if Halpha[0] > 0.:
                arr = removeBadValues(NII6583 / Halpha)
                NIIHa, eNIIHa = [np.mean(arr), np.std(arr)]
#            denN = n2.getTemDen((NII6583 + NII6548) / NII5755, tem=10000., to_eval = '(L(6584) + L(6548)) / L(5755)')
#            print('idPNMain = ',id,': denN = ',denN)
#
#             STOP
        if (OIII4363[0] > 0.) and (OIII5007[0] > 0.):
            if Hbeta[0] > 0.:
                arr = removeBadValues(OIII5007 / Hbeta)
                OIIIHb,eOIIIHb = [np.mean(arr),np.std(arr)]
#            denO = o3.getTemDen(OIII5007 / OIII4363, tem=10000., wave1=5007, wave2=4363, maxIter=1000)
#            print('idPNMain = ',id,': denO = ',denO)
        if  (SII6716[0] > 0.) and (SII6731[0] > 0.):
            if Halpha[0] > 0.:
                arr = removeBadValues((SII6716 + SII6731) / Halpha)
                SIIHa,eSIIHa = [np.mean(arr),np.std(arr)]
#            denS = s2.getTemDen(int_ratio=SII6716/SII6731,tem=10000.,wave1=6716,wave2=6731,maxIter=1000)
#            print('idPNMain = ',id,': denS = ',denS)
            if (NII5755[0] > 0.) and ((NII6548[0] > 0.) or (NII6583[0] > 0.)):
                #if NII6583 > NII6548:
                try:
                    print('trying getCrossTemDen(NII,SII)')
                    tem, den = diags.getCrossTemDen(diag_tem='[NII] 5755/6584',
                                                        diag_den='[SII] 6731/6716',
                                                        obs=obs,
                                                        #value_tem=NII5755 / NII6548,
                                                        #value_den=SII6731 / SII6716,
                                                        )
                    temNS, eTemNS = getMeanAndStd(tem)
                    denNS, eDenNS = getMeanAndStd(den)
                    print('idPNMain = ',id,': NII6583 > NII6548: temNS = ',temNS,', denNS = ',denNS)
                except:
                    print('getCrossTemDen 5755/6584 failed')
                #else:
                try:
                    print('trying getCrossTemDen(NII,SII)a')
                    tem, den = diags.getCrossTemDen(diag_tem='[NII] 5755/6548',
                                                        diag_den='[SII] 6731/6716',
                                                        obs=obs,
                                                        #value_tem=NII5755 / NII6548,
                                                        #value_den=SII6731 / SII6716,
                                                        )
                    temNSa, eTemNSa = getMeanAndStd(tem)
                    denNSa, eDenNSa = getMeanAndStd(den)
                    print('idPNMain = ',id,': NII6583 < NII6548: temNSa = ',temNSa,', denNSa = ',denNSa)
                except:
                    print('getCrossTemDen 5755/6548 failed')
                try:
                    print('trying getCrossTemDen(NII,SII)b')
                    tem, den = diags.getCrossTemDen(diag_tem='[NII] 5755/6584+',
                                                        diag_den='[SII] 6731/6716',
                                                        obs=obs,
                                                        #value_tem=NII5755 / NII6548,
                                                        #value_den=SII6731 / SII6716,
                                                        )
                    temNSb, eTemNSb = getMeanAndStd(tem)
                    denNSb, eDenNSb = getMeanAndStd(den)
                    print('idPNMain = ',id,': temNSb = ',temNSb,', denNS = ',denNSb)
                except:
                    print('getCrossTemDen 5755/6584+ failed')
                #STOP
            if (OIII4363[0] > 0.) and (OIII5007[0] > 0.):
                tem, den = diags.getCrossTemDen(diag_tem='[OIII] 4363/5007',
                                                    diag_den='[SII] 6731/6716',
                                                    obs=obs,
#                                                    #value_tem=OIII4363 / OIII5007,
#                                                    #value_den=SII6731 / SII6716,
                                                    )
                temOS,eTemOS = getMeanAndStd(tem)
                denOS,eDenOS = getMeanAndStd(den)
                print('idPNMain = ',id,': temOS = ',temOS,', denOS = ',denOS)
#                STOP
        if (Halpha[0] > 0.) and (Hbeta[0] > 0.):
            arr = removeBadValues(Halpha / Hbeta)
            HaHb, eHaHb = [np.mean(arr), np.std(arr)]
        print('csvOut.header = ',csvOut.header)
        E_BV = removeBadValues(E_BV)
        if math.isnan(np.mean(E_BV)):
            print('E_BV = ',E_BV)
            csvOut.setData('E(B-V)',idx,'0.0')
        else:
            csvOut.setData('E(B-V)',idx,'$%.2f \pm %.2f$' % (np.mean(E_BV), np.std(E_BV)))
        csvOut.setData('Target Name',idx,getNameFromIdPNMain(id,csvPaper))
        csvOut.setData('$\mathrm{(H_\\alpha)/(H_\\beta)}$',idx,'$%.2f \pm %.2f$' % (HaHb, eHaHb) if ((HaHb is not None) and (not math.isnan(HaHb))) else ' ')
        csvOut.setData('$\mathrm{(6717,31)/(H_\\alpha)}$',idx,'$%.2f \pm %.2f$' % (SIIHa,eSIIHa) if ((SIIHa is not None) and (not math.isnan(SIIHa))) else ' ')
        csvOut.setData('$\mathrm{(5007)/(H_\\beta)}$',idx,'$%.2f \pm %.2f$' % (OIIIHb,eOIIIHb) if ((OIIIHb is not None) and (not math.isnan(OIIIHb))) else ' ')
        csvOut.setData('$\mathrm{(6584)/(H_\\alpha)}$',idx,'$%.2f \pm %.2f$' % (NIIHa,eNIIHa) if ((NIIHa is not None) and (not math.isnan(NIIHa))) else ' ')
        csvOut.setData('$\mathrm{T_{e^-}([NII],[SII])}$',idx,'$%d \pm %d$' % (int(np.mean([temNS,temNSa,temNSb])),int(np.mean([eTemNS,eTemNSa,eTemNSb]))) if ((temNS is not None) and (not math.isnan(temNS))) else ' ')
#        csvOut.setData('$\mathrm{T_{e^-}([NII],[SII])a}$',idx,'$%d \pm %d$' % (int(temNSa),int(eTemNSa)) if ((temNSa is not None) and (not math.isnan(temNS))) else ' ')
#        csvOut.setData('$\mathrm{T_{e^-}([NII],[SII])b}$',idx,'$%d \pm %d$' % (int(temNSb),int(eTemNSb)) if ((temNSb is not None) and (not math.isnan(temNS))) else ' ')
        csvOut.setData('$\mathrm{T_{e^-}([OIII],[SII])}$',idx,'$%d \pm %d$' % (int(temOS),int(eTemOS)) if ((temOS is not None) and (not math.isnan(temOS))) else ' ')
        csvOut.setData('$\mathrm{\\rho_{e^-}([NII],[SII])}$',idx,'$%d \pm %d$' % (int(np.mean([denNS,denNSa,denNSb])),int(np.mean([eDenNS,eDenNSa,eDenNSb]))) if ((denNS is not None) and (not math.isnan(denNS))) else ' ')
#        csvOut.setData('$\mathrm{\\rho_{e^-}([NII],[SII])a}$',idx,'$%d \pm %d$' % (int(denNSa),int(eDenNSa)) if ((denNSa is not None) and (not math.isnan(denNS))) else ' ')
#        csvOut.setData('$\mathrm{\\rho_{e^-}([NII],[SII])b}$',idx,'$%d \pm %d$' % (int(denNSb),int(eDenNSb)) if ((denNSb is not None) and (not math.isnan(denNS))) else ' ')
        csvOut.setData('$\mathrm{\\rho_{e^-}([OIII],[SII])}$',idx,'$%d \pm %d$' % (int(denOS),int(eDenOS)) if ((denOS is not None) and (not math.isnan(denOS))) else ' ')

    sortedIndices = np.argsort(np.array([name.lower() for name in csvOut.getData('Target Name')]))
    csvOut.sort(sortedIndices)
    csvFree.writeCSVFile(csvOut,os.path.join(imPath[:imPath.rfind('/')],'lines_and_temden.csv'))
    csvFree.writeCSVFile(csvOut,os.path.join(imPath[:imPath.rfind('/')],'table_paper_sorted_with_temden.csv'),'&')

#def plots():
#    csvOut = csvFree.readCSVFile(os.path.join(imPath[:imPath.rfind('/')],'lines_and_temden.csv'))

    NIIHa = csvOut.getData('$\mathrm{(6584)/(H_\\alpha)}$')
    SIIHa = csvOut.getData('$\mathrm{(6717,31)/(H_\\alpha)}$')
    OIIIHb = csvOut.getData('$\mathrm{(5007)/(H_\\beta)}$')
    n = []
    s = []
    o = []
    for i in range(len(NIIHa)):
        print('NIIHa[',i,'] = <'+NIIHa[i]+'>, SIIHa[',i,'] = <'+SIIHa[i]+'>, OIIIHb[',i,'] = <'+OIIIHb[i]+'>')
        if (NIIHa[i].replace(' ','') != '') and (SIIHa[i].replace(' ','') != '') and (OIIIHb[i].replace(' ','') != ''):
            n.append(float(NIIHa[i][1:NIIHa[i].find(' ')]))
            s.append(float(SIIHa[i][1:SIIHa[i].find(' ')]))
            o.append(float(OIIIHb[i][1:OIIIHb[i].find(' ')]))
    plt.scatter(np.log10(np.array(n)),np.log10(np.array(o)),color='k',s=10,marker='D')
    plt.xlabel(r'log $\mathrm{(6584)/H_\alpha}$')
    plt.ylabel(r'log $\mathrm{(5007)/H_\beta}$')
    plt.savefig(os.path.join(latexPath,'images/OIIIHbetaVsNIIHalpha.pdf'),bbox_inches='tight')
    plt.show()
    plt.close()

    plt.scatter(np.log10(np.array(s)),np.log10(np.array(o)),color='k',s=10,marker='D')
    plt.xlabel(r'log $\mathrm{(6716,31)/H_\alpha}$')
    plt.ylabel(r'log $\mathrm{(5007)/H_\beta}$')
    plt.savefig(os.path.join(latexPath,'images/OIIIHbetaVsSIIHalpha.pdf'),bbox_inches='tight')
    plt.show()
    plt.close()
    return csvOut


def getNameFromIdPNMain(idPNMain,csvTablePaper):
#    print('csvTablePaper.header = ',csvTablePaper.header)
#    print('looking for idPNMain = ',type(idPNMain),': <'+idPNMain+'>')
    for i in range(csvTablePaper.size()):
#        print("csvTablePaper.getData(' HASH ID ',",i,") = <"+csvTablePaper.getData(' HASH ID ',i)+'>')
        if csvTablePaper.getData(' HASH ID ',i).replace(' ','') == idPNMain:
            print('found idPNMain at position ',i)
            return csvTablePaper.getData('Name ',i)
    print('did not find idPNMain = ',idPNMain)
    STOP

def writeFinalTable():
    lines = []

    csvTable = csvFree.readCSVFile(os.path.join(imPath[:imPath.rfind('/')],'table_paper_sorted_with_temden.csv'),'&',False)
    """remove line ratios"""
    csvTable.removeColumn('$\mathrm{(H_\\alpha)/(H_\\beta)}$')
    csvTable.removeColumn('$\mathrm{(6717,31)/(H_\\alpha)}$')
    csvTable.removeColumn('$\mathrm{(5007)/(H_\\beta)}$')
    csvTable.removeColumn('$\mathrm{(6584)/(H_\\alpha)}$')
    csvTable.removeColumn('$\mathrm{H_\\alpha}$')
    csvTable.removeColumn('$\mathrm{H_\\beta}$')
    csvTable.removeColumn('$\mathrm{[SII]_{6716}}$')
    csvTable.removeColumn('$\mathrm{[SII]_{6731}}$')
    csvTable.removeColumn('$\mathrm{[NII]_{5755}}$')
    csvTable.removeColumn('$\mathrm{[NII]_{6548}}$')
    csvTable.removeColumn('$\mathrm{[NII]_{6583}}$')
    csvTable.removeColumn('$\mathrm{[OIII]_{4363}}$')
    csvTable.removeColumn('$\mathrm{[OIII]_{5007}}$')
    csvFree.writeCSVFile(csvTable,os.path.join(imPath[:imPath.rfind('/')],'table_paper_sorted_with_temden_minus_lineRatios.csv'),'&')

    csvVRad = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/vrad.csv')

    with open(os.path.join(imPath[:imPath.rfind('/')],'table_paper_sorted_with_temden_minus_lineRatios.csv'),'r') as f:
        lines = f.readlines()
    with open(os.path.join(latexPath,'table_paper_sorted_with_temden.tex'),'w') as f:
        f.write('\\documentclass[8pt]{article}\n')
        f.write('\\usepackage[margin=2cm,a4paper]{geometry}')
        f.write('\\usepackage{graphicx}\n')
        f.write('\\usepackage{natbib}\n')
        f.write('\\usepackage{longtable}\n')
#        f.write('\\usepackage{pdflscape}\n')
        f.write('\\setcitestyle{numbers}\n')
        f.write('\\setcitestyle{square}\n')
        f.write('\\begin{document}\n')
#        f.write('\\begin{landscape}\n')
        f.write('\\clearpage\n')
        f.write('\\onecolumn\n')
        f.write('\\begin{longtable}{ | *{21}{l|} }\n')
        f.write('\\hline\n')
        f.write(lines[0].rstrip()+'\\\\\n')
        f.write('\\endhead  % header material\n')
        f.write('\\hline\endfoot  % footer material\n')
        f.write('\\hline\n')
        for line in lines[1:]:
            if line.count('pm') > 1:
                f.write(line.rstrip()+'\\\\\n')
        f.write('\\hline\n')
        f.write('\\end{longtable}\n')
#        f.write('\\end{landscape}\n')
        f.write('\\clearpage\n')
        f.write('\\twocolumn\n')
        f.write('\\end{document}\n')

def getExptimes():
    import re
    textFilesList = '/Volumes/discovery/spectra/IPHAS_GTC_DATA/textfiles.list'
    with open(textFilesList,'r') as f:
        filenames = f.readlines()
    print('filenames = ',filenames)
    csv = csvData.CSVData()
    for filename in filenames:
        with open(filename.strip(),'r') as f:
            lines = f.readlines()
        for iLine in range(len(lines)):
            lineData = lines[iLine].strip()
            if lineData[:2] == 'OB':
                print('lineData = ',lineData)
                head = re.sub(' +', ' ', lineData).split(' ')
                print('head = ',len(head),': ',head)
                if csv.header == []:
                    csv.header = head
                dat = re.sub(' +',' ',re.sub('\t',' ',lines[iLine+2].strip())).split(' ')
                print('dat = ',len(dat),': ',dat)
                csv.append(dat)
    print('csv.header = ',csv.header)
    print('csv.data = ',csv.data)

    table1 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_PNe_sorted_withPAandExpTime.tex','&',False)
    print('table1.size() = ',table1.size())
    table1.append(csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_oldPNe_sorted_withPAandExpTime.tex','&',False).data)
    print('table1.size() = ',table1.size())
    print('table1.header = ',table1.header)
    ras = table1.getData(' RA ')
    decs = table1.getData(' DEC ')
    print('ras = ',len(ras),': ',ras)
    for i in range(table1.size()):
        table1.setData('$\\mathrm{t_{exp}}$ [s]',i,'')
    #print('expTimes = ',table1.getData('$\\mathrm{t_{exp}}$ [s]'))
    #STOP

    nFound = 0
    nNotFound = 0
    for i in range(csv.size()):
        ra = csv.getData('RA',i)
        dec = csv.getData('DEC',i)
        nexp = csv.getData('NEXP',i)
        exptime = csv.getData('EXPTIME',i)
        print('RA = '+ra+', DEC = '+dec+': '+nexp+' * '+exptime)
        found = False
        minSep = 1000.
        idxMinSep = -1
        for iRA in range(len(ras)):
            sep = angularDistancePyAsl(hmsToDeg(ra),dmsToDeg(dec),hmsToDeg(ras[iRA]),dmsToDeg(decs[iRA]))*3600.
            if sep < minSep:
                minSep = sep
                idxMinSep = iRA
        if minSep < 10.:
            print('found object in ',idxMinSep)
            found = True
            nFound += 1
            if table1.getData('$\\mathrm{t_{exp}}$ [s]',idxMinSep) != '':
                table1.setData('$\\mathrm{t_{exp}}$ [s]',idxMinSep,table1.getData('$\\mathrm{t_{exp}}$ [s]',idxMinSep)+' + '+(nexp+'*' if nexp != '1' else '')+exptime)
            else:
                table1.setData('$\\mathrm{t_{exp}}$ [s]',idxMinSep,(nexp+'*' if nexp != '1' else '')+exptime)
        if not found:
            print("ERROR: could not find RA ",ra,' in table1')
            nNotFound += 1
    print('nFound = ',nFound,', nNotFound = ',nNotFound)
    print('exptimes = ',table1.getData('$\\mathrm{t_{exp}}$ [s]'))
    for i in range(table1.size()):
        table1.setData('$\\mathrm{t_{exp}}$ [s]',i,table1.getData('$\\mathrm{t_{exp}}$ [s]',i)+'\\\\')
    csvFree.writeCSVFile(table1,'/Users/azuri/daten/uni/HKU/IPHAS-GTC/table1_allPNe_sorted_withPAandExpTime.tex','&')

if __name__ == '__main__':
    print('reading table_paper_sorted')
    csvPaper = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/table_paper_sorted.csv','&',False)
    print(csvPaper.header)
#    print(csvPaper.data)

    csvTargets = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/IPHAS-GTC/Charts_IPHAS_Amateurs/Targets_GTC.cvs')
    print(csvTargets.header)
#    print(csvTargets.data)
    #    coneSearch(csvPaper)

    #fixInUseInIquote()
    if False:
        ids = findHASHid(csvPaper, csvTargets)
#        ids.append(['Ou 1','8458'])
#        ids.append(['IPHASX J055242.8+262116','9824'])
#        ids.append(['IPHASXJ190333','8506'])
#        ids.append(['BMP J1917+0200','2502'])
        print(ids)
        #STOP
        getImages(ids)
        csvPaper = csvFree.readCSVFile(os.path.join(latexPath,'table_true_names_sorted.tex'),'&',False)
#        filterIphasImagesInBashFile()
#        filterSSSImagesInBashFile()
#        filterWISE432ImagesInBashFile()
#        filterNVSSImagesInBashFile()
#        filterGalexImagesInBashFile()
        csvLines = makeSpectraTable(ids, hashPNMainFileName='/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_PNMain.csv', calculateLineIntensities = True)
#        for id in ids:
#            print('id = ',id)
#        STOP
#        for id in ids:
#            combineImages(id[1])
#        createImageTable(ids)
#        splitObservationFile(os.path.join(imPath[:imPath.rfind('/')],'observation.dat'))
#        addLinesToTable(csvPaper,obsFileName = os.path.join(imPath[:imPath.rfind('/')],'observation.dat'))#,csvLines)

#    plots()

#    writeFinalTable()
    getExptimes()
