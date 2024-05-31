import os
from astroplan import Observer
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u
from distutils.dir_util import copy_tree
import numpy as np
from PIL import Image
import shutil
import subprocess


from plot_obs_planning import plot_target
import csvFree,csvData
from myUtils import angularDistancePyAsl,hmsToDeg,dmsToDeg
from drUtils import createFindingChartFromFits

#from astroplan import download_IERS_A
#download_IERS_A()

#goodTargetsDir = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good'
goodTargetsDir = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06'
outDir = os.path.join(goodTargetsDir,'targets_SAAO_2024-06-05_priority_good/findingCharts')
fileList = os.path.join(goodTargetsDir,'allFiles_priority.list')

idPNMain_hash_no_spectrum = csvFree.readCSVFile(os.path.join(goodTargetsDir,'hash_TLPc_no_spectrum_210524.csv')).getData('idPNMain')
#idPNMain_literature_available = csvFree.readCSVFile(os.path.join(goodTargetsDir,'literature_spectrum_available.csv')).getData('idPNMain')
#idPNMain_elcat_available = idPNMain_literature_available#csvFree.readCSVFile(os.path.join(goodTargetsDir,'elcat_available.csv')).getData('idPNMain')

pnMain = csvFree.readCSVFile(os.path.join(goodTargetsDir,'hash_PNMain_210524.csv'))#'hash_TLPc_no_spectrum_010623.csv'))

fitsFiles = csvFree.readCSVFile(os.path.join(goodTargetsDir,'hash_FitsFiles_210524.csv'))

usrCommentsFile = os.path.join(goodTargetsDir,'hash_UsrComm_210524.csv')

angDiams = csvFree.readCSVFile(os.path.join(goodTargetsDir,'hash_AngDiam_210524.csv'))

names = csvFree.readCSVFile(os.path.join(goodTargetsDir,'hash_Names_210524.csv'))

allPossibleTargets = os.path.join(goodTargetsDir,'pnMain_need_spectrum.csv')

#for i in np.arange(pnMain.size()-1,-1,-1):
#    if not pnMain.getData('idPNMain',i) in idPNMain_hash_no_spectrum:
#        pnMain.removeRow(i)
#csvFree.writeCSVFile(pnMain,allPossibleTargets)
#
allPNe = csvFree.readCSVFile(allPossibleTargets)
#usrComments = csvFree.readCSVFile(usrCommentsFile)
#elcats = ['elcat_dopita.csv','elcat_ercolano.csv','elcat_frew.csv','elcat_hsia.csv','elcat_kerber.csv','elcat_kraan.csv','elcat_kwitter.csv']

observatoryName = "SAAO"
observatoryLocation = EarthLocation(lat=-32.3783*u.deg, lon=20.8105*u.deg, height=1750*u.m)
utcoffset = 2*u.hour
date = '2024-06-05'
midnight = Time(date+' 00:00:00') - utcoffset
minAltitude = 30.0

def getIdPNMainFromFileName(fName):
#    print('fName = '+fName)
    idPNTemp = fName[fName.rfind('/')+1:]
#    print('idPNTemp = <'+idPNTemp+'>')
    idPNTemp = idPNTemp[idPNTemp.find('_')+1:]
#    print('idPNTemp = <'+idPNTemp+'>')
    idPNTemp = idPNTemp[idPNTemp.find('_')+1:]
#    print('idPNTemp = <'+idPNTemp+'>')
    idPNTemp = idPNTemp[:idPNTemp.find('_')]
#    print('idPN = <'+idPNTemp+'>')
#    STOP
    return idPNTemp

def getAllIdPNMains(fNamesList):
    idPNMains = []
    for fName in fNamesList:
        id = getIdPNMainFromFileName(fName)
        if id not in idPNMains:
            idPNMains.append(id)
    return idPNMains

def getLinesForIdPNMain(lines,idPNMain):
    outLines = []
    for line in lines:
        idPNTemp = getIdPNMainFromFileName(line)
        if idPNTemp == idPNMain:
            outLines.append(line)
    return outLines

def getImagesForIdPNMain(lines,idPNMain):
    tempLines = getLinesForIdPNMain(lines,idPNMain)
    outLines = []
    for line in tempLines:
        if 'centroid' not in line:
            outLines.append(line)
    return outLines

def getObservingTimeStartAndEnd():
#    midnight = Time(date+' 00:00:00') - utcoffset
    obs = Observer.at_site(observatoryName)#, timezone='Eastern Standard Time')
    observationStartTime = obs.twilight_evening_astronomical(midnight)#Time('2019-9-5 19:10:00') - utcoffset
    observationEndTime = obs.twilight_morning_astronomical(midnight)#Time('2019-9-6 4:55:00') - utcoffset
    observationStartTime.format = 'iso'
    observationEndTime.format = 'iso'
    return [observationStartTime+utcoffset, observationEndTime+utcoffset]

def makeFindingCharts():
    with open(fileList,'r') as f:
        lines = f.readlines()
    fNamesList = [i.rstrip('\n') for i in lines]

    csvAllTargets = csvFree.readCSVFile(allPossibleTargets)

    obsStartTime, obsEndTime = getObservingTimeStartAndEnd()
    obsStartTimeHour = obsStartTime.datetime.hour
    obsEndTimeHour = obsEndTime.datetime.hour
    print('obsStartTimeHour = ',obsStartTimeHour)
    print('obsEndTimeHour = ',obsEndTimeHour)
    for hour in np.arange(obsStartTimeHour,24,1):
        tmpDir = os.path.join(outDir,
                              str(hour)
                              +':00-'
                              +(str(hour+1) if hour < 23 else '00')
                              +':00')
        if not os.path.exists(tmpDir):
            os.mkdir(tmpDir)
        print('created dir <'+tmpDir+'>')
    for hour in np.arange(0,obsEndTimeHour+1,1):
        tmpDir = os.path.join(outDir,
                              '0'
                              +str(hour)
                              +':00-0'
                              +str(hour+1)
                              +':00')
        if not os.path.exists(tmpDir):
            os.mkdir(tmpDir)
        print('created dir <'+tmpDir+'>')

    idPNMains = getAllIdPNMains(fNamesList)
    nFailed = 0
    failedidPNMains = []
    for idPNMain in idPNMains:
        print('idPNMain = ',idPNMain)
        images = getImagesForIdPNMain(fNamesList,idPNMain)
        print('images = ',images)
#        for im in images:
#            if 'findingChart' in im:
#                STOP

        pos = csvAllTargets.find('idPNMain',idPNMain,0)[0]
        print('pos = ',pos)
        print("csvAllTargets.getData('MajDiam',pos) = ",csvAllTargets.getData('MajDiam',pos))
        try:
            majDiam = float(csvAllTargets.getData('MajDiam',pos))
        except:
            majDiam = 0.
        name = csvAllTargets.getData('Name',pos)
        if (majDiam < 200) and ('Ritter' not in name) and ('Objet' not in name):
            print('outDir = <'+outDir+'>')
            outDirName = os.path.join(outDir,idPNMain)
            print('outDirName = <'+outDirName+'>')
            print('creating outDirName = <'+outDirName+'>')
            if not os.path.exists(outDirName):
                os.mkdir(outDirName)
            if not os.path.exists(outDirName):
                STOP
            for im in images:
                if not 'findingChart' in im:
                    try:
                        isGood = False
                        if os.path.exists(outDirName):
                            print('working on im = ',im)
                            if im.endswith('.png'):
                                isGood = True
                                backgroundImageName = os.path.join(outDir[:outDir.rfind('/')],im)
                                background = Image.open(backgroundImageName)
                                print('backgroundImageName <'+backgroundImageName+'> opened')
                                overLayImageName = os.path.join(outDir[:outDir.rfind('/')],im[:-4])+'_centroid.png'
                                print('opening overlayImageName <'+overLayImageName+'>')
                                overlay = Image.open(overLayImageName)

                                background = background.convert("RGBA")
                                overlay = overlay.convert("RGBA")

                                new_img = Image.new("RGBA", background.size)
                                new_img = Image.alpha_composite(new_img, background)
                                new_img = Image.alpha_composite(new_img, overlay)
                                newImName = os.path.join(outDirName,im[im.rfind('/')+1:-4]+'_findingChart.png')
#                                newImName = os.path.join(outDirName,im[:-4]+'_findingChart.png')#.replace('.0','')
                                #print('newImName = <'+newImName+'>')
                                new_img.save(newImName,"PNG")
                                #print('outDirName = <'+outDirName+'>')
                                print('saved image <'+newImName+'>')
                                #STOP
                            elif im.endswith('.fits'):
                                if im.rfind('_shs') > 0:
                                    isGood = True
                                    tmp = im[:im.rfind('_shs')]
                                    idPNMain = tmp[tmp.rfind('_')+1:]

                                    name = csvAllTargets.getData('Name',pos)
                                    print('name = <'+name+'>: MajDiam = ',majDiam)
                                    ra = csvAllTargets.getData('RAJ2000',pos)
                                    dec = csvAllTargets.getData('DECJ2000',pos)
    #                                try:
                                    new_img_name = createFindingChartFromFits(im,#os.path.join(outDir[:outDir.rfind('/')],im),
                                                                240,
                                                                majDiam,
                                                                name,
                                                                ra,
                                                                dec,
                                                                outDirName,
                                                                display=False)
                                    print('saved new_img_name = <'+new_img_name+'>')
                                    #STOP
    #                                except Exception as e:
    #                                    print(e)
    #                                    STOP
#                                    shutil.move(outDirName+'*.png',finalDir+'/')
                                elif im.rfind('_Ha') > 0:
#                                    tmp = im[:im.rfind('_Ha')]
#                                    if not isGood:
                                    print('isGood = False')
#                                        STOP
#                                        shutil.rmtree(outDirName)
                                    #raise ValueError('IPHAS image set')
                                else:
                                    print('could not identify image type')
                                    STOP
                            else:
                                raise ValueError('could not identify file type of <'+im+'>')
                            if isGood:
                                visName = im[im.rfind('/')+1:im.rfind('_')]
#                                visName = visName[:visName.rfind('_')]
                                visName = os.path.join(os.path.join(outDir,idPNMain),visName+'_visibility.png')#.replace('.0','')
                                print('creating ',visName)
                                #STOP
                                targetCoord = SkyCoord(ra=float(csvAllTargets.getData('DRAJ2000',pos))*u.deg,
                                                    dec=float(csvAllTargets.getData('DDECJ2000',pos))*u.deg,
                                                    frame='icrs')
                                maxAltitudeTime = plot_target(targetCoord,
                                                            observatoryLocation,
                                                            observatoryName,
                                                            utcoffset,
                                                            date,
                                                            False,
                                                            outFileName=visName)
                                maxAltitudeHour = maxAltitudeTime.hour
                                print('maxAltitudeHour = ',maxAltitudeHour)
                                if maxAltitudeHour < 12:
                                    if maxAltitudeHour > obsEndTimeHour:
                                        maxAltitudeHour = obsEndTimeHour
                                else:
                                    if maxAltitudeHour < obsStartTimeHour:
                                        maxAltitudeHour = obsStartTimeHour
                                finalDir = os.path.join(outDir,
                                                        (str(maxAltitudeHour) if maxAltitudeHour > 9 else ('0'+str(maxAltitudeHour)))
                                                        +':00-'
                                                        +(str(maxAltitudeHour+1) if (maxAltitudeHour+1) > 9 else ('0'+str(maxAltitudeHour+1)))
                                                        +':00')
                                print('moving outDirName = <'+outDirName+'> to finalDir = <'+finalDir+'>')
    #                                    STOP
                                #if os.path.exists(finalDir):
                                #    shutil.rmtree(finalDir)
                                #shutil.move(outDirName,finalDir)
                                command = 'mv '+outDirName+'/* '+finalDir+'/'
                                print('command = <'+command+'>')
                                os.system(command)
                                #STOP
                            else:
                                print('isGood = False')
                                #STOP
                    #        print(csvAllTargets.getData(pos))

                    except Exception as e:
                        print(e)
                        shutil.rmtree(outDirName)
                        nFailed += 1
                        failedidPNMains.append(idPNMain)
                        #STOP
                        pass

            command = 'rmdir '+outDirName
            os.system(command)
    print(nFailed,' objects failed to plot visibility plot: ',failedidPNMains)
#    STOP

def getDarkTimesUT():
    time = midnight - 8*u.hour
    times = [time]
    while time < midnight + 8*u.hour:
        time += 1*u.minute
        times.append(time)
    utTimes = Time([timeX.value for timeX in times])
    frames = AltAz(obstime=utTimes, location=observatoryLocation)
    sunaltazs = get_sun(utTimes).transform_to(frames)
    utTimes = utTimes[np.where(sunaltazs.alt < -0*u.deg)]
    return utTimes

def getTargetAltitudes(target):
    utTimes = getDarkTimesUT()
    frames = AltAz(obstime=utTimes, location=observatoryLocation)
    targetaltazs = target.transform_to(frames)
    altitudes = targetaltazs.alt.to_value()
    return altitudes

def createVisibilityDirs(path,darkTimesLocal):
    times = []
    outDirs = []
    localTimes = darkTimesLocal
    for time in localTimes:
        time = time.strftime('%H:%M')
        minute = time[time.find(':')+1:]
        if (minute == '00') or (minute == '30'):
            times.append(time)
            outDirName = os.path.join(path,'visible_at_'+time)
            outDirs.append(outDirName)
            if not os.path.exists(outDirName):
                os.mkdir(outDirName)
    return [times,outDirs]

def getIDsFromDirList(inputList):
    print('getIDsFromDirList: inputList = <'+inputList+'>')
    with open(inputList,'r') as f:
        lines = f.readlines()
    lines = [line.rstrip('\n') for line in lines]
    print('lines = ',lines)

    dirs = []
    idPNs = []
    newDir = True
    ids = []
    for line in lines:
        if line == '':
            newDir = True
            idPNs.append(ids)
        elif newDir:
            dirs.append(line[:-1])
            ids = []
            newDir = False
        else:
            ids.append(line)
    idPNs.append(ids)
    return [dirs,idPNs]

def getIDsFromFindingChartsList(inputList):
    print('getIDsFromFindingChartsList: inputList = ',inputList)
    with open(inputList,'r') as f:
        lines = f.readlines()
    lines = [line.rstrip('\n') for line in lines]
    print('lines = ',lines)

    dirs = []
    idPNs = []
    newDir = True
    for line in lines:
        if line == '':
            newDir = True
            idPNs.append(ids)
        elif newDir:
            dirs.append(line[:-1])
            ids = []
            newDir = False
        else:
            idPNMain = line[line.find('_')+1:]
            idPNMain = idPNMain[idPNMain.find('_')+1:]
            idPNMain = idPNMain[:idPNMain.find('_')]
            print('idPNMain = ',idPNMain)
            if idPNMain not in ids:
                ids.append(idPNMain)
    idPNs.append(ids)
    print('dirs = ',dirs)
    print('idPNs = ',idPNs)
    #STOP
    return [dirs,idPNs]


def findTargetsVisibleAt(inputList,outputPath):
    dirs, idPNs = getIDsFromFindingChartsList(inputList)

    print('dirs = ',dirs)
    print('idPNs = ',idPNs)
    if len(dirs) != len(idPNs):
        print('len(dirs)(=',len(dirs),') != len(idPNs(=',len(idPNs),')')
#        STOP
    for iDir in range(len(dirs)):
        print('dir = ',dirs[iDir],': targets = ',idPNs[iDir])

    darkTimesLocal = getDarkTimesUT()+utcoffset
    visTimes,outDirs = createVisibilityDirs(outputPath,darkTimesLocal)
    for iDir in range(len(dirs)):
        print('idPNs[',iDir,'] = ',idPNs[iDir])
        for id in idPNs[iDir]:
            pos = allPNe.find('idPNMain',id,0)[0]
            if pos < 0:
                print('ERROR: idPNMain <'+id+'> not found in allPNe')
                STOP
            targetCoord = SkyCoord(ra=float(allPNe.getData('DRAJ2000',pos))*u.deg,
                                   dec=float(allPNe.getData('DDECJ2000',pos))*u.deg,
                                   frame='icrs')
            altitudes = getTargetAltitudes(targetCoord)
            print('idPNMain = ',id,': altitudes = ',altitudes)
            print('idPNMain = ',id,': darkTimesLocal = ',darkTimesLocal)
            for iVis in range(len(visTimes)):
                altitude = altitudes[np.where(darkTimesLocal.strftime('%H:%M') == visTimes[iVis])]
                print('altitude at ',visTimes[iVis],' = ',altitude)
                if altitude > minAltitude:
                    print('dirs[iDir] = <'+dirs[iDir]+'>')
                    linkNameRoot = os.path.join(outDirs[iVis],'alt_%d_hashID_%s_' % (altitude,id))
                    print('linkNameRoot = <'+linkNameRoot+'>')
                    targetPath = '../'+dirs[iDir][dirs[iDir].rfind('/')+1:]
                    print('targetPath = <'+targetPath+'>')
                    command = 'ls '+os.path.join(outDirs[iVis],targetPath)+'/*_'+id+'_* > '+goodTargetsDir+'/temp.list'
                    os.system(command)
                    with open(os.path.join(goodTargetsDir,'temp.list'),'r') as g:
                        toLinks = g.readlines()
                    for toLink in toLinks:
                        toLink = toLink.strip()
                        linkName = linkNameRoot+toLink[toLink.rfind('/')+1:]
                        if not os.path.lexists(linkName):
                            os.symlink(toLink,linkName,target_is_directory=False)
            #STOP

def removeNotUsed():
    path = os.path.join(goodTargetsDir,'findingCharts.bak')
    fNameIDs = os.path.join(path,'hashIDs')
    with open(fNameIDs,'r') as f:
        lines = f.readlines()
    lines = [l.rstrip('\n') for l in lines]
    usedIDs = [l[l.rfind('_')+1:] for l in lines]
#    print('lines = ',lines)

    inputList = os.path.join(goodTargetsDir,'findingCharts.list')
    dirs, idPNs = getIDsFromDirList(inputList)
    print('dirs = ',dirs)
    print('idPNs = ',idPNs)

    dirs = [os.path.join(path,dir[dir.rfind('/')+1:]) for dir in dirs]
    print('dirs = ',dirs)

    for iDir in range(len(dirs)):
        for id in idPNs[iDir]:
            if id not in usedIDs:
                print('could not find ID ',id,' in usedIDs')
                pathToDelete = os.path.join(dirs[iDir],id)
                print('going to delete path <'+pathToDelete+'>')
                try:
                    shutil.rmtree(pathToDelete)
                except:
                    pass

def checkListObjectsRaDec():
    with open('/Users/azuri/daten/uni/HKU/observing/already_observed_May062020_RA_DEC.list','r') as f:
        lines = f.readlines()
    lines = [line.rstrip('\n') for line in lines]

    path = os.path.join(goodTargetsDir,'findingCharts.bak')
    fNameIDs = os.path.join(path,'hashIDs')
    print('fNameIDs = ',fNameIDs)
    with open(fNameIDs,'r') as f:
        ids = f.readlines()
    ids = [id.rstrip('\n') for id in ids]
    ids = [id[id.rfind('_')+1:] for id in ids]

    found = []
    for line in lines:
        ra1 = hmsToDeg(line[:line.find(' ')])
        dec1 = dmsToDeg(line[line.find(' ')+1:])

        for id in ids:
            pos = allPNe.find('idPNMain',id)[0]
            print('id = ',id,': pos = ',pos)
            ra2 = float(allPNe.getData('DRAJ2000',pos))
            dec2 = float(allPNe.getData('DDECJ2000',pos))

            if (angularDistancePyAsl(ra1, dec1, ra2, dec2) * 3600.) < 30.:
                print('found object with id '+id+'> in already observed!!!')
                found.append(id)
    print('found = ',found)

def removeListAFromB(listA,listB):
    with open(listA,'r') as f:
        idsA = f.readlines()
    idsA = [id.rstrip('\n') for id in idsA]
    idsA = [id[id.find('_')+1:] for id in idsA]
    idsA = [id[id.find('_')+1:] for id in idsA]
    idsA = [id[:id.find('_')] for id in idsA]
    print('idsA = ',idsA)

    with open(listB,'r') as f:
        linesB = f.readlines()
    linesB = [id.rstrip('\n') for id in linesB]
    idsB = [id[id.find('_')+1:] for id in linesB]
    idsB = [id[id.find('_')+1:] for id in idsB]
    idsB = [id[:id.find('_')] for id in idsB]
    print('idsB = ',idsB)

    for iLine in range(len(linesB)):
        line = linesB[iLine]
        if idsB[iLine] in idsA:
            print('line = ',line)
            rm = listB[:listB.rfind('/')+1]+'charts/'+line
            print('will remove file <'+rm+'>')
            os.remove(rm)

def separateTopMediumLowPriority():

    imageSources = os.path.join(goodTargetsDir,'findingCharts.list')
    srcDirs, srcidPNs = getIDsFromDirList(imageSources)
    print('srcDirs = ',srcDirs)
    print('srcidPNs = ',srcidPNs)

    visibleAtList = os.path.join(goodTargetsDir,'visible_at_dirs.list')
    visDirs, visidPNs = getIDsFromDirList(visibleAtList)
    print('visDirs = ',visDirs)
    print('visidPNs = ',visidPNs)

    dir = os.path.join(goodTargetsDir,'charts')
    for priority in ['top','middle','low']:
        pdir = os.path.join(dir,priority)

        #create directories for maximum altitude at time
        targetDirs = []
        for srcDir in srcDirs:
            createDir = os.path.join(pdir, 'maximumAltitudeAt_'+srcDir[srcDir.find('/')+1:])
            print('creating directory <'+createDir+'>')
            targetDirs.append(createDir)
            if not os.path.exists(createDir):
                os.mkdir(createDir)

        fileNames = [f for f in os.listdir(pdir) if os.path.isfile(os.path.join(pdir, f))]
        print('fileNames = ',fileNames)

        for fileName in fileNames:
            id = fileName[fileName.find('_')+1:]
            id = id[id.find('_')+1:]
            id = id[:id.find('_')]
            print('fileName = '+fileName+': id = '+id)

            #copy files from imageSources so we can create links
            for iSrcDir in range(len(srcDirs)):
                if id in srcidPNs[iSrcDir]:
                    print('found id '+id+' in '+srcDirs[iSrcDir])
                    copy_tree(os.path.join(os.path.join(imageSources[:imageSources.rfind('/')],srcDirs[iSrcDir]),id),os.path.join(targetDirs[iSrcDir],id))

        #create inputList for findTargetsVisibleAt(inputList)
        command = 'ls '+os.path.join(pdir,'maximum*')+' > '+os.path.join(pdir,'allTargetDirs.list')
        print('command = <'+command+'>')
        os.system(command)
        findTargetsVisibleAt(os.path.join(pdir,'allTargetDirs.list'),pdir)
#        maxAltDirsTemp, maxAltidPNsTemp = getIDsFromDirList(os.path.join(pdir,'allTargetDirs.list'))
#        for i in range(len(maxAltDirsTemp)):
#            print(maxAltDirsTemp[i],': ',maxAltidPNsTemp[i])
        #STOP

def getIDfromName(inputFileName):
    with open(inputFileName,'r') as f:
        names = f.readlines()
    names = [name.rstrip('\n') for name in names]
    allHASHObjectsFName = '/Users/azuri//daten/uni/HKU/observing/all_HASH_objects_full_table.csv'
    allHASHObjects = csvFree.readCSVFile(allHASHObjectsFName)
    allPNeNames = allHASHObjects.getData('Name')
    allPNeNames = [name.replace(' ','') for name in allPNeNames]
    allPNGNames = allHASHObjects.getData('PNG')
    idPNMain = allHASHObjects.getData('idPNMain')
    ids = []
    for name in names:
        pos = -1
        if name[0] != '#':
            if name[0:3] == 'PNG':
                print('checking PNG '+name[3:])
                for iPos in range(allHASHObjects.size()):
                    if name[3:].lower().replace('-','').replace('+','').replace('J','') == allPNGNames[iPos].lower().replace('-','').replace('+','').replace('J',''):
                        pos = iPos
#                    if idPNMain[iPos] == '11240':
#                        print('allPNGNames[',iPos,'] = <'+allPNGNames[iPos].lower().replace('-','').replace('+','')+'>, name[3:] = <'+name[3:].lower().replace('-','').replace('+','')+'>')
#                        print(allPNGNames[iPos].lower().replace('-','').replace('+','') == name[3:].lower().replace('-','').replace('+',''))
            else:
                for iPos in range(allHASHObjects.size()):
                    if name.lower().replace('-','').replace('+','').replace('J','') == allPNeNames[iPos].lower().replace('-','').replace('+','').replace('J',''):
                        pos = iPos
                    if False:
                        if idPNMain[iPos] == '2879':
                            print('allPNeNames[',iPos,'] = <'+allPNeNames[iPos].lower().replace('-','').replace('+','')+'>, name = <'+name.lower().replace('-','').replace('+','')+'>')
                            print(allPNeNames[iPos].lower().replace('-','').replace('+','') == name[3:].lower().replace('-','').replace('+',''))
            if pos == -1:
                print('ERROR: could not find name <'+name+'>')
                STOP
            else:
                ids.append(allHASHObjects.getData('idPNMain',pos))
    return ids

def createRaDecDir(directory,prior=False):
    alreadyObserved = '/Users/azuri//daten/uni/HKU/observing/already_observed_May062020.list'
    alreadyObservedIDs = []#getIDfromName(alreadyObserved)
    alreadyDone = []
    if prior:
        priorities = ['top','middle','low']
    else:
        priorities = ['']
    for priority in priorities:
        if priority == '':
            pdir = directory
        else:
            pdir = os.path.join(directory,priority)
        command = 'ls '+pdir+'/*??:00-??:00 > '+os.path.join(pdir,'allTargetDirs.list')
        print('command = <'+command+'>')
        os.system(command)
        with open(os.path.join(pdir,'allTargetDirs.list'),'r') as f:
            lines = f.readlines()
        lines = [line.rstrip('\n') for line in lines]
        print('lines = ',lines)
#        dirs, idPNs = getIDsFromDirList(os.path.join(pdir,'allTargetDirs.list'))
        dirs, idPNs = getIDsFromFindingChartsList(os.path.join(pdir,'allTargetDirs.list'))
        for iDir in range(len(dirs)):
            print('dir = ',dirs[iDir],': ids = ',idPNs[iDir])
            for id in idPNs[iDir]:
                print('id = ',id)
                #STOP
                done = False
                if id in alreadyObservedIDs:
                    alreadyDone.append(id)
                    done = True
                else:
                    command = 'ls '+dirs[iDir]+'/*'+id+'* > tmpDir'
                    print('command = <'+command+'>')
                    os.system(command)
                    with open('tmpDir','r') as f:
                        fNames = f.readlines()
                    fNames = [f.rstrip('\n') for f in fNames]
                    print('fNames = ',fNames)
                    for fName in fNames:
                        print('fName = <'+fName+'>')
                        print('dirs[',iDir,'] = ',dirs[iDir])
                        tmp = fName[fName.rfind('/')+1:]
                        #tmp = dirs[iDir][dirs[iDir].find('harts/')+6:]+'/'+id+'/'+fName
                        print('tmp = '+tmp)
                        if priority == '':
                            target = tmp
                        else:
                            target = os.path.join(priority,tmp)
                        print('target = '+target)
                        tmp = fName[fName.find('_')+1:]
                        tmp = tmp[tmp.find('_')+1:]
                        tmp = tmp[:tmp.find('_')]
                        print('tmp = '+tmp)
#                        tmpA = dirs[iDir][:dirs[iDir].rfind('/')]
#                        tmpA = tmpA[:tmpA.rfind('/')]
                        pos = allPNe.find('idPNMain',id)
                        print('found id '+id+' in allPNe at pos = ',pos)
                        ra = allPNe.getData('RAJ2000',pos[0])
                        dec = allPNe.getData('DECJ2000',pos[0])
                        print('ra = '+ra)
                        print('dec = '+dec)
                        #STOP
#                        tmpB = fName[fName.find('_')+1:]
#                        tmpB = tmpB[tmpB.find('_')+1:]
                        linkName = os.path.join(directory,'RAJ2000='+ra+'_DECJ2000='+dec+'_'+priority+'_HASH-ID='+id+'_'+target)
                        print('fName = '+fName+': linkName = '+linkName)
                        if not os.path.exists(fName):
                            print('could not create link to <'+fName+'>')
                            STOP
                        if not os.path.lexists(linkName):
                            os.symlink(fName,linkName)
#                        print('dir = <'+dirs[iDir]+'>')
#                        print('target = ',target,', linkName = <'+linkName+'>')
#                        command = 'ln -s '++' '+linkName
#                        print('command = ',command)
#                        STOP
#                        os.system('rm '+linkName)
#                        if 'visibility' not in linkName:
#                            os.system(command)
        #                    STOP
    print('alreadyDone = ',len(alreadyDone),': ',alreadyDone)

def removeContentInDirAFromDirB(dirA,dirB):
    print('dirA = <'+dirA+'>')
    print('dirB= <'+dirB+'>')
    command = 'ls '+dirA+' > temp'
    os.system(command)
    with open('temp','r') as f:
        lines = f.readlines()
    lines = [line.rstrip('\n') for line in lines]
    for line in lines:
        command = 'rm -r '+os.path.join(dirB,line)
        print('command = <'+command+'>')
        os.system(command)

def fix_UsrComments():
    with open(usrCommentsFile,'r') as f:
        lines = f.readlines()
    newLines = []
    try:
        plus = 0
        for i in range(len(lines)):
            if (i+plus) < len(lines):
                if (lines[i+plus].count('"') % 2) == 0:
                    newLines.append(lines[i+plus])
                else:
                    newLine = lines[i+plus].strip('\n')
                    print('i = ',i,', plus = ',plus,': newLine = <'+newLine+'>: newLine.count(") = ',newLine.count('"'))
                    while newLine.count('"') % 2 != 0:
                        plus += 1
                        newLine += ' '+lines[i+plus]
                        print('newLine = ',newLine)
                    newLine = newLine.replace('\n','')
                    newLine = newLine+'\n'
                    newLines.append(newLine)
    except:
        print('newLines = ',len(newLines),': ',newLines)
        print('i = ',i,', plus = ',plus,', len(lines) = ',len(lines))
        print('newLine = <'+newLine+'>')
        STOP
    print('len(lines) = ',len(lines))
    print('newLines = ',len(newLines),': ',newLines)
    with open(usrCommentsFile,'w') as f:
        for line in newLines:
            f.write(line)

def remove_HLA_from_fitsFiles():
    nRemoved = 0
    for i in np.arange(fitsFiles.size()-1,-1,-1):
        if fitsFiles.getData('setname',i) == 'HLAData':
            fitsFiles.removeRow(i)
            nRemoved += 1
    print('remove_HLA_from_fitsFiles: removed ',nRemoved,' entries')

def remove_not_TLPc_from_PNMain():
    nRemoved = 0
    for i in np.arange(pnMain.size()-1,-1,-1):
        if (pnMain.getData('PNstat',i) not in ['T','L','P','c']) or (pnMain.getData('domain',i) != 'Galaxy'):
            pnMain.removeRow(i)
            nRemoved += 1
    print('remove_not_TLPc_from_PNMain: removed',nRemoved,' entries. pnMain.size() = ',pnMain.size())

def get_IDs_of_objects_with_spectrum_in_literature():
    ids = []
    for i in range(usrComments.size()):
        comment = usrComments.getData('comment',i)
        if comment.find('iterature') >= 0:
            #print('i = ',i,': comment = <'+comment+'>')
            ids.append(usrComments.getData('idPNMain',i))
    return ids

def get_IDs_of_objects_which_need_better_spectra():
    ids = []
    for i in range(usrComments.size()):
        comment = usrComments.getData('comment',i)
        if comment.find('eeds') >= 0:
#            print('i = ',i,': comment = <'+comment+'>')
            ids.append(usrComments.getData('idPNMain',i))
    return ids

def remove_objects_with_spectra_from_PNMain(keepIDs,removeIDs):
    nRemoved = 0
    for i in range(fitsFiles.size()):
        idPNMain = fitsFiles.getData('idPNMain',i)
        if idPNMain not in keepIDs:
#            print('idPNMain = ',idPNMain,' has a spectrum, removing')
            foundAt = pnMain.find('idPNMain',idPNMain)
            if foundAt[0] >= 0:
                pnMain.removeRow(foundAt[0])
                nRemoved += 1
    for i in range(len(removeIDs)):
        if removeIDs[i] not in keepIDs:
            foundAt = pnMain.find('idPNMain',removeIDs[i])[0]
            if foundAt >= 0:
                pnMain.removeRow(foundAt)
                nRemoved += 1
    print('remove_objects_with_spectra_from_PNMain: removed ',nRemoved,'elements. pnMain.size() = ',pnMain.size())

def addAngDiams():
    pnMain.addColumn('MajDiam')
    for i in range(angDiams.size()):
        if angDiams.getData('InUse',i) == '1':
            found = pnMain.find('idPNMain',angDiams.getData('idPNMain',i))[0]
            if found >= 0:
                pnMain.setData('MajDiam',found,angDiams.getData('MajDiam',i))
                print('set MajDiam for idPNMain = ',angDiams.getData('idPNMain',i),' to ',pnMain.getData('MajDiam',found))
                if 'r' in pnMain.getData('MajDiam',found):
                    STOP

def addNames():
    pnMain.addColumn('Name')
    for i in range(names.size()):
        if names.getData('InUse',i) == '1':
            found = pnMain.find('idPNMain',names.getData('idPNMain',i))[0]
            if found >= 0:
                pnMain.setData('Name',found,names.getData('Name',i))
                print('set Name for idPNMain = ',names.getData('idPNMain',i),' to ',pnMain.getData('Name',found))

if __name__ == '__main__':
    if False:
    #    fix_UsrComments()
        usrComments = csvFree.readCSVFile(usrCommentsFile)
    #    remove_HLA_from_fitsFiles()
#        remove_not_TLPc_from_PNMain()
        keepIDs = get_IDs_of_objects_which_need_better_spectra()
        if '3103' in keepIDs:
            STOP
    #    for i in idPNMain_hash_no_spectrum:
    #        keepIDs.append(i)
        #removeIDs = get_IDs_of_objects_with_spectrum_in_literature()
        #for i in idPNMain_elcat_available:
        #    removeIDs.append(i)
        #for i in idPNMain_literature_available:
        #    removeIDs.append(i)
        #remove_objects_with_spectra_from_PNMain(keepIDs,removeIDs)
        #remove_elcats(needBetterIDs)
        truesWithoutSpectra = pnMain.find('PNstat','T')
        print('truesWithoutSpectra = ',len(truesWithoutSpectra))
        addAngDiams()
        addNames()
        csvFree.writeCSVFile(pnMain,allPossibleTargets)
        print('truesWithoutSpectra = ',len(truesWithoutSpectra),': ',pnMain.getData('idPNMain',truesWithoutSpectra))
        i = 0
        for idPNMain in keepIDs:
            i += 1
            print('INSERT INTO `needBetterSpectrum`(`idNBS`,`idPNMain`) VALUES ('+str(i)+','+idPNMain+');')
    if False:
        makeFindingCharts()
    if False:
#        subprocess.check_output(['ls', outDir+'/??\:*', '>', os.path.join(goodTargetsDir,'findingCharts.list')])
        inputList = os.path.join(goodTargetsDir,'targets_SAAO_2024-06-05_MGE','findingCharts.list')#os.path.join(outDir[:outDir.rfind('/')],'findingCharts.list')
        findTargetsVisibleAt(inputList,inputList[:inputList.rfind('.')])
    if False:
        dirs = os.listdir(outDir)
        print('dirs = ',dirs)
        for line in dirs:
            if 'png' not in line:
                subDirs = os.listdir(os.path.join(outDir,line.strip()))
                print('line = ',line,': subDirs = ',subDirs)
                for subDir in subDirs:
                    if 'png' not in subDir:
                        plots = os.listdir(os.path.join(outDir,line,subDir))
                        print('subDir = ',subDir,': plots = ',plots)
                        if len(plots) > 6:
                            STOP
                        for plot in plots:
                            if ('findingChart' in plot) and ('quot' not in plot) and ('rgb' not in plot) and ('iphas3colour' not in plot):
                                inputFile = os.path.join(outDir,line,subDir,plot)
                                outputFile = os.path.join(outDir,line,plot)
                                if os.path.exists(outputFile):
                                    os.remove(outputFile)
                                command = 'cp '+inputFile+' '+outputFile
                                print('command = ',command)
                                #os.system(command)
                                shutil.copyfile(os.path.join(outDir,line,subDir,plot),os.path.join(os.path.join(outDir,line),plot), follow_symlinks=True)
                                #os.symlink(os.path.join(os.path.join(os.path.join(outDir,line),subDir),plot),os.path.join(os.path.join(outDir,line),plot))
    if False:
#        dbs_may2008 = csvFree.readCSVFile('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2023-06-17/hash_fitsFiles_DBS_May2008.csv')
        removeIDPNMains = []
#        for i in range(dbs_may2008.size()):
#            if dbs_may2008.getData('setname',i) == 'DBS_May2008':
#                removeIDPNMains.append(dbs_may2008.getData('idPNMain',i))
        allFilesName = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/allFindingCharts.list'
        os.system('ls /Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/findingCharts/* > '+allFilesName)
        #os.system('ls /Users/azuri/daten/uni/HKU/observing/targets_SAAO_2024-05-06/targets_SAAO_2024-06-05_*/findingCharts/* > '+allFilesName)
        with open(allFilesName,'r') as f:
            allFiles = f.readlines()
        print('allFiles = ',allFiles)

        allFullPaths = []
        path = ''
        for i in range(len(allFiles)):
            line = allFiles[i].strip()
            if line != '':
                if ':' in line:
                    path = line[:line.rfind(':')]
                else:
                    allFullPaths.append(os.path.join(path,line))
        print('allFullPaths = ',allFullPaths)
        for fName in allFullPaths:
            if ('alt' in fName) and ('hashID' in fName):
                idPNMain = fName[fName.rfind('_')+1:]
                print('fName = ',fName,': idPNMain = ',idPNMain)
            elif ('moonDist' in fName):
                tmp = fName[:fName.rfind('_')]
                idPNMain = tmp[tmp.rfind('_')+1:]
                print('fName = ',fName,': idPNMain = ',idPNMain)
            elif len(fName[fName.rfind('/')+1]) < 6:
                idPNMain = fName[fName.rfind('/')+1:]
                print('fName = ',fName,': idPNMain = ',idPNMain)
            else:
                print('fName = ',fName)
                STOP
            if idPNMain in removeIDPNMains:
                print('removing file '+fName)
                if os.path.exists(fName):
                    os.system('rm -r '+fName)
                if os.path.islink(fName):
                    print('still exists, trying again')
                    os.unlink(fName)

            if fName == 'alt_81_hashID_3025':
                print('found alt_81_hashID_3025 -> please check')
                STOP
        #removeNotUsed()
        #checkListObjectsRaDec()
#        removeListAFromB(os.path.join(goodTargetsDir,'findingCharts.bak/priorityQ.list'),os.path.join(goodTargetsDir,'findingCharts/allCandidates.list'))
#        removeListAFromB(os.path.join(goodTargetsDir,'findingCharts.bak/rejectedQ.list'),os.path.join(goodTargetsDir,'findingCharts/allCandidates.list'))
#        separateTopMediumLowPriority()
        if False:
            createRaDecDir(os.path.join(goodTargetsDir,'charts'),prior=True)
            for priority in ['top',
                            'middle',
                            'low']:
                for time in ['00:00-01:00',
                            '01:00-02:00',
                            '02:00-03:00',
                            '03:00-04:00',
                            '04:00-05:00',
                            '05:00-06:00',
                            '19:00-20:00',
                            '20:00-21:00',
                            '21:00-22:00',
                            '22:00-23:00',
                            '23:00-24:00',]:
                    removeContentInDirAFromDirB(os.path.join(goodTargetsDir,'charts/'+priority+'/maximumAltitudeAt_'+time),
                                                os.path.join(goodTargetsDir,'findingCharts/'+time))
    createRaDecDir(outDir)
