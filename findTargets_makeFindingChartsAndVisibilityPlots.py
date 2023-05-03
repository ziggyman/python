import os
from astroplan import Observer
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time
import astropy.units as u
from distutils.dir_util import copy_tree
import numpy as np
from PIL import Image
import shutil

from plot_obs_planning import plot_target
import csvFree,csvData
from myUtils import angularDistancePyAsl,hmsToDeg,dmsToDeg

#from astroplan import download_IERS_A
#download_IERS_A()

goodTargetsDir = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good'
outDir = os.path.join(goodTargetsDir,'findingCharts')
fileList = os.path.join(goodTargetsDir,'allFiles.list')

allPossibleTargets = '/Users/azuri/daten/uni/HKU/observing/all_targets_PLc.csv'
allPNe = csvFree.readCSVFile(allPossibleTargets)

observatoryName = "SAAO"
observatoryLocation = EarthLocation(lat=-32.3783*u.deg, lon=20.8105*u.deg, height=1750*u.m)
utcoffset = 2*u.hour
date = '2020-05-15'
midnight = Time(date+' 00:00:00') - utcoffset
minAltitude = 30.0

def getIdPNMainFromFileName(fName):
    idPNTemp = fName[fName.find('_')+1:]
    idPNTemp = idPNTemp[idPNTemp.find('_')+1:]
    idPNTemp = idPNTemp[:idPNTemp.find('_')]
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
    for idPNMain in idPNMains:
        print('idPNMain = ',idPNMain)
        images = getImagesForIdPNMain(fNamesList,idPNMain)
        print('images = ',images)
        outDirName = os.path.join(outDir,idPNMain)
        if not os.path.exists(outDirName):
            os.mkdir(outDirName)

        for im in images:
            background = Image.open(os.path.join(goodTargetsDir,im))
            overlay = Image.open(os.path.join(goodTargetsDir,im[:-4]+'_centroid.png'))

            background = background.convert("RGBA")
            overlay = overlay.convert("RGBA")

            new_img = Image.new("RGBA", background.size)
            new_img = Image.alpha_composite(new_img, background)
            new_img = Image.alpha_composite(new_img, overlay)
            newImName = os.path.join(os.path.join(outDir,idPNMain),im[:-4]+'_findingChart.png').replace('.0','')
            new_img.save(newImName,"PNG")

        pos = csvAllTargets.find('idPNMain',idPNMain,0)[0]
        print('pos = ',pos)
#        print(csvAllTargets.getData(pos))

        visName = im[:im.rfind('_')]
        visName = visName[:visName.rfind('_')]
        visName = os.path.join(os.path.join(outDir,idPNMain),visName+'_visibility.png').replace('.0','')
        targetCoord = SkyCoord(ra=float(csvAllTargets.getData('DRAJ2000',pos))*u.deg,
                               dec=float(csvAllTargets.getData('DDECJ2000',pos))*u.deg,
                               frame='icrs')
        try:
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
            print('finalDir = <'+finalDir+'>')
            shutil.move(outDirName,os.path.join(finalDir,idPNMain))
        except:
            nFailed += 1
            pass
#        STOP
    print(nFailed,' objects failed to plot visibility plot')

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
            ids.append(line)
    idPNs.append(ids)
    return [dirs,idPNs]

def findTargetsVisibleAt(inputList,outputPath):
    dirs, idPNs = getIDsFromDirList(inputList)

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
        for id in idPNs[iDir]:
            pos = allPNe.find('idPNMain',id,0)[0]
            if pos < 0:
                print('ERROR: idPNMain '+id+' not found in allPNe')
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
                    linkName = os.path.join(outDirs[iVis],'alt_%d_hashID_%s' % (altitude,id))
                    targetPath = '../'+dirs[iDir][dirs[iDir].rfind('/')+1:]+'/'+id
                    print('targetPath = <'+targetPath+'>')
                    print('linkName = <'+linkName+'>')
                    if not os.path.lexists(linkName):
                        os.symlink(targetPath,linkName,target_is_directory=True)
            #STOP

def removeNotUsed():
    path = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts.bak'
    fNameIDs = os.path.join(path,'hashIDs')
    with open(fNameIDs,'r') as f:
        lines = f.readlines()
    lines = [l.rstrip('\n') for l in lines]
    usedIDs = [l[l.rfind('_')+1:] for l in lines]
#    print('lines = ',lines)

    inputList = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts.list'
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

    path = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts.bak'
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

    imageSources = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts.list'
    srcDirs, srcidPNs = getIDsFromDirList(imageSources)
    print('srcDirs = ',srcDirs)
    print('srcidPNs = ',srcidPNs)

    visibleAtList = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/visible_at_dirs.list'
    visDirs, visidPNs = getIDsFromDirList(visibleAtList)
    print('visDirs = ',visDirs)
    print('visidPNs = ',visidPNs)

    dir = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/charts'
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
    alreadyObservedIDs = getIDfromName(alreadyObserved)
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
        os.system(command)
        with open(os.path.join(pdir,'allTargetDirs.list'),'r') as f:
            lines = f.readlines()
        lines = [line.rstrip('\n') for line in lines]
        print('lines = ',lines)
        dirs, idPNs = getIDsFromDirList(os.path.join(pdir,'allTargetDirs.list'))
        for iDir in range(len(dirs)):
            print('dir = ',dirs[iDir],': ids = ',idPNs[iDir])
            for id in idPNs[iDir]:
                done = False
                if id in alreadyObservedIDs:
                    alreadyDone.append(id)
                    done = True
                else:
                    command = 'ls '+os.path.join(dirs[iDir],id)+' > tmpDir'
                    print('command = <'+command+'>')
                    os.system(command)
                    with open('tmpDir','r') as f:
                        fNames = f.readlines()
                    fNames = [f.rstrip('\n') for f in fNames]
                    print('fNames = ',fNames)
                    for fName in fNames:
                        tmp = dirs[iDir][dirs[iDir].find('harts/')+6:]+'/'+id+'/'+fName
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
                        tmpA = dirs[iDir][:dirs[iDir].rfind('/')]
                        tmpA = tmpA[:tmpA.rfind('/')]
                        pos = allPNe.find('idPNMain',id)
                        print('found id '+id+' in allPNe at pos = ',pos)
                        ra = allPNe.getData('RAJ2000',pos[0])
                        dec = allPNe.getData('DECJ2000',pos[0])
                        print('ra = '+ra)
                        print('dec = '+dec)
                        tmpB = fName[fName.find('_')+1:]
                        tmpB = tmpB[tmpB.find('_')+1:]
                        linkName = os.path.join(directory,'RAJ2000='+ra+'_DECJ2000='+dec+'_'+priority+'_HASH-ID='+tmpB)
                        print('fName = '+fName+': linkName = '+linkName)
#                        if not os.path.lexists(linkName):
#                            os.symlink(target,linkName)
                        print('dir = <'+dirs[iDir]+'>')
                        print('target = ',target)
                        command = 'cp '+os.path.join(directory,target)+' '+linkName
                        print('command = ',command)
#                        STOP
                        os.system('rm '+linkName)
                        if 'visibility' not in linkName:
                            os.system(command)
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

if __name__ == '__main__':
    makeFindingCharts()
    #inputList = '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts.list'
    #findTargetsVisibleAt(inputList,inputList[:inputList.rfind('.')])
    #removeNotUsed()
    #checkListObjectsRaDec()
    #removeListAFromB('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts.bak/priorityQ.list','/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts/allCandidates.list')
    #removeListAFromB('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts.bak/rejectedQ.list','/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts/allCandidates.list')
    #separateTopMediumLowPriority()
    #createRaDecDir('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/charts',prior=True)
    if False:
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
                removeContentInDirAFromDirB('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/charts/'+priority+'/maximumAltitudeAt_'+time,
                                            '/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts/'+time)
    createRaDecDir('/Users/azuri/daten/uni/HKU/observing/targets_SAAO_2020-05-15_good/good/findingCharts')
