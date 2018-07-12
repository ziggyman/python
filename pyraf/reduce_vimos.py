from datetime import date
import logging
import numpy as np
import os
import shlex
import shutil
import subprocess
import traceback
#from pyraf import iraf
#import pyfits
from astropy.io import fits as pyfits

pathIn = "/Volumes/obiwan/azuri/spectra/099_B-0215_A.bak"
pathOutGen = '/Volumes/obiwan/azuri/spectra/099_B-0215_A/%s/%s/%s/%s'#(date, filter, quadrant, TargetName)

def runCommand(command, outFileName=None):
    args = shlex.split(command)
    print 'runCommand: args = ',args
    if outFileName == None:
        p = subprocess.Popen(args)
        p.wait()
    else:
        lines = subprocess.check_output(args).splitlines()
        print 'runCommand: writing <'+outFileName+'>'
        with open(outFileName, 'w') as outFile:
            for line in lines:
                outFile.write(line)
    print 'runCommand: returncode = ',p.returncode

def findClosest(day, days, filter, quadrant, imType, fName):
    """find closest master imType in time"""
    for key, value in sorted(days.iteritems(), key=lambda (k,v): (v,k)):
        print "%s: %s" % (key, value)
        path = pathOutGen % (key, filter, quadrant, imType)
        masterImType = os.path.join(path,fName)
        if os.path.exists(masterImType):
            return masterImType
    return None

def sortAndMove():
    fitslist = os.path.join(pathIn,"allCompressedFits.list")
    moveFiles = True

    fitsFiles = []
    types = []
    exptimes = []
    dotOnes = []
    dotTwos = []

    with open(fitslist) as f:
        ff = list(f)

    for line in ff:
        fitsFiles.append(line.strip('\n'))
        print fitsFiles[len(fitsFiles)-1]

    pathOut = "/Volumes/obiwan/azuri/spectra/099_B-0215_A"
    if os.path.exists(pathOut):
        command = 'rm -r '+pathOut
        runCommand(command)

    nFails = 0
    for fFile in fitsFiles:
        pathOut = "/Volumes/obiwan/azuri/spectra/099_B-0215_A"
        datePos = fFile.find('VIMOS.')
        if datePos < 0:
            raise Exception('VIMOS. not found in fFile = <'+fFile+'>')
        datePos += 6
        day = fFile[datePos:datePos+10]
        print 'fFile = ',fFile
        print 'day = ',day
        print "fFile.rfind('/') = ",fFile.rfind('/')
        print "fFile[:fFile.rfind('/')] = ",fFile[:fFile.rfind('/')]
        path = fFile[:fFile.rfind('/')]
        print 'path = <'+path+'>'
        fileName = fFile[fFile.rfind('/')+1:]
        print 'fileName = ',fileName
        fNameCopy = os.path.join(path,'temp',fileName)
        print 'fNameCopy = <'+fNameCopy+'>'

        isOne = False
        fitsPos = fNameCopy.find('.fits')
        print 'fNameCopy[len(fNameCopy)-2:]  = ',fNameCopy[len(fNameCopy)-2:]
        if fNameCopy[len(fNameCopy)-2:] == '.1':
            isOne = True

        isTwo = False
        if fNameCopy[len(fNameCopy)-2:] == '.2':
            isTwo = True
        if isOne or isTwo:
            fNameCopy = fNameCopy[:fitsPos]+fNameCopy[fitsPos:len(fNameCopy)-2]
        print 'fNameCopy = <'+fNameCopy+'>'
        print 'isOne = ',isOne
        print 'isTwo = ',isTwo

        if os.path.exists(fNameCopy):
            os.remove(fNameCopy)
        shutil.copy(fFile, fNameCopy)

        print 'fNameCopy[len(fNameCopy)-2:] = <'+fNameCopy[len(fNameCopy)-2:]+'>'
        if fNameCopy[len(fNameCopy)-2:]=='.Z':
            fNameCopy = fNameCopy[:len(fNameCopy)-2]
            if os.path.exists(fNameCopy):
                os.remove(fNameCopy)
            print 'uncompressing <'+fNameCopy+'>'
            command = 'uncompress '+fNameCopy
            runCommand(command)
    #    shutil.move(fNameCopy, pathOut)
        try:
            hdulist = pyfits.open(fNameCopy, ignore_missing_end=True)
    #        print hdulist.info()
            print 'reading header'
            header = hdulist[0].header
    #        print 'header = ',header
            if 'OBJECT' in header:
                print 'header[OBJECT] = ',header['OBJECT']
            else:
                raise Exception('OBJECT not in header')
            if header['OBJECT'] not in types:
                types.append(header['OBJECT'])
            pathOut = os.path.join(pathOut,day)
            if isOne:
                dotOnes.append(header['OBJECT'])
            if isTwo:
                dotTwos.append(header['OBJECT'])

            if 'EXPTIME' in header:
                print 'header[EXPTIME] = ',header['EXPTIME']
            else:
                print 'keyWord <EXPTIME> not in header'
            if header['EXPTIME'] not in exptimes:
                exptimes.append(header['EXPTIME'])

            if 'NAXIS1' in header:
                print 'header[NAXIS1] = ',header['NAXIS1']
            else:
                print 'keyWord <NAXIS1> not in header'

            if 'NAXIS2' in header:
                print 'header[NAXIS2] = ',header['NAXIS2']
            else:
                print 'keyWord <NAXIS2> not in header'

            if 'ORIGFILE' in header:
                origFileName = header['ORIGFILE']
                print 'origFileName = <'+origFileName+'>'
            else:
                print 'keyword <ORIGFILE> not found in header'

            quadrantFound = 0
            filterName = 'noFilter'
            for quadrant in np.arange(1,5):
                if 'HIERARCH ESO INS FILT'+str(quadrant)+' NAME' in header:
                    quadrantFound = quadrant
                    filterName = header['HIERARCH ESO INS FILT'+str(quadrant)+' NAME']
                    print 'header[HIERARCH ESO INS FILT'+str(quadrant)+' NAME] = ',header['HIERARCH ESO INS FILT'+str(quadrant)+' NAME']
                else:
                    print 'keyWord <HIERARCH ESO INS FILT1 NAME> not in header'
            print 'quadrantFound = ',quadrantFound
            pathOut = os.path.join(pathOut,filterName)
            if quadrantFound > 0:
                pathOut = os.path.join(pathOut,str(quadrantFound))
                print 'quadrantFound = ',quadrantFound

            if 'HIERARCH ESO OBS TARG NAME' in header:
                print 'header[HIERARCH ESO OBS TARG NAME] = ',header['HIERARCH ESO OBS TARG NAME']
            else:
                print 'keyWord <HIERARCH ESO OBS TARG NAME> not in header'
            pathOut = os.path.join(pathOut,header['OBJECT'])
            pathOut = pathOut.replace(' ','_')
            print 'pathOut = <'+pathOut+'>'
            if moveFiles:
                if not os.path.exists(pathOut):
                    os.makedirs(pathOut)
                fNameOut = os.path.join(pathOut, fNameCopy[fNameCopy.rfind('/')+1:])
                if 'HIERARCH ESO OBS TARG NAME' in header:
                    fNameOut = fNameOut[:fNameOut.rfind('.')]+'_'+header['HIERARCH ESO OBS TARG NAME']+fNameOut[fNameOut.rfind('.'):]
                if isOne:
                    fNameOut = fNameOut[:fNameOut.rfind('.')]+'_1'+fNameOut[fNameOut.rfind('.'):]
                if isTwo:
                    fNameOut = fNameOut[:fNameOut.rfind('.')]+'_2'+fNameOut[fNameOut.rfind('.'):]
                fNameOut = fNameOut.replace(' ','_')
                fNameOut = fNameOut.replace(':','-')
                print 'fNameOut = <'+fNameOut+'>'
                if os.path.exists(fNameOut):
                    os.remove(fNameOut)
                shutil.move(fNameCopy, fNameOut)

        except Exception as e:
            nFails += 1
            print 'failed to read ',fNameCopy,': '
            logging.error(traceback.format_exc())
            STOP
        print ' '
        print ' '
    #    STOP
    print ' '
    print ' '
    print 'failed to read ',nFails,' files'
    print 'types = ',types
    print 'exptimes = ',exptimes
    print 'dotOnes = ',dotOnes
    print 'dotTwos = ',dotTwos
    #    for keyWord in header:
    #        print 'header[',keyWord,'] = ',header[keyWord]

def reduce(day, days, filter, quadrant):
    print 'reducing ',day,filter,quadrant

    if filter == 'Free':
        path = pathOutGen % (day, filter, quadrant, 'BIAS')
        print 'path = <'+path+'>'
        if os.path.exists(path):

            """create bias list and Set Of Frames (sof)"""
            biasList = os.path.join(path,'bias.list')
            biasSOF = os.path.join(path,'bias.sof')
            print 'biasList = <'+biasList+'>'
            command = 'ls '+path+'/*.???.fits > '+biasList
            os.system(command)

            with open(biasList) as f:
                ff = list(f)

            biasFrames = []
            for line in ff:
                biasFrames.append(line.strip('\n') + ' BIAS\n')
            if len(biasFrames) > 0:
                with open(biasSOF, 'w') as sof:
                    for line in biasFrames:
                        sof.write(line)

                """create master bias"""
                command = 'esorex --output-dir='+path+' vmbias --StackMethod=Median '+biasSOF
                runCommand(command)
            else:
                print 'no Biases found for day '+day+', filter '+filter+', quadrant '+quadrant
        else:
            print 'no Biases found for day '+day+', filter '+filter+', quadrant '+quadrant

    else:
        calendarDate = date(int(day[:day.find('-')]),
                            int(day[day.find('-')+1:day.rfind('-')]),
                            int(day[day.rfind('-')+1:]))
        print 'calendarDate = ',calendarDate
        dates = {}
        for dayTemp in days:
            print 'dayTemp = ',dayTemp
            calendarDateTemp = date(int(dayTemp[:dayTemp.find('-')]),
                                    int(dayTemp[dayTemp.find('-')+1:dayTemp.rfind('-')]),
                                    int(dayTemp[dayTemp.rfind('-')+1:]))
            print 'calendarDateTemp = ',calendarDateTemp
            dates[dayTemp] = abs((calendarDate - calendarDateTemp).days)
        print 'dates = ',dates

        """create master sky flat"""
        path = pathOutGen % (day, filter, quadrant, 'FLAT,SKY')
        print 'path = <'+path+'>'

        """create Flat list and Set Of Frames (sof)"""
        if os.path.exists(path):
            flatList = os.path.join(path,'flat.list')
            flatSOF = os.path.join(path,'flat.sof')
            print 'flatList = <'+flatList+'>'
            command = 'ls '+path+'/*_??.fits > '+flatList
            os.system(command)

            """add master bias to flatList"""
            masterBias = findClosest(day, dates, 'Free', quadrant, 'BIAS', 'master_bias.fits')
            if not masterBias:
                raise Exception('creating master Flat: ERROR: no masterBias found')
            command = 'ls '+masterBias+' >> '+flatList
            os.system(command)

            with open(flatList) as f:
                ff = list(f)

            """add tag to each image"""
            flatFrames = []
            for line in ff:
                line = line.strip('\n')
                if line[line.rfind('/')+1:] == 'master_bias.fits':
                    flatFrames.append(line + ' MASTER_BIAS\n')
                else:
                    flatFrames.append(line + ' IMG_SKY_FLAT\n')
            if len(flatFrames) > 0:
                with open(flatSOF, 'w') as sof:
                    for line in flatFrames:
                        sof.write(line)

                """create master flat"""
                command = 'esorex --output-dir='+path+' vmimflatsky --StackMethod=Median '+flatSOF
                runCommand(command)
            else:
                print 'no Flats found for day '+day+', filter '+filter+', quadrant '+quadrant
        else:
            print 'no Flats found for day '+day+', filter '+filter+', quadrant '+quadrant

        """reduce standard fields"""
        path = pathOutGen % (day, filter, quadrant, 'STD')
        print 'path = <'+path+'>'

        if os.path.exists(path):
            """create individual lists for each STD and Set Of Frames (sof)"""
            stdList = os.path.join(path,'std.list')
            print 'stdList = <'+stdList+'>'
            command = 'ls '+path+'/*.fits > '+stdList
            os.system(command)

            """add master bias to stdList"""
            masterBias = findClosest(day, dates, 'Free', quadrant, 'BIAS', 'master_bias.fits')
            if not masterBias:
                raise Exception('reducing standard star: ERROR: no masterBias found')

            """add master flat to stdList"""
            masterFlat = findClosest(day, dates, filter, quadrant, 'FLAT,SKY', 'img_master_sky_flat.fits')
            if not masterFlat:
                raise Exception('reducing standard star: ERROR: no master Flat found')

            with open(stdList) as f:
                ff = list(f)

            """add tag to each image"""
            stdFrames = []
            for line in ff:
                line = line.strip('\n')
                stdFrames.append(line + ' IMG_STANDARD\n')
            if len(stdFrames) > 0:
                stdSOF = ' '
                for frame in stdFrames:
                    stdSOF = os.path.join(path,frame[frame.rfind('/')+1:frame.rfind('.')]+'.sof')
                    with open(stdSOF, 'w') as sof:
                        sof.write(frame)
                        sof.write(masterBias+' MASTER_BIAS\n')
                        sof.write(masterFlat+' IMG_MASTER_SKY_FLAT\n')
                        sof.write('/Users/azuri/eso-pipelines/vimos/calib/vimos-3.2.3/cal/phstd_stetson.tfits PHOTOMETRIC_CATALOG\n')
                        sof.write('/Users/azuri/eso-pipelines/vimos/calib/vimos-3.2.3/cal/badpixel.'+quadrant+'.tfits CCD_TABLE\n')
                        sof.write('/Users/azuri/eso-pipelines/vimos/calib/vimos-3.2.3/cal/ipc_'+filter+'.'+quadrant+'.tfits PHOTOMETRIC_TABLE')

                """reduce standard field"""
                command = 'esorex --output-dir='+path+' vmimstandard --CleanBadPixel=TRUE '+stdSOF
                runCommand(command)

                """create phot_coeff_table.fits"""
                photSOF = os.path.join(path,'phot.sof')
                with open(photSOF,'w') as sof:
                    sof.write(os.path.join(path,'img_star_match_table.fits')+' IMG_STAR_MATCH_TABLE\n')
                    sof.write('/Users/azuri/eso-pipelines/vimos/calib/vimos-3.2.3/cal/ipc_'+filter+'.'+quadrant+'.tfits PHOTOMETRIC_TABLE')
                command = 'esorex --output-dir='+path+' vmimcalphot '+photSOF
                runCommand(command)
                stdPath = path
            else:
                print 'no standard stars found for day '+day+', filter '+filter+', quadrant '+quadrant
        else:
            print 'no standard stars found for day '+day+', filter '+filter+', quadrant '+quadrant

        """reduce science frames"""
        path = pathOutGen % (day, filter, quadrant, 'Kathryns_Wheel')
        print 'path = <'+path+'>'

        if os.path.exists(path):
            """create individual lists for each STD and Set Of Frames (sof)"""
            scienceList = os.path.join(path,'science.list')
            print 'scienceList = <'+scienceList+'>'
            command = 'ls '+path+'/*l.fits > '+scienceList
            os.system(command)

            """add master bias to scienceList"""
            masterBias = findClosest(day, dates, 'Free', quadrant, 'BIAS', 'master_bias.fits')
            if not masterBias:
                raise Exception('reducing science frames: ERROR: no masterBias found')
            print 'closest master Bias = <'+masterBias+'>'

            """add master flat to scienceList"""
            masterFlat = findClosest(day, dates, filter, quadrant, 'FLAT,SKY', 'img_master_sky_flat.fits')
            if not masterFlat:
                raise Exception('reducing science frames: ERROR: no master Flat found')
            print 'closest master Flat = <'+masterFlat+'>'

            photometricTable = findClosest(day, dates, filter, quadrant, 'STD', 'photometric_table.fits')
            if not photometricTable:
                raise Exception('reducing science frames: ERROR: no photometric table found')
            print 'closest photometricTable = <'+photometricTable+'>'

            with open(scienceList) as f:
                ff = list(f)

            """add tag to each image"""
            scienceFrames = []
            for line in ff:
                line = line.strip('\n')
                scienceFrames.append(line + ' IMG_SCIENCE\n')

            if len(scienceFrames) > 0:

                scienceSOF = os.path.join(path,frame[frame.rfind('/')+1:frame.rfind('.')]+'.sof')
                with open(scienceSOF, 'w') as sof:
                    for frame in scienceFrames:
                        sof.write(frame)
                    sof.write(masterBias+' MASTER_BIAS\n')
                    sof.write(masterFlat+' IMG_MASTER_SKY_FLAT\n')
                    sof.write('/Users/azuri/eso-pipelines/vimos/calib/vimos-3.2.3/cal/badpixel.'+quadrant+'.tfits CCD_TABLE\n')
                    sof.write(photometricTable+' PHOTOMETRIC_TABLE\n')
                    #sof.write('/Users/azuri/eso-pipelines/vimos/calib/vimos-3.2.3/cal/ipc_R.1.tfits PHOTOMETRIC_TABLE')

                outDir = os.path.join(path,'calibrated')
                if not os.path.exists(outDir):
                    os.makedirs(outDir)

                """reduce science frames with standard calibration"""
                command = 'esorex --output-dir='+outDir+' vmimobsjitter --CleanBadPixel=TRUE '+scienceSOF
                runCommand(command)

                scienceSOF = os.path.join(path,frame[frame.rfind('/')+1:frame.rfind('.')]+'.sof')
                with open(scienceSOF, 'w') as sof:
                    for frame in scienceFrames:
                        sof.write(frame)
                    sof.write(masterBias+' MASTER_BIAS\n')
                    sof.write(masterFlat+' IMG_MASTER_SKY_FLAT\n')
                    sof.write('/Users/azuri/eso-pipelines/vimos/calib/vimos-3.2.3/cal/badpixel.'+quadrant+'.tfits CCD_TABLE\n')
                    #sof.write(os.path.join(stdPath,'photometric_table.fits PHOTOMETRIC_TABLE\n'))
                    sof.write('/Users/azuri/eso-pipelines/vimos/calib/vimos-3.2.3/cal/ipc_R.1.tfits PHOTOMETRIC_TABLE')

                outDir = os.path.join(path,'uncalibrated')
                if not os.path.exists(outDir):
                    os.makedirs(outDir)

                """reduce science frames"""
                command = 'esorex --output-dir='+outDir+' vmimobsjitter --CleanBadPixel=TRUE '+scienceSOF
                runCommand(command)
            else:
                print 'No Science frames found for day '+day+', filter '+filter+', quadrant '+quadrant
        else:
            print 'No Science frames found for day '+day+', filter '+filter+', quadrant '+quadrant

sortAndMove()
days = ['2017-08-14', '2017-08-15', '2017-08-16', '2017-08-19']
filters = ['B', 'V']
quadrants = ['1', '2', '3', '4']
for day in days:
    for quadrant in quadrants:
        reduce(day, days, 'Free', quadrant)
        for filter in filters:
            reduce(day, days, 'B', quadrant)
