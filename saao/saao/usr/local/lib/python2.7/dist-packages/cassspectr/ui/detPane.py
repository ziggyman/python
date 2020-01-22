import time
from PyQt4 import QtCore, QtGui

import logging as log
from cassspectr.controller.controller import Controller
import os
import json
import numpy as np
from astropy.io import fits

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


class DetPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl
        self.lastCCDStatus = 0

    def frameNumPrefixChange(self):
        self.settings['propID'] = str(self.ui.lineEditPropID.text())
        self.ui.displFrameNumPrefix.setText('a'+self.settings['propID'])
        # use this to trigger reset of frame numbers (if this run num is different from that stored in settings file):
##        self.logger.info('New run number entered. Starting image number reset to 1001!')


    def confirmNewRun(self):
        # on the first night of the run, this will reset the frame numbering to 1001
        # warn the user with a confirm box to make sure they know that this will
        # reset the number, potentially overwriting data
        check_msg = """Is this the first night of your run?

Answering yes will reset the image numbering.
If you have already taken data in this directory
you may potentially overwrite it!!
"""
        reply = QtGui.QMessageBox.question(None, 'Removal', check_msg,
                                           QtGui.QMessageBox.Yes | QtGui.QMessageBox.No,
                                           QtGui.QMessageBox.No)
        if (reply == QtGui.QMessageBox.Yes):
            self.logger.info('New run number entered. Starting image number reset to 1001!')
            self.ui.displFrameNum.setText('1001')
        else:
            pass


    def displDataDir(self):
        self.ui.displdatadir.setText('datadir:    '+self.datadir)


    def diskUsage(self):
        rootdatadir='/home/'
        foo = os.popen('df -h %s'%rootdatadir)
        f = foo.readline()
        f = foo.readline()
        t = f.split()[-2]
        percentageUsed = float(t.replace('%',''))
        avail = f.split()[-3]
        self.ui.progressBarDiskUsage.setProperty("value", percentageUsed)
        self.ui.labTextDiskUsage_3.setText(avail.replace(',','.'))

    
    def setBinning(self):
        # if binning pars have changed, update in dict:
        # get binning:
        binString = self.ui.comboCCDbinning.currentText()
        xbin = binString.split('x')[0]
        ybin = binString.split('x')[1]
        if (self.dbglvl>2):
            print 'binning: %s x %s'%(xbin,ybin)
        self.settings['nColBin'] = int(xbin)
        self.settings['nRowBin'] = int(ybin)

        binStr = ' %s x %s'%(self.settings['nColBin'],self.settings['nRowBin'])
        self.ui.labCCDBinningState.setText(str(binStr))


    def setCCDmode(self):
        # if CCD mode (gain or noise) changes, update in dict (and set on detector):
        self.settings['CCDmodeGainIndx'] = self.ui.comboCCDmodeGain.currentIndex()
        self.settings['CCDmodeNoiseIndx'] = self.ui.comboCCDmodeNoise.currentIndex()

        self.con.dc.set_gain_index(self.settings['CCDmodeGainIndx'])
        self.con.dc.set_readout_speed_index(self.settings['CCDmodeNoiseIndx'])
        # update text on GUI:
        gains_table = self.con.dc.get_gains()
        #print("The gain is {}".format(gains_table[ro][gn]["gain"]))
        #print("The noise in electrons is {}".format(gains_table[ro][gn]["noise_e"]))
        gn = self.settings['CCDmodeGainIndx']
        ro = self.settings['CCDmodeNoiseIndx']
        gainval = gains_table[ro][gn]["gain"]
        noiseval = gains_table[ro][gn]["noise_e"]
        self.ui.displTextCCDMode.setText('gain= %.2f; RN=%.2fe-'%(gainval,noiseval))
        self.ui.labCCDmodeStateGain.setText(str(self.ui.comboCCDmodeGain.currentText()))
        self.ui.labCCDmodeStateNoise.setText(str(self.ui.comboCCDmodeNoise.currentText()))


    def setCCDwindowing(self):
        if(self.ui.comboCCDwindowing.currentText()[0:3]=='ENG'):
            self.ccdFullFrame = True
            self.logger.info('''CCD set to full frame ENGINEERING mode.
            Please remember to change back after use.''')
        else: #SCIENCE:
            self.ccdFullFrame = False
            if(self.dbglvl>1):
                print 'CCD windowing set to SCIENCE mode'


    # This routine currently handles polling of CCD status and
    # essentially handles end of exposure, etc.
    def drawExp(self):
        try: 
            self.lastCCDStatus = self.lastCCDStatus*1
        except:
            self.lastCCDStatus=-1

        if(self.dbglvl>2):
            print self.lastCCDStatus

        if(self.dbglvl>4): print 'd1: self.con.dc.get_exposure_status()'
        try:
            status = self.con.dc.get_exposure_status()
        except Exception as e:
            print("Caught detector controller exception: {}".format(e))
            self.ui.frameCCD.setLineWidth(5)
            self.ui.frameCCD.setPalette(self.ui.pWarn)
            return

        if(self.dbglvl>2):
            print 'status = %s'%status

        # draw exposure time-related indicators:
        if (status==0): # idle
            self.ui.lcdNumberExpTime.display(float(self.settings['expTime']))

            if (self.lastCCDStatus==2):
                # we've just finished a readout:
                self.doneReadout()
                self.detInterlock=1 # safe
            ###remaining_time = float(self.ui.displExpTime.text())
            remaining_time = self.settings['expTime']
            self.sw_remaining_time = float(self.settings['expTime'])
            ###

            # abort and stop buttons should be disabled:
#            self.ui.pushButAbort.setEnabled(False)
#            self.ui.pushButStop.setEnabled(False)

            percentageRem=100.0
            self.lastCCDStatus = 0

        elif (status == 1): # exposing
            self.detInterlock=0 # LOCK!
            remaining_time = self.con.dc.get_remaining_time()
            if(self.dbglvl>2):
                print remaining_time, self.polltime, self.startTime
            if(remaining_time<=0.0 ):
                # we have just finished exposing!
                self.doneExposure()

            percentageRem = remaining_time/float(self.settings['expTime'])*100.
            self.lastCCDStatus = 1

##            self.sw_remaining_time = float(self.sw_remaining_time) - float(self.pollIncrSec)
#--
            self.sw_remaining_time = remaining_time
#--
            if(self.sw_remaining_time<0.):
                print 'sw_remaining_time has gone wrong - attempting to fix'
                self.sw_remaining_time = self.settings['expTime']
                # Assume it went wrong at the start. We could maybe
                # check difference between start time and now to
                # correct, but this should be better than allowing to
                # run -ve

            elapsedTime = np.floor(float(self.settings['expTime'])-self.sw_remaining_time)
#            self.ui.lcdNumberExpTime.display(np.floor(self.sw_remaining_time))
            self.ui.lcdNumberExpTime.display(np.floor(self.sw_remaining_time))

            self.ui.pushButAbort.setEnabled(True)

        elif (status==2): # readout
            self.detInterlock=1 # should be safe to move PLC settings during
                                # readout (but don't request CCD temp!)
            if ( (self.lastCCDStatus==1) | (self.lastCCDStatus==0) ):
                # we've just started reading out:
                self.readoutStarted()

            nPixRead = self.con.dc.get_pixel_count()

            if (self.dbglvl>2):
                #*** Hmm, don't know what's going on here!
                print nPixRead, np.sum(self.con.dc.get_pixel_count())
                
            nPixTot = self.settings['nRows'] * self.settings['nCols']
            percentageRem = 100.*float(nPixRead)/float(nPixTot)
            self.lastCCDStatus = 2
            # estimate readout time based on past-averaged speed:
            tt = time.time()-self.ROstartTime
            percentageDone = percentageRem # actually!
            total_time = 100./percentageDone * tt  
            remaining_time = total_time - tt
            # won't be quite perfectly in sync with elapsed time, but good enough for display purposes:

        try:
            self.ui.displTimeRem.setText("%.0f"%(np.floor(elapsedTime)))
        except:
            self.ui.displTimeRem.setText("0")
        self.ui.progressBarExpProgress.setProperty("value", percentageRem)
        

    # called by pressing the EXPOSE button:
    def doFirstExposure(self):
        # -- get number of exposures when EXPOSE button is pressed:
        self.NexpRemaining = int(self.ui.displNexp.value())
        if (self.NexpRemaining>1):
            self.ui.pushButCancelSequence.setEnabled(True)
        print '%s exposures requested'%self.NexpRemaining 
        self.ui.displNexpRem.setText(str(self.NexpRemaining))
        self.startExposure()


    def startExposure(self):
        # set some info (for header keywords) at start of exposure:
        headDictStart = self.createStartHeader()
        # Add settings as fits_info:
        for k,v in headDictStart.items():
            self.con.fits_info[k] = v
        # this should still be there for writing fits_info at END of exposure

        print('Exposing!')
        self.logger.info('exposing...   %s'%self.settings['Nexp'])

        exptime = float(self.settings['expTime'])
        nRows = self.settings['nRows']
        nRowBin = self.settings['nRowBin']
        nRowCen = self.settings['nRowCen']
        nCols = self.settings['nCols']
        nColBin = self.settings['nColBin']
        nColCen = self.settings['nColCen']

        # switch to use getTime() routine, which can be changed as
        # required to make more accurate later:
        self.getTime()
        self.expStartTime = self.timeSec # need to copy??***
        if(self.ccdFullFrame):
            print 'FULL FRAME exposure. DO NOT USE FOR SCIENCE'
            self.con.dc.start_engineering_exposure(exptime)
        else: # science mode!
            self.con.start_exposure(exptime, nRows, nCols, nRowBin, nColBin, nRowCen, nColCen, supp_fits_info = self.con.fits_info)

        # set states at start of exposure: 
        # un-grey 'abort' button
        # grey out Expose button
        # start count-down timer, etc.
        if(self.dbglvl>2):
            print "abort button enabled"
###        self.ui.pushButStop.setEnabled(True)
        self.ui.pushButAbort.setEnabled(True)
        self.ui.pushButExpose.setEnabled(False)

        # Grab certain one-time info at START of exposure
        self.startTime = self.polltime
        self.reqNumExps = int(self.ui.displNexp.value())
        startUT = time.strftime("%H:%M:%S", time.gmtime())
        self.ui.displUTstart.setText(str(startUT))

        self.sw_remaining_time = exptime # keep track of expected remaining time in software.
                                         # Check against CCD output for problems
        self.ui.pushButAbort.setEnabled(True)
        self.ui.pushButExpose.setEnabled(False)
        # prevent exposure time being edited

        if (self.dbglvl>1):
            print("""Command sent:
start_exposure ({}{}{}{}{}{}{})""".format(
    exptime, nRows, nCols, nRowBin, nColBin, nRowCen, nColCen))
        if (self.dbglvl>2):
            print '**',self.con.dc.get_exposure_status()


    def cancelSequence(self):
        # cancel the remaining queued exposures. Don't do anything to current exposure.
        self.NexpRemaining = 0
        self.logger.info('Cancelled remaining exposures in sequence.')
        self.ui.displNexpRem.setText(str(self.NexpRemaining))


    def PlayAlert(self):
        # play sound effect at end of exposure. Could make sound file user-choosable...
#        alertFile = '/usr/share/sounds/purple/alert.wav'
        alertFile = '/usr/share/sounds/freedesktop/stereo/onboard-key-feedback.oga'
        try:
            os.system('play '+alertFile)
        except:
            pass # might not have sox installed (for play command)
                 # or might not be able to find alertFile
        if(self.dbglvl>0):
            print 'playing %s'%alertFile


    def readoutStarted(self):
        # for short exposures, lastCCDstatus may never enter "1"
        ### this is effectively doneExposure()
        # let's leave here for now:
        if(self.NexpRemaining>0): # careful: if we abort sequence during readout,
                                  # we'll end up with -1 exps remaining!
            self.NexpRemaining -= 1
        
        self.ui.displNexpRem.setText(str(self.NexpRemaining))
        # disable abort and pause buttons:
        self.ui.pushButAbort.setEnabled(False)
        self.ui.pushButStop.setEnabled(False)

        self.getTime()
        self.expEndTime = self.timeSec #** copy??
        self.expTimeActual = self.expEndTime - self.expStartTime
        if (self.dbglvl>1):
            # this is probably not as accurate as requested time(?) unless exposure is stopped:
            print 'actual exposure time was: %.1f'%self.expTimeActual
        startUT = time.strftime("%H:%M:%S", time.gmtime())

        self.ui.displUTstart.setText(str(startUT))
        self.ROstartTime = time.time()
        self.ui.labTextExpStarted.setText('Readout started')
        self.logger.info('CCD reading out...')
        self.lastT = 0
        self.doneExposure()


    def doneExposure(self):
        if(self.settings['soundEndExposure']):
            self.PlayAlert()


    def doneReadout(self): # is really "done_readout"
        # after receiving signal that exposure has finished, return to previous state
#        self.ui.displExpTime.setText(str(self.reqExpTime))
        self.ui.displExpTime.setText(str(self.settings['expTime']))
        self.ui.pushButAbort.setEnabled(False)
        self.ui.pushButExpose.setEnabled(True)
        self.ui.pushButCancelSequence.setEnabled(True)
        self.ui.labTextExpStarted.setText('exposure started')

        self.sw_remaining_time = 0.

        if (self.settings['expType']=='TEST'):
            if (self.dbglvl>1):
                print 'TEST image'
            self.ui.displUTstart.setText('')
            self.logger.info('exposure complete. CCD ready')

            # -- write image:
            image = self.con.dc.get_image_data()
            outimage = self.datadir+'/test.fits'
        else:
            # increment frame number:
            framenum = int( self.ui.displFrameNum.text() )
            framenum +=1

            self.ui.displFrameNum.setText(str(framenum))
            self.ui.displUTstart.setText('')
            self.logger.info('exposure complete. CCD ready')

            # set displayed exptime back to original value (otherwise done
            # when IDLE; sequence won't reach idle):
            self.ui.lcdNumberExpTime.display(float(self.settings['expTime']))

            self.HartmannExposureWaiting=False

            # -- write image:
            image = self.con.dc.get_image_data()
            outimage = self.datadir+'a%s%s.fits'%(self.settings['propID'],self.settings['frameNum'])
            print 'writing to %s'%outimage
            # Check that file does not already exist - PREVENT OVERWRITING:
            #if os.path.exists(outimage):
            #    self.logger.warning('Output image %s already exists!')
#                outimage = outimage.replace()
            #self.logger.info('written a%s%s.fits'%(self.settings['propID'],self.settings['frameNum']))

        # Add settings as fits_info:
        headDict = self.createHeader()
        for k,v in headDict.items():
            self.con.fits_info[k] = v
        self.con.gen_fits_with_keys(image, outimage)

        print image.shape
        self.image = image
        
        if(1):
            self.quickView()
        
        if(self.settings['soundEndReadout']):
            self.PlayAlert()
        
        if (self.NexpRemaining>1):
            self.ui.pushButCancelSequence.setEnabled(True)
        else:
            self.ui.pushButCancelSequence.setEnabled(False)

        # check if we need to start another exposure:
        if (self.NexpRemaining>0):
            self.startExposure()

        
    def stopExposure(self):
        #*** stop not implemented in Carel's code yet??
        self.con.dc.finish_exposure()
        self.logger.info('exposure stop request sent')

        self.sw_remaining_time = 0.
        self.ui.pushButExpose.setEnabled(False)
        self.ui.pushButStop.setEnabled(False)


    def doAbortExposure(self):
        print 'aborting!'
        # *** send an abort command to Steve's code and wait for reply ****
        self.con.dc.abort_exposure()
        self.logger.info('abort request sent')
        # disable check for sw_remaining_time (or will report failure):
        self.sw_remaining_time = 0.
        self.abortPollTime = self.polltime
        self.ui.pushButExpose.setEnabled(False)
        if(self.dbglvl>2):
            status = self.con.dc.get_exposure_status()
            print 'status = %s'%status
        # wait for successful completion, then:
        self.NexpRemaining = self.NexpRemaining -1
        self.doneAbort()
        
    def doneAbort(self):
        # after receiving signal that abort has completed,
        # disable abort button
        # re-enable Expose button
        self.ui.pushButAbort.setEnabled(False)
        self.ui.pushButStop.setEnabled(False)
        self.ui.pushButExpose.setEnabled(True)
        # turn off interlock
        self.detInterlock = 1 

        self.logger.info('exposure aborted')
        self.doneExposure()

    def doAbortCleanup(self):
        # Suppress CCD Failure warning and set NexpRem = 0 after a scuucessful abort:
        self.ui.frameCCD.setLineWidth(0)
        self.ui.frameCCD.setPalette(self.ui.pNorm)
        self.drawDetIndicators()
        self.detInterlock=1 # safe

    def flushCCDPane(self):
        self.doAbortExposure()
        self.doneAbort()
        self.ui.CCDtabWidget.setCurrent(0)
        self.doAbortCleanup()
        self.detInterlock=1 # safe


    def setExpTime(self):
        if(self.dbglvl>0):
            print 'expType changed'
            print str(self.ui.comboExpType.currentText())
        # automatically set exposure time (and Nexp) based on exposure type:
        if(self.ui.comboExpType.currentText()=='BIAS'):
            self.ui.displExpTime.setText('0.0')
            self.ui.displNexp.setValue(self.settings['lastBiasNexp'])
        elif(self.ui.comboExpType.currentText()=='ARC'):
            self.ui.displExpTime.setText('%.1f'%self.settings['lastArcExpTime'])
            self.ui.displNexp.setValue(self.settings['lastArcNexp'])
            self.setupArc()
        elif(self.ui.comboExpType.currentText()=='SCIENCE'):
            self.ui.displExpTime.setText('%.1f'%self.settings['lastScienceExpTime'])
            self.ui.displNexp.setValue(self.settings['lastScienceNexp'])
            self.unsetArc()
        elif(self.ui.comboExpType.currentText()=='FLAT'):
            self.ui.displExpTime.setText('%.1f'%self.settings['lastFlatExpTime'])
            self.ui.displNexp.setValue(self.settings['lastFlatNexp'])
            self.unsetArc()
        elif(self.ui.comboExpType.currentText()=='HARTMANN'):
            #self.ui.displExpTime.setText('%.1f'%self.settings['lastFlatExpTime'])
            # set exptime in setupHartmann (need to check which lamp first)
            self.ui.displNexp.setValue(1)
            self.setupHartmann()
        elif(self.ui.comboExpType.currentText()=='TEST'):
            self.ui.displNexp.setValue(1.0) # assume we only want 1 test exposure
                                            # leave exptime unchanged
            # I think settings['expType'] is changed to 'TEST' automatically.
            # Use this as check later when writing frame


    def expTimeChanged(self):
        # store the most recent value of exptime for this exptype:
        print ';',self.ui.comboExpType.currentText()
        settings = self.settings
        if(self.ui.comboExpType.currentText()=='SCIENCE'):
            self.settings['lastScienceExpTime'] = float(self.ui.displExpTime.text())

        elif(self.ui.comboExpType.currentText()=='ARC'):
            settings['lastArcExpTime'] = float(self.ui.displExpTime.text())

            if(self.dbglvl>2):
                print  self.settings['lastArcExpTime']
                print settings
        elif(self.ui.comboExpType.currentText()=='FLAT'):
            self.settings['lastFlatExpTime'] = float(self.ui.displExpTime.text())

        self.settings = settings


    def nExpChanged(self):
        # store the most recent value of Nexp for this exptype:
        print ';',self.ui.comboExpType.currentText()
        print '::',self.ui.displNexp.value()
        settings = self.settings
        if(self.ui.comboExpType.currentText()=='BIAS'):
            self.settings['lastBiasNexp'] = float(self.ui.displNexp.value())
        if(self.ui.comboExpType.currentText()=='SCIENCE'):
            self.settings['lastScienceNexp'] = float(self.ui.displNexp.value())
        elif(self.ui.comboExpType.currentText()=='ARC'):
            self.settings['lastArcNexp'] = float(self.ui.displNexp.value())
        elif(self.ui.comboExpType.currentText()=='FLAT'):
            self.settings['lastFlatNexp'] = float(self.ui.displNexp.value())
        self.settings = settings


    def qlView(self):
        if(1):
            try:
                image = self.con.dc.get_image_data()
                keys = self.con.dc.get_exposure_keys()
                fits.writeto('Img.fits',image,clobber=True)
#                self.con.dc.gen_fits(image, "Img.fits")
###                self.con.gen_fits_with_keys(image, "Img.fits")
                print image.shape
                self.image = image
            except:
#                fits.writeto('tt.fits',self.image,clobber=True)
                return
            # -- test code to try plotting to quick look viewer
#            self.ui.image = np.random.random((255,2048))
            self.ui.qlv.canvas.ax.clear()
            self.ui.qlv.canvas.ax.imshow(image, interpolation='nearest', cmap='gist_gray', \
                                         origin='lower', aspect='auto', vmin=0., vmax=image.max())
            non0 = np.where(image > 0.0)
            inds=np.array(non0)
            self.ui.qlv.canvas.ax.plot(inds[1],inds[0],'or',alpha=0.4)
            self.ui.qlv.canvas.draw()
            self.ui.qlv.show()
            # --

    def quickView(self):
        return #***
        print 'displaying...'
        self.ui.qlv.canvas.ax.clear()
        self.ui.qlv.canvas.ax.imshow(self.image, origin='lower', interpolation='nearest', cmap='gist_gray')
        #self.ui.qlv.canvas.ax.imshow(self.ui.image, interpolation='nearest', cmap='gist_gray')
        self.ui.qlv.show()
        print 'max val = %.2f'%self.image.max()
        print 'min val = %.2f'%self.image.min()

