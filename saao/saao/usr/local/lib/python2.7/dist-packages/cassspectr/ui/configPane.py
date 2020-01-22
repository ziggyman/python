import time, datetime
from PyQt4 import QtCore, QtGui

import logging as log

from cassspectr.controller.controller import Controller
import numpy as np

import os,json


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

#dbglvl=2




class ConfigPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl


    def createDefaultConfig(self):
        # check if configdir exists:
        if not os.path.isdir(self.configdir):
            os.mkdir(self.configdir)

        settings = {"FiltPos5": "5", "FiltPos4": "4", "FiltPos7": "7", "FiltPos6": "6", "FiltPos1": "OPEN", "observer": "Smith", \
                    "FiltPos2": "2", "expType": "FLAT", "expName": "", "Nexp": "1", "FiltPos8": "8", "lastArcNexp": 1.0, "lastFlatNexp": 1.0, \
                    "lastArcExpTime": 2.0, "FiltPos3": "3", "soundEndReadout": True, "lastScienceNexp": 6.0, "lastFlatExpTime": 5.0, \
                    "soundEndExposure": False, "expTime": "5.0", "lastBiasnExp": 11.0, "ArcLamp1ID": "CuAr", "ArcLamp2ID": "CuNe", \
                    "propID": "001", "lastBiasNexp": 11.0, "ccdBinning": "1x1", "comment": "", "frameNum": "1001", \
                    "lastScienceExpTime": 1200.0}

        if(1): # add settings previously in ccdpars:
            ccdpars = {}
            ccdpars['nRows'] = 256
            ccdpars['nCols'] = 2148
            binString = settings['ccdBinning']
            xbin = binString.split('x')[0]
            ybin = binString.split('x')[1]
            if (self.dbglvl>2):
                print 'binning: %s x %s'%(xbin,ybin)
            ccdpars['nColBin'] = int(xbin)
            ccdpars['nRowBin'] = int(ybin)

#            ccdpars['nRowBin'] = 1
#            ccdpars['nColBin'] = 1
            ccdpars['nRowCen'] = 384
            ccdpars['nColCen'] = 1124
#            ccdpars['gainMode'] = 'TODO'
            ccdpars['CCDmodeGainIndx'] = 0
            ccdpars['CCDmodeNoiseIndx'] = 0
#            ccdpars['exptime'] = 10.0
#            ccdpars['frameNum'] = 10001 # *** should probably check for existing images to make sure not over-written!
#            ccdpars['EXPTYPE'] = 'SCIENCE'

        for k,v in ccdpars.items():
            settings[k] = v

        self.settings = settings
        self.writeConfigFile()


    def readConfigFile(self):
        if not os.path.isfile('%ssettings.json'%self.configdir):
            if(self.dbglvl>1):
                print 'settings.json file not found. creating'
            self.createDefaultConfig()
#        settings = np.load('settings.npy')
        f=open('%ssettings.json'%self.configdir,'r') #*** change path to home dir or somewhere later? ***
        settings=json.load(f)
        f.close()
        self.settings = settings

    def writeConfigFile(self):
        #np.save('%ssettings.npy'%self.configdir,self.settings)
        f=open('%ssettings.json'%self.configdir,'w')
        json.dump(self.settings,f)
        f.close()
        if self.dbglvl>3: print 'written settings.json'
#        TypeError: PyQt4.QtCore.QString(u'5') is not JSON serializable
#        TypeError: closeEvent() takes exactly 1 argument (2 given)


    def readConfig(self):
        # read through variables in config panel and populate relevant tabs, boxes, etc.
        settings = self.settings
#        settings = {} # don't overwite. Always read from config file at newcontrol.py startup!
        settings['ArcLamp1ID'] = str(self.ui.lineEditArc1ID.text())
        settings['ArcLamp2ID'] = str(self.ui.lineEditArc2ID.text())
        settings['FiltPos1'] = str(self.ui.lineEditFiltPos1ID.text())
        settings['FiltPos2'] = str(self.ui.lineEditFiltPos2ID.text())
        settings['FiltPos3'] = str(self.ui.lineEditFiltPos3ID.text())
        settings['FiltPos4'] = str(self.ui.lineEditFiltPos4ID.text())
        settings['FiltPos5'] = str(self.ui.lineEditFiltPos5ID.text())
        settings['FiltPos6'] = str(self.ui.lineEditFiltPos6ID.text())
        settings['FiltPos7'] = str(self.ui.lineEditFiltPos7ID.text())
        settings['FiltPos8'] = str(self.ui.lineEditFiltPos8ID.text())


        # also read relevant GUI settings from main GUI (to store last-used state):
        settings['frameNum'] = str(self.ui.displFrameNum.text())
        settings['ccdBinning'] = str(self.ui.comboCCDbinning.currentText())
        settings['expType'] = str(self.ui.comboExpType.currentText())
        settings['expTime'] = str(self.ui.displExpTime.text())
        settings['Nexp'] = str(self.ui.displNexp.value())
        # header:
        settings['observer'] = str(self.ui.lineEditObserver.text())
        settings['propID'] = str(self.ui.lineEditPropID.text())
        settings['comment'] = str(self.ui.lineEditComment.text())

        settings['soundEndExposure'] = (self.ui.checkBoxSoundEndExposure.isChecked())
        settings['soundEndReadout'] = (self.ui.checkBoxSoundEndReadout.isChecked())


        self.settings = settings
        
        
    def applyConfig(self):
        self.ui.labArcLamp1.setText('Arc 1: %s'%self.settings['ArcLamp1ID'])
        self.ui.labArcLamp2.setText('Arc 2: %s'%self.settings['ArcLamp2ID'])
        filtList = [self.settings['FiltPos1'],\
                    self.settings['FiltPos2'],\
                    self.settings['FiltPos3'],\
                    self.settings['FiltPos4'],\
                    self.settings['FiltPos5'],\
                    self.settings['FiltPos6'],\
                    self.settings['FiltPos7'],\
                    self.settings['FiltPos8']]
        
        self.ui.comboBoxFilter.clear()            
        self.ui.comboBoxFilter.addItems(filtList)
        
        # apply to config tab too!
        self.ui.lineEditArc1ID.setText(self.settings['ArcLamp1ID'])
        self.ui.lineEditArc2ID.setText(self.settings['ArcLamp2ID'])
#        self.ui.lineEditManFiltID.setText()
        self.ui.lineEditFiltPos1ID.setText(self.settings['FiltPos1'])
        self.ui.lineEditFiltPos2ID.setText(self.settings['FiltPos2'])
        self.ui.lineEditFiltPos3ID.setText(self.settings['FiltPos3'])
        self.ui.lineEditFiltPos4ID.setText(self.settings['FiltPos4'])
        self.ui.lineEditFiltPos5ID.setText(self.settings['FiltPos5'])
        self.ui.lineEditFiltPos6ID.setText(self.settings['FiltPos6'])
        self.ui.lineEditFiltPos7ID.setText(self.settings['FiltPos7'])
        self.ui.lineEditFiltPos8ID.setText(self.settings['FiltPos8'])


#        try:
        # also relevant GUI settings from main GUI (to store last-used state):
        self.ui.displFrameNum.setText(str(self.settings['frameNum']))
        index = self.ui.comboCCDbinning.findText(str(self.settings['ccdBinning']))
        self.ui.comboCCDbinning.setCurrentIndex(index)
        index = self.ui.comboExpType.findText(str(self.settings['expType']))
        self.ui.comboExpType.setCurrentIndex(index)
        #self.ui.comboExpType.setText(self.settings['expType'])
        self.ui.displExpTime.setText(self.settings['expTime'])
        self.ui.displNexp.setValue(int(self.settings['Nexp']))
        # header:
        self.ui.lineEditObserver.setText(self.settings['observer'])
        self.ui.lineEditPropID.setText(self.settings['propID'])
        self.ui.lineEditComment.setText(self.settings['comment'])

        self.ui.comboCCDmodeGain.setCurrentIndex(int(self.settings['CCDmodeGainIndx']))
        self.ui.comboCCDmodeNoise.setCurrentIndex(int(self.settings['CCDmodeNoiseIndx']))


        self.ui.checkBoxSoundEndExposure.setCheckState(self.settings['soundEndExposure'])
        self.ui.checkBoxSoundEndReadout.setCheckState(self.settings['soundEndReadout'])

#        except:
#            pass

    def readLastDataDir(self):
        f = open('%sdatadir.json'%self.configdir,'r')
        self.datadir = json.load(f)
        print self.datadir
        f.close()



    def OLDautoDataDir(self):
        # create a datadir automatically for tonight's data. Switch over at 10am (probably a reasonable time for noone to be using GUI)

        # current hour:
        lt = time.localtime()
        if (time.strftime("%H", lt) < 10.):
            # get yesterday's date.
            ld = datetime.date(lt.tm_year, lt.tm_mon, lt.tm_mday)-datetime.timedelta(1.)
        else:
            ld = datetime.date(lt.tm_year, lt.tm_mon, lt.tm_mday)

        fmtdate = ld.strftime("%Y%m%d")
        self.datadir = self.rootdatadir+fmtdate+'/'
        f=open('%sdatadir.json'%self.configdir,'w')
        json.dump(self.datadir,f)
        f.close()

        self.ui.lineEditDataDir.setText(self.datadir)

        if not os.path.isdir(self.datadir):
            os.mkdir(self.datadir)
            #*** not got this bit working:
            self.settings['frameNum']=1001
#            self.settings['datadir'] = self.datadir
            self.writeConfigFile()

    def autoDataDir(self):
        # create a datadir automatically for tonight's data. Switch over at 10am (probably a reasonable time for noone to be using GUI)

#        Now just find the name of what it should be and populate the lineEdit under Advanced

        # current hour:
        lt = time.localtime()
        if (time.strftime("%H", lt) < 10.):
            # get yesterday's date.
            ld = datetime.date(lt.tm_year, lt.tm_mon, lt.tm_mday)-datetime.timedelta(1.)
        else:
            ld = datetime.date(lt.tm_year, lt.tm_mon, lt.tm_mday)

        fmtdate = ld.strftime("%Y%m%d")
        #self.
        Autodatadir = self.rootdatadir+fmtdate+'/'
        #f=open('%sdatadir.json'%self.configdir,'w')
        #json.dump(self.datadir,f)
        #f.close()

        self.ui.lineEditDataDir.setText(Autodatadir)

        #if not os.path.isdir(self.datadir):
        #    os.mkdir(self.datadir)
        #    #*** not got this bit working:
        #    self.settings['frameNum']=1001
#       #     self.settings['datadir'] = self.datadir
        #    self.writeConfigFile()


    def manualDataDir(self):
        # read from lineEdit in Advanced tab
        self.datadir = str(self.ui.lineEditDataDir.text())
        #--
        # check for trailing slash:
        if self.datadir[-1] <> '/':
            self.datadir = self.datadir+'/'
            self.ui.lineEditDataDir.setText(self.datadir)
        #--
        f=open('%sdatadir.json'%self.configdir,'w')
        json.dump(self.datadir,f)
        f.close()
        self.settings['frameNum']=1001
        self.writeConfigFile()
        if not os.path.isdir(self.datadir):
            os.mkdir(self.datadir)
#            #*** not got this bit working:
#            self.settings['frameNum']=1001
##            self.settings['datadir'] = self.datadir
            self.writeConfigFile()

#    def setStartFrameNumber(self):
#        self.ui.displFrameNum.setText(str(self.ui.lineEditStartFrameNumber.text()))


    # def setBrake(self):
    #     if self.gratingBrakeOn:
    #         # Brake is on. Turn it off and enable grating angle changes
    #         self.gratingBrakeOn = False
    #         self.ui.brakeButton.setText("Apply Grating Brake")
    #         self.ui.pushButGratingAngleGo.setEnabled(True)
    #         self.ui.lineEditGratingAngleDeg.setEnabled(True)
    #         self.ui.pushButGratingAngleInit.setEnabled(True)
    #         self.ui.pushButInitAll.setEnabled(True)
    #     else:
    #         # Brake is off. Turn it on and disable grating angle changes 
    #         self.gratingBrakeOn = True
    #         self.ui.brakeButton.setText("Deactivate Grating Brake")
    #         self.ui.pushButGratingAngleGo.setEnabled(False)
    #         self.ui.lineEditGratingAngleDeg.setEnabled(False)
    #         self.ui.pushButGratingAngleInit.setEnabled(False)
    #         self.ui.pushButInitAll.setEnabled(False)
    #         self.ui.pushButResetPCIController.setEnabled(False)

    def unlockBox(self):
        if (self.configLocked):
            self.showDialog()
        else:
            self.doLockPanel()
        """
        self.btn = QtGui.QPushButton('Dialog', self)
        self.btn.move(20, 20)
        self.btn.clicked.connect(self.showDialog)
        
        self.le = QtGui.QLineEdit(self)
        self.le.move(130, 22)
        
        self.setGeometry(300, 300, 290, 150)
        self.setWindowTitle('Input dialog')
        self.show()
        """
        
    def doUnlockPanel(self):
        # unlock all:
        self.ui.lineEditManGratingName.setEnabled(True)
        self.ui.lineEditFiltPos1ID.setEnabled(True)
        self.ui.lineEditFiltPos2ID.setEnabled(True)
        self.ui.lineEditFiltPos3ID.setEnabled(True)
        self.ui.lineEditFiltPos4ID.setEnabled(True)
        self.ui.lineEditFiltPos5ID.setEnabled(True)
        self.ui.lineEditFiltPos6ID.setEnabled(True)
        self.ui.lineEditFiltPos7ID.setEnabled(True)
        self.ui.lineEditFiltPos8ID.setEnabled(True)
        self.ui.lineEditArc1ID.setEnabled(True)
        self.ui.lineEditArc2ID.setEnabled(True)
        #        self.ui.brakeButton.setEnabled(True)
        # self.ui.pushButResetPCIController.setEnabled(True)
        # self.ui.pushButResetDetectorController.setEnabled(True)
        # self.ui.pushButSlitWidthReset.setEnabled(True)
        # self.ui.pushButGratingAngleReset.setEnabled(True)
        # self.ui.pushButFWReset.setEnabled(True)
        # self.ui.pushButCameraFocusReset.setEnabled(True)
        # self.ui.pushButResetAll.setEnabled(True)
        # set button name to 'Lock'
        self.ui.pushButUnlockConfig.setText('Lock')
        self.configLocked=False
        
    def doLockPanel(self):
        # update settings:
        
        self.readConfig()
        self.applyConfig()
        self.writeConfigFile()
        
        # unlock all:
        self.ui.lineEditManGratingName.setEnabled(False)
        self.ui.lineEditFiltPos1ID.setEnabled(False)
        self.ui.lineEditFiltPos2ID.setEnabled(False)
        self.ui.lineEditFiltPos3ID.setEnabled(False)
        self.ui.lineEditFiltPos4ID.setEnabled(False)
        self.ui.lineEditFiltPos5ID.setEnabled(False)
        self.ui.lineEditFiltPos6ID.setEnabled(False)
        self.ui.lineEditFiltPos7ID.setEnabled(False)
        self.ui.lineEditFiltPos8ID.setEnabled(False)
        self.ui.lineEditArc1ID.setEnabled(False)
        self.ui.lineEditArc2ID.setEnabled(False)
        #        self.ui.brakeButton.setEnabled(False)
        # self.ui.pushButResetDetectorController.setEnabled(False)
        # self.ui.pushButSlitWidthReset.setEnabled(False)
        # self.ui.pushButGratingAngleReset.setEnabled(False)
        # self.ui.pushButFWReset.setEnabled(False)
        # self.ui.pushButCameraFocusReset.setEnabled(False)
        # self.ui.pushButResetAll.setEnabled(False)
        # set button name to 'Lock'
        self.ui.pushButUnlockConfig.setText('Unlock')
        self.configLocked=True
        
        
        
        
    def showDialog(self):
        
        text, ok = QtGui.QInputDialog.getText(self, 'Input Dialog', 
            'Enter password to unlock config panel:')
        
        # check if already unlocked:
        if(not self.configLocked):
            self.doLockPanel()
#            pass
        
        if (text=='SpUP'):
            #print 'accepted'
            self.doUnlockPanel()
        #else:
        #    print 'rejected'
        
    
    
    
    """
    def closeEvent(self):

        quit_msg = "Are you sure you want to quit?"
        #reply = QtGui.QMessageBox.question(self, 'Message', quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        reply = QtGui.QMessageBox.question(None, 'Removal', quit_msg, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if (reply == QtGui.QMessageBox.Yes):
#            QtGui.QApplication.quit() # Nope!
            QtGui.qApp.closeAllWindows()
        else:
            pass
    """
