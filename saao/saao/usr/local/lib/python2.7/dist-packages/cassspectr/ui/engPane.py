import time
from PyQt4 import QtCore, QtGui

import logging as log

from cassspectr.controller.controller import Controller

#from mkClickable import clickable


plc_host = 'localhost'
plc_port = 9090
detector_host = 'localhost'
detector_port = 9091

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s




class EngPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl
        
    def showPLCStatus(self):
        # Write the status of all PLC variables into engineering tab
        
        #plcStatus="fred"

        # Stolen from Carel's CLI:
        status_string = """
======================================================  STATUS =================================================================

FilterwheelPosition   : {FilterwheelPosition}        FilterInit            : {FilterInit}       FilterCentred         : {FilterCentred}
FilterMoving          : {FilterMoving}        FilterFailure         : {FilterFailure}

SlitWidthPosition     : {SlitWidthPosition}        SlitIllumination      : {SlitIllumination}       SlitIlluminationValue : {SlitIlluminationValue}
SlitWidthInitPos      : {SlitWidthInitPos}        SlitWidthInitReq      : {SlitWidthInitReq}       SlitWidthMoving       : {SlitWidthMoving}        SlitWidthFailure      : {SlitWidthFailure}
SlitShutter           : {SlitShutter}        SlitShutterFailure    : {SlitShutterFailure}

GratingAngleSteps     : {GratingAngleSteps:06d}   GratingAngle          : {GratingAngle:03.2f}   GratingAngleInit      : {GratingAngleInit}
AngleInitReq          : {AngleInitReq}        GratingAngleMoving    : {GratingAngleMoving}       GratingAngleLimit1    : {GratingAngleLimit1}        GratingAngleLimit2    : {GratingAngleLimit2}
GratingID             : {GratingID:05d}    GratingInserted       : {GratingInserted}       GratingHatchClosed    : {GratingHatchClosed}        GratingAngleFailure   : {GratingAngleFailure}

CameraFocusInit       : {CameraFocusInit}        CameraFocusMoving     : {CameraFocusMoving}       CameraFocusLimit1     : {CameraFocusLimit1}        CameraFocusLimit2     : {CameraFocusLimit2}
FocusAtPosition       : {FocusAtPosition}        FocusPosition         : {FocusPosition}     FocusLVDT             : {FocusPositionPot}

GMCentred             : {GMCentred}        GMInbeam              : {GMInbeam}       GMMoving              : {GMMoving}        GMFailure             : {GMFailure}

RearOfSlitMirror      : {RearOfSlitMirror}        RoSMirrorFailure      : {RoSMirrorFailure}

ARCMirror             : {ARCMirror}        ARC1                  : {ARC1}       ARC2                  : {ARC2}        ARCMirrorFailure      : {ARCMirrorFailure}

HartmanA              : {HartmanA}        HartmanB              : {HartmanB}       HartmanFailure        : {HartmanFailure}

TopCrateInterlock     : {TopCrateInterlock}        FilterInterlock       : {FilterInterlock}       PneumaticsInterlock   : {PneumaticsInterlock}         ARCInterlock          : {ARCInterlock}
SlitWidthInterlock    : {SlitWidthInterlock}        BottomSignalInterlock : {BottomSignalInterlock}       BottomDriveInterlock  : {BottomDriveInterlock}

================================================================================================================================"""

        status_string_2=""




        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return


        plcStatus = (status_string.format(**status))

        #***
        font = QtGui.QFont("Monospace", 8)#, QtGui.QFont.Bold)
        self.ui.labelPLCstatus.setFont(font)
###        self.ui.labelPLCstatus_2.setFont(font)
        #***
        
        self.ui.labelPLCstatus.setText(QtGui.QApplication.translate("MainWindow", plcStatus, None, QtGui.QApplication.UnicodeUTF8))

###        plcStatus_2 = (status_string_2.format(**status))
        
###        self.ui.labelPLCstatus_2.setText(QtGui.QApplication.translate("MainWindow", plcStatus_2, None, QtGui.QApplication.UnicodeUTF8))

# this has now been moved out of engineering tab to appear on all tabs, but code is kept here for historical reasons


    

    
    def showGraphicView(self):
        
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return


        # should this be done once at start of newcontrol.py??
#        gv = self.ui.graphicsView
        #
#        gv.scene = QtGui.QGraphicsScene(self)
#        scene = gv.scene
#        scene.addPixmap(QtGui.QPixmap("%ssch_bgd.png"%self.pngdir))

        scene = self.ui.scene

        scene.clear()        
        scene.addPixmap(QtGui.QPixmap("%ssch_bgd.png"%self.pngdir))
        scene.addPixmap(QtGui.QPixmap("%sguides.png"%self.pngdir))
        
        if(status['GMInbeam']==True): # out of beam
            scene.addPixmap(QtGui.QPixmap("%sGMacq.png"%self.pngdir))
        else: # in beam
            scene.addPixmap(QtGui.QPixmap("%sGMsci.png"%self.pngdir))
 
           
        if(status['ARCMirror']==False): # out of beam
            scene.addPixmap(QtGui.QPixmap("%sAMout.png"%self.pngdir))
        else:
            scene.addPixmap(QtGui.QPixmap("%sAMin.png"%self.pngdir))
    
        if(status['HartmanA']==True):
            scene.addPixmap(QtGui.QPixmap("%sHM1in.png"%self.pngdir))
        else:
            scene.addPixmap(QtGui.QPixmap("%sHM1out.png"%self.pngdir))
        if(status['HartmanB']==True):
            scene.addPixmap(QtGui.QPixmap("%sHM2in.png"%self.pngdir))
        else:
            scene.addPixmap(QtGui.QPixmap("%sHM2out.png"%self.pngdir))
    
        if(status['RearOfSlitMirror']==True):
            scene.addPixmap(QtGui.QPixmap("%sRoSin.png"%self.pngdir))
        else:
            scene.addPixmap(QtGui.QPixmap("%sRoSout.png"%self.pngdir))
            


        if(status['SlitIllumination']==True):
            scene.addPixmap(QtGui.QPixmap("%sSlitIllum.png"%self.pngdir))
    
#        if(status['SlitShutter']==1): # 1=OPEN
        scene.addPixmap(QtGui.QPixmap("%sSlitShutterOut.png"%self.pngdir))
#        else:
#            scene.addPixmap(QtGui.QPixmap("%sSlitShutterIn.png"%self.pngdir))
    
    
        if(status['GratingHatchClosed']==1):
            scene.addPixmap(QtGui.QPixmap("%sGHclosed.png"%self.pngdir))
        else:
            scene.addPixmap(QtGui.QPixmap("%sGHopen.png"%self.pngdir))
        
        
        # light rays:        
        if(status['GMInbeam']==True):
            scene.addPixmap(QtGui.QPixmap("%srayTop.png"%self.pngdir))
            scene.addPixmap(QtGui.QPixmap("%srayTopRight.png"%self.pngdir))
        
#        if( (status['RearOfSlitMirror']==True) & (status['GMInbeam']==False) &  (status['ARCMirror']==False) & (status['SlitShutter']==1) ):
        if( (status['RearOfSlitMirror']==True) & (status['GMInbeam']==False)  &  (status['ARCMirror']==False) ):
            scene.addPixmap(QtGui.QPixmap("%srayTop.png"%self.pngdir))
            scene.addPixmap(QtGui.QPixmap("%sraysRoS.png"%self.pngdir))
    
        if( (status['RearOfSlitMirror']==False) & (status['GMInbeam']==False) & (status['ARCMirror']==False) & \
            (status['RearOfSlitMirror']==False)  ):
            scene.addPixmap(QtGui.QPixmap("%srayTop.png"%self.pngdir))
            scene.addPixmap(QtGui.QPixmap("%srayMid.png"%self.pngdir))
            scene.addPixmap(QtGui.QPixmap("%sraysToCCD.png"%self.pngdir))
    
        if(status['SlitIllumination']==1):
            scene.addPixmap(QtGui.QPixmap("%sslitIllum.png"%self.pngdir))
            if( (status['GMInbeam']==False) ):
                scene.addPixmap(QtGui.QPixmap("%srayTopRight.png"%self.pngdir))
                scene.addPixmap(QtGui.QPixmap("%srayTopLeft.png"%self.pngdir))
    
    
        if(status['ARC1']==True):
            scene.addPixmap(QtGui.QPixmap("%sarc1.png"%self.pngdir))
            if (status['ARCMirror']==True):
                scene.addPixmap(QtGui.QPixmap("%srayArcScreenToAM.png"%self.pngdir))
                if(status['RearOfSlitMirror']==False):
                    scene.addPixmap(QtGui.QPixmap("%sraysToCCD.png"%self.pngdir))
                else:
                    scene.addPixmap(QtGui.QPixmap("%srayMid.png"%self.pngdir))
                    scene.addPixmap(QtGui.QPixmap("%sraysRoS.png"%self.pngdir))


        
        if(status['ARC2']==True):
            scene.addPixmap(QtGui.QPixmap("%sarc2.png"%self.pngdir))
            if (status['ARCMirror']==True):
                scene.addPixmap(QtGui.QPixmap("%srayArcScreenToAM.png"%self.pngdir))
                if(status['RearOfSlitMirror']==False):
                    scene.addPixmap(QtGui.QPixmap("%sraysToCCD.png"%self.pngdir))
                else:
                    scene.addPixmap(QtGui.QPixmap("%srayMid.png"%self.pngdir))
                    scene.addPixmap(QtGui.QPixmap("%sraysRoS.png"%self.pngdir))




        """
        if(status['ARC1']==True):
            scene.addPixmap(QtGui.QPixmap("%sarc1.png"%self.pngdir))
            if (status['ARCMirror']==True):
                scene.addPixmap(QtGui.QPixmap("%srayArcScreenToAM.png"%self.pngdir))
                scene.addPixmap(QtGui.QPixmap("%sraysToCCD.png"%self.pngdir))

        
        if(status['ARC2']==True):
            scene.addPixmap(QtGui.QPixmap("%sarc2.png"%self.pngdir))
            if (status['ARCMirror']==True):
                scene.addPixmap(QtGui.QPixmap("%srayArcScreenToAM.png"%self.pngdir))
                scene.addPixmap(QtGui.QPixmap("%sraysToCCD.png"%self.pngdir))
        """
        # check for any warnings:

        try:
            checkErrors = self.con.plcc.get_errors()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return


        if len(checkErrors)>0:
#            print checkErrors
            scene.addPixmap(QtGui.QPixmap("%swarning.png"%self.pngdir))
        """
        if( (status['HartmanFailure']==1) | (status['GMFailure']==1) | (status['ARCMirrorFailure']==1) | \
            (status['FilterFailure']==1) | (status['RoSMirrorFailure']==1) | (status['SlitShutterFailure']==1) | \
            (status['GratingAngleFailure']==1) ): #*** is that everything? NO CAMREA FOCUS FAILURE flag???***
            scene.addPixmap(QtGui.QPixmap("%swarning.png"%self.pngdir))
#            scene.addPixmap(QtGui.QPixmap("%scalib.png"%self.pngdir))
        """

        # check for calibration mode:
        #***
        
        if( (status['GMInbeam']==1) | (status['ARCMirror']==1) | (status['RearOfSlitMirror']==1) | \
            (status['HartmanA']==1) | (status['HartmanB']==1) ): #*** lots more to add
            scene.addPixmap(QtGui.QPixmap("%scalib.png"%self.pngdir))
#                scene.addPixmap(QtGui.QPixmap("%scalib.png"%self.pngdir))
        
        # include a 'moving' graphic??


        self.ui.gv.setScene(scene)        
    
        # write labels:
        self.ui.labSchFilterPos.setText(str(status['FilterwheelPosition']))
        self.ui.labSchSlitWidthPos.setText(str(status['SlitWidthPosition']))
        #self.ui.labSchGratingAngle.setText(str(status['GratingAngle']))
        self.ui.labSchGratingAngle.setText('{:03.2f}'.format(status['GratingAngle']))
        self.ui.labSchGratingID.setText(str(status['GratingID']))
        #self.ui.labSchSlitIllumLev.setText(str(status['']))
        self.ui.labSchFocusPos.setText('{:03.3f}'.format(status['FocusPosition']))
    

            
    def showCCDstate(self):
        # these should really be done from dictionary values, but read from GUI for now:
        lastCCDTemp = self.ui.labCCDTemp.text()
        lastCFTemp = self.ui.labColdFingerTemp.text()
        ageofTempinfo = self.polltime - self.lastTempPollTime


        # if no exposure has been taken yet:
        try:
            lastCCDStatus = self.lastCCDStatus
            if lastCCDStatus == 0:
                wordCCDStatus = 'IDLE'
            elif lastCCDStatus == 1:
                wordCCDStatus = 'EXPOSING'
            elif lastCCDStatus == 2:
                wordCCDStatus = 'READING'

        except:
            lastCCDStatus = 'N/A'
            wordCCDStatus = 'N/A'

        try:
            ageofCCDinfo = self.polltime - self.startTime
#            ageofCCDinfo = self.startTime-self.polltime
        except:
            ageofCCDinfo = 0.0

#            print lastCCDStatus, ageofCCDinfo

        status_string = """ ========= DETECTOR STATUS ========

  CCD state = %s  (%s)
  age of last info = %.0f sec

  CCD temp. = %s
  Cold finger temp. = %s
  age of last poll = %.0f sec
"""%(lastCCDStatus, wordCCDStatus, ageofCCDinfo, lastCCDTemp,  lastCFTemp, ageofTempinfo)
        #plcStatus = (status_string.format(**status))

        #***
        font = QtGui.QFont("Monospace", 8)#, QtGui.QFont.Bold)
        self.ui.labCCDstate.setFont(font)
###        self.ui.labelPLCstatus_2.setFont(font)
        #***

#        print ageofTempinfo

#        self.ui.labCCDstate.setText(QtGui.QApplication.translate("MainWindow", status_string, None, QtGui.QApplication.UnicodeUTF8))
        self.ui.labCCDstate.setText(status_string)

