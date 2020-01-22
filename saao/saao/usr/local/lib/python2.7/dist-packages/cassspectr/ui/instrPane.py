import time
from PyQt4 import QtCore, QtGui

import logging as log
import numpy as np
import os

from cassspectr.controller.controller import Controller

from cassspectr.interface.spectrograph_interface.ttypes import PLCThriftException

plc_host = 'localhost'
plc_port = 9090
detector_host = 'localhost'
detector_port = 9091

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

#dbglvl=2

class InstrPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl
        
    def ArcMirrorToggle(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        try:
            if(self.dbglvl>4): print '2: self.con.plcc.get_status()'
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if(self.dbglvl>1):
            print
            print status['ARCMirror']
            #self.logger.debug(status['ARCMirror'])
#            print self.AMstatus
        # for some reason this only works with true/false on some builds rather than 0/1!
        # current status is 0 (OUT of beam) move IN beam
        if(status['ARCMirror']==False):
            # change to in beam
#            self.AMstatus="In"
#            self.ui.labArcMirrorStatus.setText(str("In beam"))
            self.con.plcc.set_arc_mirror(True)
            
        # current status is 1 (IN beam) move OUT of beam
        elif(status['ARCMirror']==1):
            self.con.plcc.set_arc_mirror(False)

        self.logger.info('ARCMirror current status: %s; requested change'%status['ARCMirror'])
            
            
    def slitIllumToggle(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '3: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if(self.dbglvl>1):
            #print status
            print
            print status['SlitIllumination']
        
        if (status['SlitIllumination']==1):  # slit illumated, turn off:
            self.logger.info('slit illumination request: deactivate')
            self.con.plcc.set_slit_illumination_status(0)
        else:
            self.logger.info('slit illumination request: activate')
            self.slitIllumLevel() # set level (in case value has not been changed yet)
            self.con.plcc.set_slit_illumination_status(1)
            
    def slitIllumLevel(self):
        # currently n ofeedback from PLC(!) so just assume desired level is correctly achieved
        lev = int(self.ui.sliderSlitIllum.value())
        if(self.dbglvl>1):
            print 'illumination level %s'%lev
        self.con.plcc.set_slit_illumination_value(lev)
            
            
    def GMToggle(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '4: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if(self.dbglvl>1):
            #print status
            print
            print status['GMInbeam']
            print status['GMCentred']
            print status['GMMoving']
            print status['GMFailure']
        
        
        # I assume GM in beam means ACQ! and GMcentred=science *** CHECK ***
        # -- Currently SCIENCE, toggle to ACQ  
#        if(self.GMstatus=="Sci"):
        if (status['GMCentred']):    
            # change to in-beam[?]/acquisition
#            self.GMstatus="Acq"
            self.con.plcc.set_gm_inbeam()

            self.logger.info('Guide Mirror request: in beam')

#            self.ui.colIndicatorGM_ACQ.setColor(QtGui.QColor(0, 192, 0)) # GREEN
#            self.ui.colIndicatorGM_SCI.setColor(QtGui.QColor(160, 160, 160)) # GREY
            
        # -- Currently ACQ, toggle to SCIENCE 
        elif (status['GMInbeam']):
            #self.GMstatus="Sci"
            
            self.con.plcc.set_gm_centred()
#            self.ui.colIndicatorGM_ACQ.setColor(QtGui.QColor(160, 160, 160)) # GREY
#            self.ui.colIndicatorGM_SCI.setColor(QtGui.QColor(0, 192, 0)) # GREEN
            self.logger.info('Guide Mirror request: centred')

            
        elif (status['GMMoving']):
            #pass
#            self.ui.colIndicatorGM_ACQ.setColor(QtGui.QColor(160, 160, 160)) # GREY
#            self.ui.colIndicatorGM_SCI.setColor(QtGui.QColor(160, 160, 160)) # GREY
            if (self.dbglvl>1): print 'Moving'
            #*** Need to set text indicator
            self.logger.info('Guide Mirror moving')

        
        elif (status['GMFailure']):
            self.logger.info('Guide Mirror FAILURE!!!')
#            pass 
            # deal with warnings! ***
            
        
        
    def slitShutterToggle(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '5: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if(self.dbglvl>1):
            #print status
            print
            print status['SlitShutter']
        # 1=OPEN
        # if CLOSED send to OPEN:

        if (status['SlitShutter']==0): 
            self.con.plcc.set_slit_shutter(1)
        # if OPEN change to CLOSED
        elif (status['SlitShutter']==1):
            self.con.plcc.set_slit_shutter(0)
            
        self.logger.info('Slit shutter current status %s; requested change'%status['SlitShutter'])




    def HM1Toggle(self):
        # Want to set Hartmann status to A(/1) (no longer TOGGLE!)
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return
        if(self.dbglvl>4): print '6: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if (status['HartmanFailure']==1):
            print 'HM failure'
            return
        if(self.dbglvl>1):
            #print status
            print
            print status['HartmanA']
        # if CLOSED send to OPEN:
        if (status['HartmanA']==False): # *** is 0/1 meaning correct?
            self.con.plcc.set_hartman_a(True)
        else:
            self.con.plcc.set_hartman_a(False)

        self.logger.info('Hartmann 1 current status: %s; requested change'%status['HartmanA'])

    def HM2Toggle(self):
        # Want to set Hartmann status to B(/2) (no longer TOGGLE!)
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '7: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if (status['HartmanFailure']==1):
            print 'HM failure'
            return
        if(self.dbglvl>1):
            #print status
            print
            print status['HartmanB']
        # if CLOSED send to OPEN:
        if (status['HartmanB']==False): # *** is 0/1 meaning correct?
            self.con.plcc.set_hartman_b(True)
        else:
            self.con.plcc.set_hartman_b(False)
        self.logger.info('Hartmann 2 current status: %s; requested change'%status['HartmanB'])

    def HMOpen(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '8: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if (status['HartmanFailure']==1):
            print 'HM failure'
            return
        # Set Hartmann shutter to CLEAR/Open by setting both positions to false:
        self.logger.info('Requested Hartmann shutter open')
        if (status['HartmanB']==True): # *** is 0/1 meaning correct?
            self.con.plcc.set_hartman_b(False)
        elif (status['HartmanA']==True): # *** is 0/1 meaning correct?
            self.con.plcc.set_hartman_a(False)
        else:
#            print 'HM already open?'
            pass
            self.logger.info('Hartmann shutter already open?')


    def RoSToggle(self):

        if(self.dbglvl>4): print '9: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if(self.dbglvl>1):
            #print status
            print
            print status['RearOfSlitMirror']
        # if CLOSED send to OPEN:

        if (self.detInterlock):
            if (status['RearOfSlitMirror']==False): # *** is 0/1 meaning correct?
                self.con.plcc.set_rear_of_slit_mirror(True)
                
            elif (status['RearOfSlitMirror']==True): # *** is 0/1 meaning correct?
                self.con.plcc.set_rear_of_slit_mirror(False)
                
            self.logger.info('Rear of Slit Mirror current status: %s; requested change'%status['RearOfSlitMirror'])

    # it seems that Carel's code (or at the PLC level), the arc mirror is automatically set to be put in beam if
    # the arc lamp is turned on.
    #
    # ***
    # We need some sanity checks here to prevent arc lamps being turned on while taking an exposure (since there is no shutter)
    # Maybe also force AM out of beam when lamp is turned off??? but maybe not 
    # ***
    #
    # Hmm, there might be a bug here: 
    # If lamp is left on while arc mirror is taken out, turning lamp off then re-inserts arc mirror into beam.    
        
        
    def ArcLamp1Toggle(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '10: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if(self.dbglvl>1):
            #print status
            print
            print status['ARC1']
        # if CLOSED send to OPEN:
        if (status['ARC1']==False): 
            self.con.plcc.set_arc_lamp_1(True)
        elif (status['ARC1']==True): 
            self.con.plcc.set_arc_lamp_1(False)
        self.logger.info('Arc Lamp (1) current status: %s; requested change'%status['ARC1'])

    def ArcLamp2Toggle(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '11: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        if(self.dbglvl>1):
            #print status
            print
            print status['ARC2']
        # if CLOSED send to OPEN:
        if (status['ARC2']==False): 
            self.con.plcc.set_arc_lamp_2(True)
        elif (status['ARC2']==True): 
            self.con.plcc.set_arc_lamp_2(False)
        self.logger.info('Arc Lamp (2) current status: %s; requested change'%status['ARC1'])
        
        
        
    def changeSlitWidth(self): #*** need to sanity check values and possibly convert to "?
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '12: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

#        reqSWIndx = self.ui.comboBoxSlitWidth.currentIndex() + 1 # 1-28
        reqSWIndx = self.ui.comboBoxSlitWidth.currentIndex()  # 0-27

        if(self.dbglvl>1):
            #print status
            print
            print 'current slit width position %s'%status['SlitWidthPosition']
            
##            print 'requested: %s = %s arcsec'%(reqSWIndx,self.ui.swidths[reqSWIndx-1])
            print 'requested: %s = %s arcsec'%(reqSWIndx,self.ui.swidths[reqSWIndx])
#        if(status['SlitWidthPosition'])
        
        if(self.detInterlock):
            self.con.plcc.set_slit_width(reqSWIndx)
            #*** need to double check this. above print statements were tested with real PLC
            self.logger.info('Slit Width current status: position %s (%s);\nrequested change to position %s (%s)'\
                         %(status['SlitWidthPosition'],self.ui.swidths[status['SlitWidthPosition']],reqSWIndx,self.ui.swidths[reqSWIndx]))
        
        
    def filterGo(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '13: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        reqFiltIndx = self.ui.comboBoxFilter.currentIndex() + 1
        print 'reqFiltIndx %s'%reqFiltIndx
        if(self.dbglvl>1):
            #print status
            print
#            print 'current filter position %s'%status['FilterwheelPosition']
            
        reqFiltKeyword = 'FiltPos%s'%(reqFiltIndx)
        print 'requested: %s = %s'%(reqFiltIndx,self.settings[reqFiltKeyword])
        try:
            currFiltKeyword = 'FiltPos%s'%(status['FilterwheelPosition'])
        except:
            pass
#                currFiltKeyword = 'FiltPos1' #*** workaround for moment!
##            if (status['FilterwheelPosition']==0): 
##                # do sthg***!
##                currFiltKeyword = 'FiltPos1' #??

        self.con.plcc.set_fw_move(reqFiltIndx)
        self.logger.info('Filter wheel current status: position %s (%s);\nrequested change to position %s (%s)'\
                         %(status['FilterwheelPosition'],self.settings[currFiltKeyword],reqFiltIndx,self.settings[reqFiltKeyword]))

    def autoFilterInsert(self,status):
        # For certain grating angle and grating combos, we want to put in a filter:
        # gr5, angle<0: GG 495
        # gr5, angle>0: BG 39
        # gr4,6,7: OPEN
        # gr8 1st order GG 495
        # gr9 1st order GG 495
        # gr8,9 2nd order BG 39
        # gr10,11,12 GG 495

        gratingAngle = status['GratingAngle']
        grName = self.settings['grName']
        if ( (grName=='gr4') | (grName=='gr6') | (grName=='gr7') ):
            reqFilter = 'OPEN'
        elif (grName=='gr5'):
            if(gratingAngle<0.):
                reqFilter = 'GG 495'
            else:
                reqFilter = 'BG 39'
        elif ( (grName=='gr10') | (grName=='gr11') | (grName=='gr12') ):
                reqFilter = 'GG 495'
        else:
            self.logger.info('auto filter not implemented for %s yet'%grName)
        # insert required filter:
        reqFiltIndx = self.ui.comboBoxFilter.findText(reqFilter)
        self.logger.info('auto-inserting filter %s for this setup'%reqFilter)
        self.con.plcc.set_fw_move(reqFiltIndx+1)


    def gratingAngleGo(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '14: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        currGratingAngle = status['GratingAngle']
        # need to build angle in decimal degs from degs and mins
        # Careful of sign!
        
        reqGratingAngle = float(self.ui.lineEditGratingAngleDeg.text())
        self.logger.info('Current grating angle %sdeg, requested %sdeg'%(currGratingAngle,reqGratingAngle))

        try: 
            self.con.plcc.grating_angle_move_abs(reqGratingAngle)
        except PLCThriftException, e:
            self.logger.warning("Exception thrown by set_grating_angle():\n {}".format(e.message))
            print e.message

        # Use this step to check if user needs a filter 
        try:
            self.autoFilterInsert(status)
        except Exception as e:
            print e
            print 'could not auto insert filter'

    """
    def gratingAngleIncrementGo(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '15: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        currGratingAngle = status['GratingAngle']
        
        reqGratingAngleIncrement = float(self.ui.lineEditGratingAngleIncrementDeg.text())
        self.logger.info('Current grating angle %sdeg, requested increment %sdeg'%(currGratingAngle,reqGratingAngleIncrement))

        try: 
            self.con.plcc.grating_angle_move_rel(reqGratingAngleIncrement)
        except PLCThriftException, e:
            self.logger.warning("Exception thrown by grating_angle_rel():\n {}".format(e.message))
            print e.message
    """





    def cameraFocusGo(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '16: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

#        currFocusPos = status['FocusAtPosition']
        currFocusPos = status['FocusPosition']
        reqFocusPos = float(self.ui.lineEditFocus.text())
        if (self.dbglvl>1):
            print reqFocusPos
        # Hmm, don't understand Carel's syntax
#        intval = 39999 ; fracval = reqFocusPos/float(intval)
        value = reqFocusPos
#        intval = int(np.floor(value))
#        fracval = 100 * (value - intval)
        
        self.logger.info('Current focus position %s, requested %s'%(currFocusPos,reqFocusPos))
        try:
            #self.con.plcc.set_camera_focus_value(reqFocusPos)
            #self.con.plcc.set_camera_focus_value(reqFocusPos)
            #self.con.plcc.camera_focus_move(reqFocusPos)
            self.con.plcc.camera_focus_move_abs(value)
        except PLCThriftException, e:
            self.logger.warning("Exception thrown by camera_focus_move_abs():\n {}".format(e.message))



    """
    def cameraFocusIncrementGo(self):
        if(self.detInterlock==0): 
            if(self.dbglvl>1):
                print 'locked!'
            return

        if(self.dbglvl>4): print '17: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return

        currFocusPos = status['FocusPosition']
        reqFocusIncrementPos = float(self.ui.lineEditFocusIncrement.text())
        if (self.dbglvl>1):
            print reqFocusIncrementPos
        value = reqFocusIncrementPos
        
        self.logger.info('Current focus position %s, requested increment %s'%(currFocusPos,reqFocusIncrementPos))
        try:
            self.con.plcc.camera_focus_move_rel(value)
        except PLCThriftException, e:
            self.logger.warning("Exception thrown by camera_focus_move_rel():\n {}".format(e.message))
    """


    def setupArc(self):
        # set up for taking an arc, once exptype is selected:
        # should be no harm in setting these changes even if already at requested setting
        # Piet's PLC logic will take care of this
        self.con.plcc.set_arc_mirror(True)
        # Need to check which lamp we were using before (almost certainly 2):
        #        if(self.settings['
        # Actually, just assume 2 (CuAr) for everything except gr5 and gr-angle<0. (CuNe):
        if ( (self.ui.labGratingName.text()=='gr5') & (float(self.ui.labGratingAngleCurrent.text())<0.) ):
            self.con.plcc.set_arc_lamp_2(False)
            self.con.plcc.set_arc_lamp_1(True) # CuNe
        else:
            self.con.plcc.set_arc_lamp_2(True) # CuAr
            self.con.plcc.set_arc_lamp_1(False)

        # also, let's make sure RoS, GM, Hartmann are all in correct posns:
        self.con.plcc.set_rear_of_slit_mirror(False)
        self.con.plcc.set_gm_centred()
        self.con.plcc.set_hartman_a(False)
        self.con.plcc.set_hartman_b(False)

    def unsetArc(self):
        self.con.plcc.set_arc_mirror(False)
        # Need to check which lamp we were using before (almost certainly 2):
        #        if(self.settings['
        # Actually, just assume 2:
        self.con.plcc.set_arc_lamp_2(False)
        self.con.plcc.set_arc_lamp_1(False)
        # also, let's make sure RoS, GM, Hartmann are all in correct posns:
        self.con.plcc.set_rear_of_slit_mirror(False)
        self.con.plcc.set_gm_centred()
        self.con.plcc.set_hartman_a(False)
        self.con.plcc.set_hartman_b(False)


    def setupHartmann(self):
        # As per Lisa's request, if user selects Hartmann, setup for arc with:
        # "It'd be good if opting to run a Hartmann sequence would put the arc mirror in the beam,
        # turn on the relevant lamp & set the exposure time (10s for CuAr & 1s for CuNe) -& remove the rear-of-slit mirror if it happens to be in the beam."
        #
        # DGG: did we also want to set the slit width to be wide/narrow? I thought we'd done this somewhere?!?
        self.setupArc()
        # for CuAr, set exptime = 10s
        # for CuNe, set exptime = 1s
        # taken from abov. Need to make sure consistent if setupArc() changes!
        if ( (self.ui.labGratingName.text()=='gr5') & (float(self.ui.labGratingAngleCurrent.text())<0.) ):
            # CuNe:
            self.ui.displExpTime.setText('1.')
        else: # CuAr
            self.ui.displExpTime.setText('10.0')




    def translateHartmannPos(self,status):
        hartState = ''
        if(status['HartmanFailure']==True):
            hartState='FAIL'
        else:
            if (status['HartmanA']==True):
                hartState='A'
            elif (status['HartmanB']==True):
                hartState='B'
            elif ( (status['HartmanA']==False) & (status['HartmanB']==False) ):
                hartState='Open'
        return hartState



    def changeFocusAndHartmannandWait(self,requestedFocusPos,requestedHartmannPos):
        if(self.dbglvl>4): print '18: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except Exception as e:
            self.logger.warn("Caught exception: {}".format(e))
            return
        hartState = self.translateHartmannPos(status)    
        focusPos = status['FocusPosition']
        # send movement commands to both:
#        self.focus***
        try:
            self.con.plcc.camera_focus_move_abs(requestedFocusPos)
        except PLCThriftException, e:
            self.logger.warning("Exception thrown by camera_focus_move_abs():\n {}".format(e.message))

        # hart:
        if(requestedHartmannPos=='A'):
            self.con.plcc.set_hartman_a(True)
        elif(requestedHartmannPos=='B'):
            self.con.plcc.set_hartman_b(True)
        else:
            print 'do not understand hartmann pos: %s'%requestedHartmannPos
        while ( (hartState<>requestedHartmannPos) | (np.abs(focusPos-requestedFocusPos)>0.003) ):
            # stay here till we're at requested positions:
            if(self.dbglvl>4): print '19: self.con.plcc.get_status()'
            try:
                status = self.con.plcc.get_status()
            except Exception as e:
                self.logger.warn("Caught exception: {}".format(e))
                return
            hartState = self.translateHartmannPos(status)    
            focusPos = status['FocusPosition']
        print 'READY!',focusPos
        self.drawInstrIndicators()
        self.showGraphicView()
        self.showPLCStatus()
        self.showGraphicView()
        self.showTargets()
        self.showTime()
            
        # take exposure here and wait until it finishes:
        self.drawInstrIndicators()
        self.showGraphicView()

        self.HartmannExposureWaiting=True
        self.doFirstExposure()
        # wait for exposure to finish
        while(self.HartmannExposureWaiting):
            print 'waiting'

            #stay here
            #self.drawExp()
            self.drawDetIndicators()
            self.drawInstrIndicators()
            self.showGraphicView()
            self.showPLCStatus()
            self.showGraphicView()
            self.showTargets()
            self.showTime()

            ## save GUI settings constantly
            self.readConfig()



        

    def HartmannSequence(self):

        # setup hartmann:
        self.setupHartmann()

        # read values from GUI and run sequence:

        # first, error check start and end focus points:
        # *** TODO
        self.logger.info('Starting Hartmann sequence')

        # switch to 1" slit for routine and switch back after:
        lastSlitText = self.ui.labSlitWidthCurrent.text() # 10: 1.50"
        self.lastSlitWidth = int(lastSlitText.split(":")[0])
        #self.ui.comboBoxSlitWidth.setCurrentIndex(7) # 1.05"
        self.con.plcc.set_slit_width(7)

        nsteps = int(self.ui.spinBoxHartmannNsteps.value())
        focus0 = float(self.ui.lineEditHartmannInitialFocusPos.text())
        focusStep = float(self.ui.lineEditHartmannFocusIncrement.text())
        print 'Hartmann sequence...'
        print nsteps
        print focus0
        lastExpType = self.settings['expType']
        lastNExp = self.settings['Nexp']
        # set to Hartmann:
        self.settings['expType']='HARTMANN'
        # Need to force update:
        #        self.drawDetIndicators()
        indx = self.ui.comboExpType.findText('HARTMANN')
        self.ui.comboExpType.setCurrentIndex(indx)

        self.settings['Nexp']=1 # make sure to only take 1 exposure at each position

        self.settings['HartmannSeq'] = '0/0'
        for i in range(nsteps):
            self.logger.debug('in Hartmann loop, i=%s/%s\n'%(i,nsteps))
            self.settings['HartmannSeq'] = '%s/%s'%(i+1,nsteps) # 1-based numbering
            tfocus = focus0 + float(i)*focusStep
            
            print 'change focus here'
            if i%2==0:
                print tfocus, i,  self.settings['HartmannSeq'], 'B'                
                self.changeFocusAndHartmannandWait(tfocus,'B')
                print tfocus, i,  self.settings['HartmannSeq'], 'A'                
                self.changeFocusAndHartmannandWait(tfocus,'A')

            else:
                print tfocus, i,  self.settings['HartmannSeq'], 'A'                
                self.changeFocusAndHartmannandWait(tfocus,'A')
                print tfocus, i,  self.settings['HartmannSeq'], 'B'                
                self.changeFocusAndHartmannandWait(tfocus,'B')

        # reset:
        self.settings['HartmannSeq'] = '0/0'
        self.logger.debug('last slitwidth = %.3f'%self.lastSlitWidth)
        self.con.plcc.set_slit_width(self.lastSlitWidth)
            



        self.logger.debug('Hartmann sequence finished successfully. Resetting values to previous states')

        # return values to previous:
        self.settings['expType']=lastExpType
        self.settings['Nexp']=lastNExp
        self.con.plcc.set_hartman_a(False)
        self.con.plcc.set_hartman_b(False)

        self.logger.info('Hartmann sequence completed')

        # -- look for best-fit focus written by qlspupnic:
        try:
            requestedFocusPos = np.loadtxt('/tmp/tmp_focus.dat')
            self.logger.info('Focus position from ql-spupnic: %.3f'%requestedFocusPos)
            try:
                self.con.plcc.camera_focus_move_abs(requestedFocusPos)
            except PLCThriftException, e:
                self.logger.warning("Exception thrown by camera_focus_move_abs():\n {}".format(e.message))

            try:
                status = self.con.plcc.get_status()
            except Exception as e:
                self.logger.warn("Caught exception: {}".format(e))
                return
            focusPos = status['FocusPosition']
            # delete focus file so we don't accidnetally pick up old values next time:
            os.system('rm -f /tmp/tmp_focus.dat')
        except:
            self.logger.info('Could not find focus results from ql-spupnic')


    def scienceSequence(self):
        # set back to Main CCD tab:
        self.ui.CCDtabWidget.setCurrentIndex(0)


        # take a sequence specified in sequence tab. Likely to be a blocking call like HartmannSequence:
        Narc1 = self.ui.displNexpSeqArc1.value()
        Nloops = self.ui.displSeqNLoops.value()
        Nsci = self.ui.displNexpSeqSci.value()
        Narc2 = self.ui.displNexpSeqArc2.value()
        exptimeArc1 = float(self.ui.displExpTimeSeqArc1.text())
        exptimeSci = float(self.ui.displExpTimeSeqSci.text())
        exptimeArc2 = float(self.ui.displExpTimeSeqArc2.text())

        #self.ui.displNexp.setText(str(self.Narc1exp))
        # do multiple exposures as explicit loops here (rather than by setting Nexp explicitly)
        # then we can set self.HartmannWaiting=True at start of each exposure and use it to see when we've finished.

        # -- confirm sequence:
        nInSeqTot = Narc1 + Nloops*(Nsci+Narc2)
        print 'sequence requested'
        self.logger.info('sequence requested')
        for iii in range(Narc1):
            print '%ss arc'%exptimeArc1
            self.logger.info('%ss arc'%exptimeArc1)
        for iii in range(Nloops):
            for jjj in range(Nsci):
                print '%ss arc'%exptimeSci
                self.logger.info('%ss science'%exptimeSci)
            for jjj in range(Narc2):
                print '%ss arc'%exptimeArc2
                self.logger.info('%ss arc'%exptimeArc2)

        nInSeq = 0
        for iii in range(Narc1):
            nInSeq +=1
            print 'step %s / %s in sequence'%(nInSeq,nInSeqTot)
            self.logger.info('step %s / %s in sequence'%(nInSeq,nInSeqTot))
            # do initial arc:
            self.setupArc()
            # set exposure type, time and Narc1exp:
            self.settings['expType']='ARC'
            reqExpTypeIndx = self.ui.comboExpType.findText('ARC')
            self.ui.comboExpType.setCurrentIndex(reqExpTypeIndx)
            self.settings['expTime'] = exptimeArc1
            self.ui.displExpTime.setText(str(self.settings['expTime']))
            time.sleep(5) # wait for setup to complete [crudely]

            self.HartmannExposureWaiting=True
            self.doFirstExposure()
            # wait for exposure to finish
            while(self.HartmannExposureWaiting):
                print 'waiting'

                #stay here
                #self.drawExp()
                self.drawDetIndicators()
                self.drawInstrIndicators()
                self.showGraphicView()
                self.showPLCStatus()
                self.showGraphicView()
                self.showTargets()
                self.showTime()

                ## save GUI settings constantly
                self.readConfig()
        # switch to science (and arc) loops:
        for iii in range(Nloops):
            # science frame(s)
            for jjj in range(Nsci):
                nInSeq +=1
                print 'step %s / %s in sequence'%(nInSeq,nInSeqTot)
                self.logger.info('step %s / %s in sequence'%(nInSeq,nInSeqTot))
                self.logger.info('science %ss'%exptimeSci)

                self.unsetArc()
                self.settings['expType']='SCIENCE'
                reqExpTypeIndx = self.ui.comboExpType.findText('SCIENCE')
                self.ui.comboExpType.setCurrentIndex(reqExpTypeIndx)
                self.settings['expTime'] = exptimeSci
                self.ui.displExpTime.setText(str(self.settings['expTime']))
                time.sleep(5) # wait for setup to complete [crudely]

                self.HartmannExposureWaiting=True
                self.doFirstExposure()
                # wait for exposure to finish
                while(self.HartmannExposureWaiting):
                    print 'waiting'

                    #stay here
                    #self.drawExp()
                    self.drawDetIndicators()
                    self.drawInstrIndicators()
                    self.showGraphicView()
                    self.showPLCStatus()
                    self.showGraphicView()
                    self.showTargets()
                    self.showTime()

                    ## save GUI settings constantly
                    self.readConfig()

            for jjj in range(Narc2):
                nInSeq +=1
                print 'step %s / %s in sequence'%(nInSeq,nInSeqTot)
                self.logger.info('step %s / %s in sequence'%(nInSeq,nInSeqTot))

                self.setupArc()
                self.settings['expType']='ARC'
                reqExpTypeIndx = self.ui.comboExpType.findText('ARC')
                self.ui.comboExpType.setCurrentIndex(reqExpTypeIndx)
                self.settings['expTime'] = exptimeArc2
                self.ui.displExpTime.setText(str(self.settings['expTime']))
                time.sleep(5) # wait for setup to complete [crudely]

                self.HartmannExposureWaiting=True
                self.doFirstExposure()
                # wait for exposure to finish
                while(self.HartmannExposureWaiting):
                    print 'waiting'

                    #stay here
                    #self.drawExp()
                    self.drawDetIndicators()
                    self.drawInstrIndicators()
                    self.showGraphicView()
                    self.showPLCStatus()
                    self.showGraphicView()
                    self.showTargets()
                    self.showTime()

                    ## save GUI settings constantly
                    self.readConfig()


        print 'Sequence done.'
        self.unsetArc()

