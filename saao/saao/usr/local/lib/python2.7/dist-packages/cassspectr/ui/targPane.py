import time
from PyQt4 import QtCore, QtGui

import logging as log

from cassspectr.controller.controller import Controller

from cassspectr.utils.tcs_pull import read_remote_file, get_last_line, get_tcs_info, write_remote_file

plc_host = 'localhost'
plc_port = 9090
detector_host = 'localhost'
detector_port = 9091

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


class TargPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl

    """
    def showTargets(self):
        # update target coords in target list
        if(self.targetsLoaded == False):
            return
        currTargInd = self.ui.targetNamesList.currentIndex()
#        if (dbglvl>1):
#            print currTargInd
        self.ui.displTargRA.setText(str(self.targetListCoords[0,currTargInd]))
        self.ui.displTargDec.setText(str(self.targetListCoords[1,currTargInd]))

    def showTime(self):
        # Might need to replace this with something more accurate.
        # currently uses machine time?? 
        self.sast = time.strftime("%H:%M:%S", time.localtime())
        self.gmt = time.strftime("%H:%M:%S", time.gmtime())
        #'Thu, 28 Jun 2001 14:17:15 +0000'
        self.ui.displSAST.setText(self.sast)
        self.ui.displUT.setText(self.gmt)
    """
    def passTargettoExposureName(self):
        # pass the target name (from the targets window) to CCD control.
        # just do it once at change to allow user to modify name (e.g. "<TARGET> direct image")

#        self.showTargets() # need to make sure reloaded changes!

        if(self.targetsLoaded == False):
            return
        currTargInd = self.ui.targetNamesList.currentIndex()
        self.settings['TARG-NAME'] = str(self.ui.targetNamesList.currentText())

        self.ui.displExpName.setText(self.settings['TARG-NAME'])

    def writeTargetTCS(self,targetLine):
        # write coords to target.dat
        # format: rah ram ras decd decm decs
        #def write_remote_file(server, user, password, remotepath, lines):
        server = 'tcs74v3.suth.saao.ac.za'
        user='ccd'
        passwd='Saaoccd'
        remotepath="/home/ccd/target.dat"
#        remotepath="/home/ccd/targetFromSpup.dat"
        try:
            write_remote_file(server, user, passwd, remotepath,targetLine)
            print 'written target.dat on tcs'
        except Exception as e:
            print("Caught exception: {}".format(e))


        pass

    def showTargets(self):
        # update target coords in target list
        if(self.targetsLoaded == False):
            return
        currTargInd = self.ui.targetNamesList.currentIndex()
#        if (dbglvl>1):
#            print currTargInd
        self.ui.displTargRA.setText(str(self.targetListCoords[0,currTargInd]))
        self.ui.displTargDec.setText(str(self.targetListCoords[1,currTargInd]))

        self.settings['TARG-RA'] = str(self.targetListCoords[0,currTargInd])
        self.settings['TARG-DEC'] = str(self.targetListCoords[1,currTargInd])
        self.settings['TARG-NAME'] = str(self.ui.targetNamesList.currentText())

        self.passTargettoExposureName()
        #--
        targetCooNoColons = (str(self.targetListCoords[0,currTargInd])+' '+str(self.targetListCoords[1,currTargInd])).replace(':',' ')
        print targetCooNoColons
        # round to nearest second so TCS doesn't throw a fit!:
        sline = targetCooNoColons.split(' ')
        targetTCSformat = '%s %s %.0f %s %s %.0f'%(sline[0],sline[1],float(sline[2]),sline[3],sline[4],float(sline[5]))
#        self.writeTargetTCS(targetCooNoColons+'\n')
        self.writeTargetTCS(targetTCSformat+'\n')
        #--


    def getTime(self):
        # Might need to replace this with something more accurate.
        self.sast = time.strftime("%H:%M:%S", time.localtime())
        self.gmt = time.strftime("%H:%M:%S", time.gmtime())
        from astropy.time import Time
        self.hjd = Time(Time.now()).mjd + 0.5
        self.ut = self.gmt
        self.timeSec = time.time()
        self.utdate = time.strftime("%Y-%m-%d", time.gmtime())

    def showTime(self):
        # currently uses machine time??
        self.getTime()
        #'Thu, 28 Jun 2001 14:17:15 +0000'
        self.ui.displSAST.setText(self.sast)
        self.ui.displUT.setText(self.gmt)
        self.ui.DisplDate.setText(self.utdate)

    def get_telescope_pos(self):
        server = 'tcs74v2.suth.saao.ac.za'
        user='ccd'
        passwd='Saaoccd'
        remotepath="/home/ccd/live.dat"
        remotepath="/home/ccd/archive/070715.live"
        localpath="live.dat"
        #from astropy.coordinates import Angle

        # need to fix angle problem and stop output from being so verbose!
        try:
            lines = read_remote_file(server, user, passwd, remotepath)
            lastline = get_last_line(lines)
            tcs = get_tcs_info(lastline)
            show_tcs_info(tcs)

        except Exception as e:
            print("Caught exception: {}".format(e))

