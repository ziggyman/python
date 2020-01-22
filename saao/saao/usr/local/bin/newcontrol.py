#!/usr/bin/python
# -*- coding: utf-8 -*-

import time
from PyQt4 import QtCore, QtGui

import logging, StringIO, time
import numpy as np

from cassspectr.controller.controller import Controller
from cassspectr.interface.spectrograph_interface.ttypes import DetectorException,PLCThriftException

import argparse


# let's set the most used options via command line arguments:
parser = argparse.ArgumentParser(description="")

#plc_host = 'localhost'
#plc_host = '172.17.0.133'
#plc_host = 'spupnic'
plc_port = 9090
#detector_host = 'localhost'
#detector_host = '172.17.0.133'
#detector_host = 'spupnic'
detector_port = 9091



parser.add_argument('--plc_server',
                        help="Specify the PLC server.", default='local',
                        choices=['local','remote'])
parser.add_argument('--det_server',
                        help="Specify the detector server.", default='local',
                        choices=['local','remote'])
parser.add_argument('--dbglvl',
                        help="Specify the debug level", default='2',
                        choices=['0','1','2','3','4','5'])
args = parser.parse_args()


if (args.plc_server=='local'):
    plc_host = 'localhost'
else:
#    plc_host = '172.17.0.133'
    plc_host = 'spupnic'

if (args.det_server=='local'):
    detector_host = 'localhost'
else:
    #detector_host = '172.17.0.133'
    detector_host = 'spupnic'

dbglvl = int(args.dbglvl)

print args
print detector_host
#stop()

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


#from newGUI import Ui_MainWindow
from cassspectr.ui.newGUI import Ui_MainWindow

from cassspectr.ui.instrPane import InstrPane
from cassspectr.ui.detPane import DetPane
from cassspectr.ui.targPane import TargPane
from cassspectr.ui.headerPane import HeaderPane
from cassspectr.ui.LogPane import logPane,logBuffer
from cassspectr.ui.engPane import EngPane
from cassspectr.ui.configPane import ConfigPane

from cassspectr.ui.mkClickable import clickable

# allow smart completion. following:
# http://stackoverflow.com/questions/4827207/how-do-i-filter-the-pyqt-qcombobox-items-based-on-the-text-input
#from PyQt4.QtGui import QComboBox, QApplication, QCompleter, QSortFilterProxyModel, QStandardItemModel, QStandardItem
#from PyQt4.QtGui import QCompleter,QSortFilterProxyModel
from PyQt4.QtCore import Qt
from PyQt4.QtGui import QCompleter, QComboBox, QSortFilterProxyModel


import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
#from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
##from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


#9AB87C
# Override default Ubuntu Unity ugly& confusing progress bar style:
# Linux Mint style
PROGRESSBAR_STYLE = """
QProgressBar{
    border: 1px solid black;
    border-radius: 2px;
    text-align: center
}

QProgressBar::chunk {
    background-color: #9AB87C;
    width: 10px;
    margin: 0px;
}
"""


#configdir = '/home/gilbank/Proj/CassSpect/sw/cassspectr/cassspectr/'
# assume we always run from cassspectr/ i.e. ui/../
###configdir = './config/'
configdir = '/home/ccd/.spup/'
###pngdir = '%s../ui/qtcreated/'%configdir
#pngdir = './ui/qtcreated/'
pngdir = '/etc/cassspectr/'

#rootdatadir = '/home/ccd/cassspectr/cassspectr/'
rootdatadir = '/home/ccd/data/'

#self.con.dc.initialize(512, 2148)
ccd_ny = 512
ccd_nx = 2148


# ----
# http://matplotlib.1069221.n5.nabble.com/key-press-events-in-matplotlib-embedded-in-pyqt4-td27958.html
class MplCanvas(FigureCanvas):
 
    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
 
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
 
 
class matplotlibWidget(QtGui.QWidget):
 
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        #--
#        self.mpl_toolbar = NavigationToolbar(self.canvas, None) 
#        self.vbl.addWidget(self.mpl_toolbar)
        #--
        self.setLayout(self.vbl)
# ----



"""
class MyPopup(QtGui.QWidget):
    def __init__(self,Picture1):
        QtGui.QWidget.__init__(self)

#        Picture1 = QtGui.QPixmap(  "/home/gilbank/tmp/uber3.jpg" )
        self.label = QtGui.QLabel( "Quick Look Viewer", self )
        self.label.setPixmap( Picture1 )
"""    



class ExtendedComboBox(QComboBox):
    def __init__(self, parent=None):
        super(ExtendedComboBox, self).__init__(parent)

        self.setFocusPolicy(Qt.StrongFocus)
        self.setEditable(True)

        # add a filter model to filter matching items
        self.pFilterModel = QSortFilterProxyModel(self)
        self.pFilterModel.setFilterCaseSensitivity(Qt.CaseInsensitive)
        self.pFilterModel.setSourceModel(self.model())

        # add a completer, which uses the filter model
        self.completer = QCompleter(self.pFilterModel, self)
        # always show all (filtered) completions
        self.completer.setCompletionMode(QCompleter.UnfilteredPopupCompletion)
        self.setCompleter(self.completer)

        # connect signals
        self.lineEdit().textEdited[unicode].connect(self.pFilterModel.setFilterFixedString)
        self.completer.activated.connect(self.on_completer_activated)


    # on selection of an item from the completer, select the corresponding item from combobox
    def on_completer_activated(self, text):
        if text:
            index = self.findText(text)
            self.setCurrentIndex(index)
            self.activated[str].emit(self.itemText(index))

    # on model change, update the models of the filter and completer as well
    def setModel(self, model):
        super(ExtendedComboBox, self).setModel(model)
        self.pFilterModel.setSourceModel(model)
        self.completer.setModel(self.pFilterModel)


    # on model column change, update the model column of the filter and completer as well
    def setModelColumn(self, column):
        self.completer.setCompletionColumn(column)
        self.pFilterModel.setFilterKeyColumn(column)
        super(ExtendedComboBox, self).setModelColumn(column)




class CSGui(QtGui.QMainWindow, InstrPane, DetPane, TargPane, HeaderPane, logPane, EngPane, ConfigPane):
    def __init__(self, parent=None, infile=None):
        QtGui.QWidget.__init__(self, parent)
        
        #set up the main UI
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.dbglvl = dbglvl

        self.targetsLoaded = False
        self.targetListCoords = None

        #--
        #self.ui.targetNamesList = QtGui.QComboBox(self.ui.gridLayoutWidget_6)
#        self.ui.targetNamesList = ExtendedCombo(self.ui.gridLayoutWidget_6)
        self.ui.targetNamesList = ExtendedComboBox(self.ui.targetNamesList)
        #--

        self.configLocked=True
#        self.gratingBrakeOn=False
        self.populateSlitWidth()

        self.detInterlock = 1 # detector software interlock. Don't allow any insturment movement unless detector is IDLE
                                # 1=OK -- Piet's convention

        self.configdir = configdir
        self.pngdir = pngdir
        self.rootdatadir = rootdatadir
###        self.readConfig() # populate results from config panel
        self.readConfigFile() # populate results from config panel
        self.applyConfig()

        self.autoDataDir()
###         self.manualDataDir()
        self.readLastDataDir()

        self.setBinning() # update display with binning, etc. labels first time
        self.ccdFullFrame = False # use windowed half frame for SCIENCE by default
###        self.getCCDPars()
###        print self.ccdpars['nRows']

#        self.logView = QtGui.QTextEdit(self.ui.labScrollLog)
#        self.logView.setReadOnly(True)



#        self.ui.image = np.random.random((255,2048))

        # ---- Quick Look Viewer:
#        image = np.random.random((255,2048))        
        self.ui.qlv = matplotlibWidget()
        self.ui.qlv.setGeometry(QtCore.QRect(100, 100, 1000, 500))
#        self.w.show()
#        self.w.canvas.ax.clear()
#        self.w.canvas.ax.imshow(image, interpolation='nearest', cmap='gist_gray')
#        self.w.show()
        
        gv = self.ui.graphicsView        
        gv.scene = QtGui.QGraphicsScene(self)
        self.ui.scene = gv.scene
        self.ui.gv = gv

        # Force update at startup
        self.frameNumPrefixChange()

        self.readDetectorINI()
#        self.ui.frameLog.setVisible(False)


        # Create the controller class
        #log.debug("Instantiating Controller.")    
#        try:
# Should probably just let the GUI fail to start if there's nothing to connect to!
        self.con = Controller(plc_host, plc_port, detector_host, detector_port)        
#        except Exception as e:
#            print("Caught PLC controller exception: {}".format(e))
            #return

#        self.con.dc.initialize(512, 512)
#        self.con.dc.initialize(512, 2098)
#        try:
#            ###self.con.dc.initialize(512, 2148)
        self.con.dc.initialize(ccd_ny, ccd_nx)
#        except Exception as e:
#            print("Caught PLC controller exception: {}".format(e))

        # create a (constant) timer object which will poll for events: 
        self.polltimer = QtCore.QTimer()

        self.lastTempPollTime = 0.0 - 600.0 # force check at startup

        # taken from: http://www.rkblog.rk.edu.pl/w/p/qtimer-making-timers-pyqt4/
        # start constant timer
#        self.pollIncrSec = 0.1
##        self.pollIncrSec = 1.0 
        self.pollIncrSec = 0.5
#        self.pollIncrSec = 10.0

        self.polltime=0.0
#        self.polltimer.start(100) # 100ms 
        self.polltimer.start(1000.0 * self.pollIncrSec) # 100ms 
        

        self.abortPollTime = -999. # keep track of when last detector abort command was issued!

        # copying from example in remote_plc_cli.py
#        log.debug("Instantiating PLC controller.")
        
        # PLCController is already instantiated as self.con.plcc


        # ---- Set up logging
        # see: https://docs.python.org/2/howto/logging-cookbook.html
        
#        logfile = 'CassSpect%s.log'%(time.strftime("%Y%m%d", time.gmtime()))
        ###logfile = '/var/log/cassspectr-ui%s.log'%(time.strftime("%Y%m%d", time.gmtime()))
        logfile = '%s/cassspectr-ui%s.log'%(self.datadir,time.strftime("%Y%m%d", time.gmtime()))
        logging.basicConfig(filename=logfile,level=logging.DEBUG,format='%(asctime)s %(message)s')
        logFormatter = logging.Formatter('%(asctime)s   %(levelname)-8s %(message)s', datefmt='%H:%M:%S')
        #logFormatter = logging.Formatter('%(name)-12s %(asctime)s   %(levelname)-8s %(message)s', datefmt='%H:%M:%S')
        self.logBuffer = logBuffer()
        self.logBuffer.bufferMessage.connect(self.on_logBuffer_bufferMessage)

        logHandler = logging.StreamHandler(self.logBuffer)
        logHandler.setFormatter(logFormatter)

        # add some colours:
#        logging.addLevelName( logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
#        logging.addLevelName( logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))
#        logging.addLevelName( logging.INFO, "\x1b[31m%s" % logging.getLevelName(logging.INFO)) #**


        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(logHandler)
        
        self.logtext = '' 
        #self.ui.scrollAreaWidgetContents = QtGui.QTextEdit()

        self.logger.info('GUI started')
        self.logger.info('Data dir: %s'%self.datadir)
        # ----


#        self.ui.exposing=False ; self.ui.finished=False ; self.ui.exped=0.0


        # ---- Set some palette colours
        self.ui.pNorm = QtGui.QPalette()
        self.ui.pNorm.setColor(QtGui.QPalette.Foreground,QtCore.Qt.black)
        self.ui.pCal = QtGui.QPalette()
        self.ui.pCal.setColor(QtGui.QPalette.Foreground,QtCore.Qt.blue)
        myOrange = QtGui.QColor(242,155,13)
#        self.ui.pCal.setColor(QtGui.QPalette.Foreground,myOrange)
#        self.ui.pCal.setColor(QtGui.QPalette.Background,QtCore.Qt.black)
        self.ui.pWarn = QtGui.QPalette()
        self.ui.pWarn.setColor(QtGui.QPalette.Foreground,QtCore.Qt.red)
        self.ui.pWait = QtGui.QPalette()
        self.ui.pWait.setColor(QtGui.QPalette.Foreground,myOrange)

        self.ui.lInvalid = QtGui.QPalette()
        self.ui.lInvalid.setColor(QtGui.QPalette.Text, QtCore.Qt.red)
        self.ui.lValid = QtGui.QPalette()
        self.ui.lValid.setColor(QtGui.QPalette.Text, QtCore.Qt.black)

        self.ui.progressBarDiskUsage.setStyleSheet(PROGRESSBAR_STYLE)
        self.ui.progressBarExpProgress.setStyleSheet(PROGRESSBAR_STYLE)


        # replace all below with status grabbed from plcc:

        QtCore.QObject.connect(self.ui.pushButArcMirrorToggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.ArcMirrorToggle)
        QtCore.QObject.connect(self.ui.pushButGMToggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.GMToggle)
        QtCore.QObject.connect(self.ui.pushButSlitIllumToggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.slitIllumToggle)

        QtCore.QObject.connect(self.ui.lineEditPropID, QtCore.SIGNAL(_fromUtf8("textChanged(QString)")), self.frameNumPrefixChange)

        QtCore.QObject.connect(self.ui.pushButExpose, QtCore.SIGNAL(_fromUtf8("pressed()")), self.doFirstExposure)
        QtCore.QObject.connect(self.ui.pushButAbort, QtCore.SIGNAL(_fromUtf8("pressed()")), self.doAbortExposure)
        QtCore.QObject.connect(self.ui.pushButStop, QtCore.SIGNAL(_fromUtf8("pressed()")), self.stopExposure)
        QtCore.QObject.connect(self.ui.pushButCancelSequence, QtCore.SIGNAL(_fromUtf8("pressed()")), self.cancelSequence)
        QtCore.QObject.connect(self.ui.comboExpType, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(QString)")), self.setExpTime)
        QtCore.QObject.connect(self.ui.displExpTime, QtCore.SIGNAL(_fromUtf8("textChanged(QString)")), self.expTimeChanged)
        QtCore.QObject.connect(self.ui.displNexp, QtCore.SIGNAL(_fromUtf8("valueChanged(QString)")), self.nExpChanged)
        QtCore.QObject.connect(self.ui.comboCCDmodeGain, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(QString)")), self.setCCDmode)
        QtCore.QObject.connect(self.ui.comboCCDbinning, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(QString)")), self.setBinning)
        QtCore.QObject.connect(self.ui.comboCCDmodeNoise, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(QString)")), self.setCCDmode)
        QtCore.QObject.connect(self.ui.comboCCDwindowing, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(QString)")), self.setCCDwindowing)
        QtCore.QObject.connect(self.ui.pushButFlushCCDPane, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(QString)")), self.flushCCDPane)

        QtCore.QObject.connect(self.ui.pushButSeqStart, QtCore.SIGNAL(_fromUtf8("pressed()")), self.scienceSequence)

        QtCore.QObject.connect(self.ui.pushButHartmannRun, QtCore.SIGNAL(_fromUtf8("pressed()")), self.HartmannSequence)

###        QtCore.QObject.connect(self.ui.pushButSlitShutterToggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.slitShutterToggle)
        QtCore.QObject.connect(self.ui.pushButHM1Toggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.HM1Toggle)
        QtCore.QObject.connect(self.ui.pushButSlitWidthGo, QtCore.SIGNAL(_fromUtf8("pressed()")), self.changeSlitWidth)
        QtCore.QObject.connect(self.ui.pushButFilterGo, QtCore.SIGNAL(_fromUtf8("pressed()")), self.filterGo)
        QtCore.QObject.connect(self.ui.pushButHM2Toggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.HM2Toggle)
        QtCore.QObject.connect(self.ui.pushButHMOpen, QtCore.SIGNAL(_fromUtf8("pressed()")), self.HMOpen)
        QtCore.QObject.connect(self.ui.pushButArcLamp1Toggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.ArcLamp1Toggle)
        QtCore.QObject.connect(self.ui.pushButArcLamp2Toggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.ArcLamp2Toggle)
        #QtCore.QObject.connect(self.ui.pushButSnapshot, QtCore.SIGNAL(_fromUtf8("pressed()")), self.createHeader)
        QtCore.QObject.connect(self.ui.pushButRoSToggle, QtCore.SIGNAL(_fromUtf8("pressed()")), self.RoSToggle)
        QtCore.QObject.connect(self.ui.pushButFocusGo, QtCore.SIGNAL(_fromUtf8("pressed()")), self.cameraFocusGo)
#        QtCore.QObject.connect(self.ui.pushButCameraFocusIncrementGo, QtCore.SIGNAL(_fromUtf8("pressed()")), self.cameraFocusIncrementGo)
        QtCore.QObject.connect(self.ui.pushButGratingAngleGo, QtCore.SIGNAL(_fromUtf8("pressed()")), self.gratingAngleGo)
#        QtCore.QObject.connect(self.ui.pushButGratingAngleIncrementGo, QtCore.SIGNAL(_fromUtf8("pressed()")), self.gratingAngleIncrementGo)
        QtCore.QObject.connect(self.ui.sliderSlitIllum, QtCore.SIGNAL(_fromUtf8("sliderReleased()")), self.slitIllumLevel)

#        QtCore.QObject.connect(self.ui.pushButGraphArcLamp1Toggle, QtCore.SIGNAL(_fromUtf8("clicked()")), self.ArcLamp1Toggle)

        QtCore.QObject.connect(self.ui.displExpTime, QtCore.SIGNAL(_fromUtf8("textChanged(QString)")), self.validateExpTime)
        #QtCore.QObject.emit(self.ui.displExpTime, QtCore.SIGNAL(_fromUtf8("textChanged()"))) # trigger once to set valid number
        #QtCore.QObject.emit(self.ui.displNexp, QtCore.SIGNAL(_fromUtf8("textChanged"))) #
        QtCore.QObject.connect(self.ui.displNexp, QtCore.SIGNAL(_fromUtf8("textChanged(QString)")), self.validateNexp)
        
        QtCore.QObject.connect(self.ui.lineEditGratingAngleDeg, QtCore.SIGNAL(_fromUtf8("textChanged(QString)")), self.validateGratingAngle)
        QtCore.QObject.connect(self.ui.lineEditFocus, QtCore.SIGNAL(_fromUtf8("textChanged(QString)")), self.validateFocus)
        

        QtCore.QObject.connect(self.ui.pushButFWInit, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doFilterWheelInit)    
        QtCore.QObject.connect(self.ui.pushButFWReset, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doFilterWheelReset)
        QtCore.QObject.connect(self.ui.pushButSlitWidthInit, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doSlitWidthInit)    
        QtCore.QObject.connect(self.ui.pushButSlitWidthReset, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doSlitWidthReset)
        # Camera focus init/reset to be implemented
        QtCore.QObject.connect(self.ui.pushButCameraFocusInit, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doCameraFocusInit)    
        QtCore.QObject.connect(self.ui.pushButCameraFocusReset, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doCameraFocusReset)
        QtCore.QObject.connect(self.ui.pushButGratingAngleInit, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doGratingAngleInit)    
        QtCore.QObject.connect(self.ui.pushButGratingAngleReset, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doGratingAngleReset)
        QtCore.QObject.connect(self.ui.pushButResetAll, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doAllReset)
        QtCore.QObject.connect(self.ui.pushButInitAll, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doAllInit)

        QtCore.QObject.connect(self.ui.pushButResetPCIController, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doResetPCIController)
        QtCore.QObject.connect(self.ui.pushButResetDetectorController, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doResetDetectorController)
        QtCore.QObject.connect(self.ui.pushButResetDetectorServer, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doResetDetectorServer)
        QtCore.QObject.connect(self.ui.pushButResetPLCServer, QtCore.SIGNAL(_fromUtf8("clicked()")), self.doResetPLCServer)

        QtCore.QObject.connect(self.ui.pushButPollTemp, QtCore.SIGNAL(_fromUtf8("clicked()")), self.pollTemperature)


#        QtCore.QObject.connect(self.ui.brakeButton, QtCore.SIGNAL(_fromUtf8("clicked()")), self.setBrake)
        QtCore.QObject.connect(self.ui.pushButUnlockConfig, QtCore.SIGNAL(_fromUtf8("clicked()")), self.unlockBox)
        QtCore.QObject.connect(self.polltimer, QtCore.SIGNAL("timeout()"), self.constantUpdate)
        QtCore.QObject.connect(self.ui.actionLoad_Target_Catalogue, QtCore.SIGNAL(_fromUtf8("activated()")), self.loadTargets)    
        QtCore.QObject.connect(self.ui.pushButLoadTargetCat, QtCore.SIGNAL(_fromUtf8("clicked()")), self.loadTargets)    
        QtCore.QObject.connect(self.ui.actionFilter_wheel_init, QtCore.SIGNAL(_fromUtf8("activated()")), self.doFilterWheelInit)    
        QtCore.QObject.connect(self.ui.actionFilter_wheel_reset, QtCore.SIGNAL(_fromUtf8("activated()")), self.doFilterWheelReset)
        QtCore.QObject.connect(self.ui.actionSlit_width_init, QtCore.SIGNAL(_fromUtf8("activated()")), self.doSlitWidthInit)    
        QtCore.QObject.connect(self.ui.actionSlit_width_reset, QtCore.SIGNAL(_fromUtf8("activated()")), self.doSlitWidthReset)
        QtCore.QObject.connect(self.ui.actionGrating_angle_init, QtCore.SIGNAL(_fromUtf8("activated()")), self.doGratingAngleInit)    
        QtCore.QObject.connect(self.ui.actionGrating_angle_reset, QtCore.SIGNAL(_fromUtf8("activated()")), self.doGratingAngleReset)
        QtCore.QObject.connect(self.ui.actionRESET_ALL, QtCore.SIGNAL(_fromUtf8("activated()")), self.doAllReset)
        QtCore.QObject.connect(self.ui.actionINIT_ALL, QtCore.SIGNAL(_fromUtf8("activated()")), self.doAllInit)
        QtCore.QObject.connect(self.ui.actionQuit, QtCore.SIGNAL(_fromUtf8("activated()")), self.closeEvent)


#        QtCore.QObject.connect(self.ui.targetNamesList, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(QString)")), self.passTargettoExposureName)
        QtCore.QObject.connect(self.ui.targetNamesList, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(QString)")), self.showTargets)


        QtCore.QObject.connect(self.ui.actionShow_hide_Log_Frame, QtCore.SIGNAL(_fromUtf8("activated()")), self.showHideLogFrame)
        QtCore.QObject.connect(self.ui.actionShow_hide_Target_Frame, QtCore.SIGNAL(_fromUtf8("activated()")), self.showHideTargetFrame)
        QtCore.QObject.connect(self.ui.actionShow_hide_Instrument_Setup_Frame, QtCore.SIGNAL(_fromUtf8("activated()")), self.showHideInstInfo)

        QtCore.QObject.connect(self.ui.pushButMakeDataDir, QtCore.SIGNAL(_fromUtf8("clicked()")), self.manualDataDir)
#        QtCore.QObject.connect(self.ui.pushButSetStartFrameNumber, QtCore.SIGNAL(_fromUtf8("clicked()")), self.setStartFrameNumber)
        QtCore.QObject.connect(self.ui.pushButConfirmNewRun, QtCore.SIGNAL(_fromUtf8("clicked()")), self.confirmNewRun)


        # schematic view graphics - make labels clickable
        clickable(self.ui.pushLabelGraphArcLamp1Toggle).connect(self.ArcLamp1Toggle)
        clickable(self.ui.pushLabelGraphArcLamp2Toggle).connect(self.ArcLamp2Toggle)
        clickable(self.ui.pushLabelGraphArcMirrorToggle).connect(self.ArcMirrorToggle)
        clickable(self.ui.pushLabelGraphRoSMirrorToggle).connect(self.RoSToggle)
        clickable(self.ui.pushLabelGraphHM1Toggle).connect(self.HM1Toggle)
        clickable(self.ui.pushLabelGraphHM2Toggle).connect(self.HM2Toggle)
        clickable(self.ui.pushLabelGraphGuideMirrorToggle).connect(self.GMToggle)
###        clickable(self.ui.pushLabelGraphSlitShutterToggle).connect(self.slitShutterToggle)
        clickable(self.ui.pushLabelGraphSlitIllumToggle).connect(self.slitIllumToggle)

        #--
        self.readConfigFile() # only read config FILE once at start [danger of overwriting if not everything is saved!]
        self.setCCDmode() # load gain/rdnoise vals at startup
        print self.settings
        self.lastErrors = [] # keep a list of error states generated
        #--

        # Disable grating angle changing by default. Enable in the config tab
#        self.ui.lineEditGratingAngleDeg.setEnabled(False)
#        self.ui.pushButGratingAngleGo.setEnabled(False)
#        self.ui.pushButGratingAngleInit.setEnabled(False)
	# Disable reset and init buttons by default. Enable in the config tab
        self.ui.pushButInitAll.setEnabled(True)
        self.ui.pushButResetPCIController.setEnabled(True)
        self.ui.pushButResetDetectorController.setEnabled(True)
        self.ui.pushButSlitWidthReset.setEnabled(True)
        self.ui.pushButGratingAngleReset.setEnabled(True)
        self.ui.pushButFWReset.setEnabled(True)
        self.ui.pushButCameraFocusReset.setEnabled(True)
        self.ui.pushButResetAll.setEnabled(True)
        #**** grey out some none-implemented things:
#        self.ui.pushButGratingAngleIncrementGo.setEnabled(False)
#        self.ui.lineEditGratingAngleIncrementDeg.setEnabled(False)
#        self.ui.pushButStop.setEnabled(False)

        # just turn off on start up. Once first exposure is taken, drawDet... will take care of this
        self.ui.pushButCancelSequence.setEnabled(False)

    def constantUpdate(self):
        """
        slot for constant timer timeout
        """
###        self.polltime = self.polltime+0.1 # keep track of time in seconds
        self.polltime = self.polltime + self.pollIncrSec # keep track of time in seconds
        if(dbglvl>2):
            self.logger.debug(self.polltime)
#            self.logger.info(self.polltime)
#        print '.'
        
        # update important stuff every visit:
        self.drawDetIndicators()
        self.showTime()
        if self.polltime % 1 == 0:  # only do every 1 second
            self.drawInstrIndicators()
            self.showPLCStatus()
            self.showGraphicView()
##            self.showTargets()
            self.showCCDstate()


        self.checkErrors()
###        self.get_telescope_pos()


        ## save GUI settings constantly
        if self.polltime % 1 ==0: # only do every 1 second
            self.readConfig()
            self.writeConfigFile()
        ##

        # update other stuff less frequently
        # every 2 mins:
#        pollIncrMilliSec = self.pollIncrSec/1000.0
        #if ( (self.polltime % 120.0) < 0.12): # need to be slightly more than 100ms for initial poll
        #if ( (self.polltime % 120.0) < (self.pollIncrSec+0.02)): # need to be slightly more than 100ms for initial poll
        if ( (self.polltime % 120.0) < (self.pollIncrSec*1.2)): # need to be slightly more than 100ms for initial poll
            self.diskUsage()

#        print self.polltime

        # every 10 mins:
#        if ( (self.polltime % 600.0) <0.12):
#        if ( (self.polltime % 600.0) < (self.pollIncrSec+0.02)):

        dTempPollTime = self.polltime - self.lastTempPollTime
#        print dTempPollTime
#        if ( (self.polltime % 600.0) < (self.pollIncrSec*1.2)):
        #doTpoll = 0
        # *** not thoroughly tested yet:

        #--
        if (dTempPollTime > 600.0):
            if self.dbglvl>2:
                print 'attempting to poll temperature'
            self.pollTemperature()
        #--

    def readDetectorINI(self):
        # get tempearture ranges from Carel's detector.ini file.
        # using ini readers, like ConfigParser, does not seem to work.
        # Let's just do this the dumb way:
        inifile='/usr/local/etc/detector.ini'
        f=open(inifile,'r')
        self.mincfTemp=81
        self.maxcfTemp=89
        self.minCCDTemp=165
        self.maxCCDTemp=169
        return

        for line in f.readlines():
            sline = line.split()
            if sline[0]=='cold_finger_min':
                self.mincfTemp=float(sline[2])
            elif sline[0]=='cold_finger_max':
                self.maxcfTemp=float(sline[2])
            elif sline[0]=='ccd_min':
                self.minCCDTemp=float(sline[2])
            elif sline[0]=='ccd_max':
                self.maxCCDTemp=float(sline[2])
        f.close()


    def checkErrors(self):
        # run get_errors and see if any new error has appeared since previous run:
        try:
            checkErrors = self.con.plcc.get_errors()
        except PLCThriftException, e:
            print("Caught PLC controller exception: {}".format(e))
            return
        except Exception, e:
            print("Caught exception: {}".format(e))
            return
        newErrors = [x for x in checkErrors if x not in self.lastErrors]
        for newError in newErrors:
            self.logger.warn(newError)
        if self.dbglvl>2:
            print newErrors
        # update erors. This should deal with errors which have disappeared, and also not assign if nothing new happened
        if ( (len(newErrors)>0) | (len(checkErrors)<>len(self.lastErrors)) ):
        #if len(newErrors) >0:
            self.lastErrors = checkErrors 

    def pollTemperature(self):
        # moved to separate routine to try to make easier to debug:


        # check we are not reading out, nor about to within N seconds
        N=2.0

        if self.lastCCDStatus==2:
            isReadingOut = True
        else:
            isReadingOut = False

        # **** Now change so that temp. can't be polled even while exposing (not just reading out). DBC's suggestion 26/10/2015:
        if(self.detInterlock==0):
            # exposing! locked:
            isReadingOut = True # Actually exposing, but still lock

        if(self.dbglvl>3):
            print 'isReadingOut = %s'%isReadingOut

        if ( not isReadingOut):
            # CCD not reading out.
            # safe to continue:
            try:
                ccdTemp = self.con.dc.get_ccd_temperature()
            except Exception as e:
                print("Caught PLC controller exception: {}".format(e))
                return
            try:
                cfTemp = self.con.dc.get_cold_finger_temperature()
            except Exception as e:
                print("Caught PLC controller exception: {}".format(e))
                return
            self.ui.labCCDTemp.setText('%.1fK'%(ccdTemp))
            self.ui.labColdFingerTemp.setText('%.1fK'%(cfTemp))
            self.logger.info('CCD temp. = %.2f'%((ccdTemp)))
            self.logger.info('cold finger temp. = %.2f'%((cfTemp)))
            self.lastTempPollTime = np.copy(self.polltime)

            # check temperatures are okay:
            if ( (ccdTemp>self.minCCDTemp) & (ccdTemp<self.maxCCDTemp) ):
                self.ui.labCCDTemp.setPalette(self.ui.pNorm)
            else:
                self.logger.warn('CCD temperature is outside specified limits: %.1f--%.1f'%(self.minCCDTemp,self.maxCCDTemp))
                self.ui.labCCDTemp.setPalette(self.ui.pWarn)
            
            if ( (cfTemp>self.mincfTemp) & (cfTemp<self.maxcfTemp) ):
                self.ui.labColdFingerTemp.setPalette(self.ui.pNorm)

            else:
                self.logger.warn('cold finger temperature is outside specified limits: %.1f--%.1f'%(self.mincfTemp,self.maxcfTemp))
                self.ui.labColdFingerTemp.setPalette(self.ui.pWarn)


    def populateSlitWidth(self):
        widths=np.linspace(0.,4.19,29) # from spmanv6.1.pdf pg 28 (Piet says there's a 0 width) 
        swidths = ['%.0f:  %.2f"'%(i,x) for i,x in enumerate(widths)]
        self.ui.swidths = swidths # store for later translation
        self.ui.comboBoxSlitWidth.clear()
        self.ui.comboBoxSlitWidth.addItems(swidths)


    def drawInstrIndicators(self):
        # display status of Instrument panel indicators based on PLC status object:

        if(self.dbglvl>4): print '1: self.con.plcc.get_status()'
        try:
            status = self.con.plcc.get_status()
        except PLCThriftException as e:
            print("Caught PLC Thrift exception: {}".format(e))
            self.logger.exception(e)
            return
        except Exception as e:
            print("Caught exception: {}".format(e))
            return

        # ---- Guide Mirror ----
        # this is not an else as it happens in conjunction with in/out status
        if (status['GMFailure']):  # Failure
            self.ui.labGMstate.setText(str('FAIL'))
            self.ui.labGMstate.setPalette(self.ui.pWarn)
        else:
            # no failure, check for status:
            if (status['GMCentred']):  # SCIENCE
                self.ui.labGMstate.setText(str('SCIENCE'))
                self.ui.labGMstate.setPalette(self.ui.pNorm)
                
            elif (status['GMInbeam']):  # ACQ
                self.ui.labGMstate.setText(str('ACQUISITION'))
                self.ui.labGMstate.setPalette(self.ui.pCal)
    
    
            elif (status['GMMoving']):  # Moving
                
                # **** Deactivate button while moving so user cannot press again!
                self.ui.labGMstate.setText(str('MOVING'))
                self.ui.labGMstate.setPalette(self.ui.pWait)
            
          
        
        
        # ---- RoS Mirror ----
        if (status['RoSMirrorFailure']==1):
            RoSMirrorState='FAIL'
            self.ui.labRoSMirrorState.setPalette(self.ui.pWarn)
        else:

            if (status['RearOfSlitMirror']==0):  # out of beam
                RoSMirrorState='Out of beam'
                self.ui.labRoSMirrorState.setPalette(self.ui.pNorm)
            elif (status['RearOfSlitMirror']==1):  # in beam
                self.ui.labRoSMirrorState.setPalette(self.ui.pCal)
                RoSMirrorState='In beam'

        self.ui.labRoSMirrorState.setText(str(RoSMirrorState))

        
        # ----
        
        # ---- Slit width
        if (status['SlitWidthFailure']==1):
            self.ui.labSlitWidthCurrent.setPalette(self.ui.pWarn)
            self.ui.labSlitWidthCurrent.setText('FAIL')
        elif (status['SlitWidthMoving']==1):
            self.ui.labSlitWidthCurrent.setPalette(self.ui.pWait)
        elif (status['SlitWidthPosition']==0):
            self.ui.labSlitWidthCurrent.setPalette(self.ui.pCal)
        else:
            self.ui.labSlitWidthCurrent.setPalette(self.ui.pNorm)

        # ----


        # ---- slit illumination
        if (status['SlitIllumination']==1):  # calib
            self.ui.labSlitIllumState.setPalette(self.ui.pCal)
            slitIllumState='ON'
        elif (status['SlitIllumination']==0): # off, science
            self.ui.labSlitIllumState.setPalette(self.ui.pNorm)
            slitIllumState='OFF'
        self.ui.labSlitIllumState.setText(slitIllumState)
        
        # ----
        
        
        # ---- Arc Mirror ----
        if(status['ARCMirrorFailure']==1):
            arcMirrorState='FAIL'
            self.ui.labArcMirrorStatus.setPalette(self.ui.pWarn)
        else:
            if (status['ARCMirror']==0):  # out of beam
                arcMirrorState='Out of beam'
                self.ui.labArcMirrorStatus.setPalette(self.ui.pNorm)
            elif (status['ARCMirror']==1):  # in beam
                self.ui.labArcMirrorStatus.setPalette(self.ui.pCal)
                arcMirrorState='In beam'
            
        self.ui.labArcMirrorStatus.setText(str(arcMirrorState))


        # -- Arc Lamps --
        if (status['ARC1']==False):
            self.ui.labArcLamp1State.setText(str('OFF'))
            self.ui.labArcLamp1State.setPalette(self.ui.pNorm)
        elif (status['ARC1']==True):
            self.ui.labArcLamp1State.setText(str('ON'))
            self.ui.labArcLamp1State.setPalette(self.ui.pCal)

        if (status['ARC2']==False):
            self.ui.labArcLamp2State.setText(str('OFF'))
            self.ui.labArcLamp2State.setPalette(self.ui.pNorm)
        elif (status['ARC2']==True):
            self.ui.labArcLamp2State.setText(str('ON'))
            self.ui.labArcLamp2State.setPalette(self.ui.pCal)


        # -- check if an arc is configured and update its name for the header:
        if (status['ARCMirror']==1):
            self.settings['ARC']=''
            # allow for combination in case both are on??
            if (status['ARC1']==1):
                self.settings['ARC']+=str(self.ui.lineEditArc1ID.text())
            if (status['ARC2']==1):
                self.settings['ARC']+=str(self.ui.lineEditArc2ID.text())
        else:
            self.settings['ARC']='OFF'



        # changed how this works. The state can only be A/B or Open really
        # -- Hartmann shutter
        if(status['HartmanFailure']==True):
                self.ui.labHartmann1State.setText(str('FAIL'))
                self.ui.labHartmann1State.setPalette(self.ui.pWarn)            
        else:
            if (status['HartmanA']==True):
                self.ui.labHartmann1State.setText(str('A'))
                self.ui.labHartmann1State.setPalette(self.ui.pCal)
            elif (status['HartmanB']==True):
                self.ui.labHartmann1State.setText(str('B'))
                self.ui.labHartmann1State.setPalette(self.ui.pCal)
            elif ( (status['HartmanA']==False) & (status['HartmanB']==False) ):
                self.ui.labHartmann1State.setText(str('Open'))
                self.ui.labHartmann1State.setPalette(self.ui.pNorm)

        self.settings['HartmannState'] = str(self.ui.labHartmann1State.text())

        # illuminate frame when hartmann statuses are hidden behind tab:
        if ((status['HartmanA']==True) | (status['HartmanB']==True)):
            self.ui.frameHartmann.setPalette(self.ui.pCal)
        elif (status['HartmanFailure']==True):
            self.ui.frameHartmann.setPalette(self.ui.pWarn)
        else:
            self.ui.frameHartmann.setPalette(self.ui.pNorm)


        # slit width:
        swidthAS = self.ui.swidths[ status['SlitWidthPosition'] ] # check indexing later
        self.ui.labSlitWidthCurrent.setText(str(swidthAS))
    
        # ---- grating 
        if (status['GratingInserted']==False):
            #grName = 'None'
            self.settings['grName']='None'
        else:
            # ID grating
            gratingID = status['GratingID']
            self.getGratingName(gratingID)

        self.ui.labGratingName.setText(self.settings['grName'])
        self.ui.labGratingAngleCurrent.setText('{:03.2f}'.format(status['GratingAngle']))
        # ----
        
        # ---- Filter wheel
        currFiltKeyword = 'FiltPos%s'%(status['FilterwheelPosition'])
#        if (status['FilterwheelPosition']==-1):
        if (status['FilterwheelPosition']<1): #***
            # do sthg***!
            currFiltKeyword = 'FiltPos1' #??
        self.ui.labFilterStateCurrent.setText(self.settings[currFiltKeyword])
        if(status['FilterMoving']==1):
            self.ui.labFilterStateCurrent.setText('MOVING')
            self.ui.labFilterStateCurrent.setPalette(self.ui.pWait)
        else:
            self.ui.labFilterStateCurrent.setPalette(self.ui.pNorm)
        # ----
        
        # ---- Focus Position
        self.settings['focusPosition'] = status["FocusPosition"]
        if(status['CameraFocusMoving']==1):
            self.ui.labFocusCurrent.setText('MOVING')
            self.ui.labFocusCurrent.setPalette(self.ui.pWait)
        else:
            self.ui.labFocusCurrent.setPalette(self.ui.pNorm)
            self.ui.labFocusCurrent.setText('{:03.3f}'.format(self.settings['focusPosition']))
        # ---
    
    

    
    
    def drawDetIndicators(self):
        #self.diskUsage()
        self.displDataDir()
        self.drawExp()

        ##** for now - DISABLED
        self.ui.pushButStop.setEnabled(False)
        ##**

        # -- exposure status:
        try:
            status = self.con.dc.get_exposure_status()
        except Exception as e:
            self.logger.exception(e)
            return

        if (status == 0):
            self.ui.labCCDState.setText('IDLE')
            self.ui.pushButAbort.setEnabled(False)
            self.ui.pushButStop.setEnabled(False)
        elif (status ==1):
            self.ui.labCCDState.setText('EXPOSING')
            # check time remaining from CCD and c.f. sw prediction
            # if they differ by more than 5 sec shout CCD FAILURE!
            remaining_time = self.con.dc.get_remaining_time()
            expTimeDiff = np.abs(self.sw_remaining_time - remaining_time)

            """
##            if( (expTimeDiff>5.0) & (np.abs(self.polltime-self.abortPollTime)>3.0) ):
            if( (expTimeDiff>5.0) & (np.abs(self.polltime-self.abortPollTime)>3.0) & (self.settings['expType']<>'HARTMANN') ):
                #*** UNLESS AN ABORT COMMAND WAS JUST ISSUED!
                #*** AND PROVIDED WE'RE NOT IN A HARTMANN SEQUENCE WHERE THINGS GO SCREWY WHEN THE SOFTWARE CLOCK GOES OUT OF SYNC FROM BEING TRAMPLED ALL OVER!

                print 'CCD failure!'
                #*** do sthg! ***
                
                self.logger.warning('CCD problem: expected detector to report %.1fs remaining\nit reported %.1fs remaining\n'\
                                    %(self.sw_remaining_time,remaining_time))

                #*** is this the only place a CCD failure is detected? What about propagating exceptions from Carel's code?***
                self.ui.frameCCD.setLineWidth(5)
                self.ui.frameCCD.setPalette(self.ui.pWarn)
            else:
                self.ui.frameCCD.setLineWidth(0)
                self.ui.frameCCD.setPalette(self.ui.pNorm)
            """

            if(self.dbglvl>2):
                print expTimeDiff
                
        elif (status==2):
            self.ui.labCCDState.setText('READING OUT')
    
        
    
    @QtCore.pyqtSlot(str)
    def on_logBuffer_bufferMessage(self, message):
#        self.logtext = self.logtext+'\n'+message
        self.logtext = message

#        self.logView = QtGui.QTextEdit(self.ui.labScrollLog)
#        self.logView.setReadOnly(True)

        # switching to the method I used for qlgui.py
        # solves the problem with poor auto-scrolling, etc.
        cursor = self.ui.logView.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(self.logtext)
        self.ui.logView.setTextCursor(cursor)
        self.ui.logView.ensureCursorVisible()


    def showHideLogFrame(self):
        if self.ui.frameLog.isVisible():
            self.ui.frameLog.setVisible(False)
        else:
            self.ui.frameLog.setVisible(True)

    def showHideTargetFrame(self):
        if self.ui.frameTarget.isVisible():
            self.ui.frameTarget.setVisible(False)
        else:
            self.ui.frameTarget.setVisible(True)

    def showHideInstInfo(self):
        if self.ui.frameInstInfo.isVisible():
            self.ui.frameInstInfo.setVisible(False)
        else:
            self.ui.frameInstInfo.setVisible(True)

    def loadTargets(self):
#        print "foo"
#        fname = QtGui.QFileDialog.getOpenFileName(None,'Open file', str(self.ui.displCatName))
        #fname = QtGui.QFileDialog.getOpenFileName(None,'Open file', self.ui.displCatName)
        fname = QtGui.QFileDialog.getOpenFileName(None,'Open file', '/home/ccd/')

#        self.ui.displCatName.setText(displCatName)

        shortname = (fname.split('/'))[-1]
        if(dbglvl>1):
            #print 'trying to read: %s'%fname
            self.logger.debug('trying to read: %s'%fname)
            self.logger.debug(shortname)
        f = open(fname, 'r')
 
        try:
            targname,targra,targdec = np.loadtxt(f,dtype='a',unpack=True,usecols=(0,1,2))
        except:
            print 'could not load %s'%shortname
            return
        
        # build dictionary of coords for each name:
        
#        targDictRA = dict(zip(targname,targra))
#        targDictDec = dict(zip(targname,targdec))        
        
        
        # -- Update target info:        
        self.ui.displCatName.setText(shortname)
        self.ui.labelLoadedFrom.setText('loaded from %s'%shortname)
        self.ui.targetNamesList.clear()
        self.ui.targetNamesList.insertItems(0,targname) 
        ##--
        # allow smart completion. following:
        # http://stackoverflow.com/questions/4827207/how-do-i-filter-the-pyqt-qcombobox-items-based-on-the-text-input
        """
        model = QStandardItemModel()
        for iii,word in enumerate( targname ):
            item = QStandardItem(word)
            model.setItem(iii, 0, item)
        self.ui.targetNamesList.setModel(model)
        self.ui.targetNamesList.setModelColumn(0)
        self.ui.targetNamesList.show()
        """
        combo = self.ui.targetNamesList
        combo.addItems(targname)
        # or use another model
        #combo.setModel(QStringListModel(string_list))
        combo.resize(300, 40)
        combo.show()

        ##--

        self.targetsLoaded=True
        self.targetListCoords = np.array((targra, targdec))

        self.showTargets()




    # Reset menu items
    def doFilterWheelInit(self):
        self.logger.info('Sending filter wheel INIT')
        self.con.plcc.set_fw_init()

    def doFilterWheelReset(self):
        self.logger.info('Sending filter wheel RESET')
        self.con.plcc.set_fw_reset()
    
    def doSlitWidthInit(self):
        self.logger.info('Sending slit INIT')
        self.con.plcc.set_slit_init()

    def doSlitWidthReset(self):
        self.logger.info('Sending slit RESET')
        self.con.plcc.set_slit_reset()
            
    def doCameraFocusInit(self):
        self.logger.info('Sending camera focus INIT')
        self.con.plcc.set_camera_focus_init()

    def doCameraFocusReset(self):
        self.logger.info('Sending camera focus RESET')
        self.con.plcc.set_camera_focus_reset()

    def doGratingAngleInit(self):
        self.logger.info('Sending grating angle INIT')
        self.con.plcc.set_grating_angle_init()

    def doGratingAngleReset(self):
        self.logger.info('Sending grating angle RESET')
        self.con.plcc.set_grating_angle_reset()


    def doAllInit(self):
        self.logger.info('Sending ALL INIT')
        self.doFilterWheelInit()
        self.doSlitWidthInit()
        self.doGratingAngleInit()

    def doAllReset(self):
        self.logger.info('Sending ALL RESET')
        self.doFilterWheelReset()
        self.doSlitWidthReset()
        self.doGratingAngleReset()

    def doResetPCIController(self):
        self.logger.info('Sending restore connection to Detector Server')
        self.con.dc.reset_device()

    def doResetDetectorController(self):
        self.logger.info('Sending restore connection to Detector Server')
        self.con.dc.reset_controller()



    def doResetDetectorServer(self):
        self.logger.info('Sending restore connection to Detector Server')
        print 'reset_connection'
        try:
            self.con.dc.reset_connection()
            #        self.con.dc.initialize(512, 2148)
            print 'initialise'
            self.con.dc.initialize(ccd_ny, ccd_nx)

            # change panels back to normal:
            self.ui.frameCCD.setLineWidth(0)
            self.ui.frameCCD.setPalette(self.ui.pNorm)
            # return exposure buttons to idle state
            self.doneAbort()
            self.logger.info('Detector server connection restored')
        except Exception as e:
            
            self.ui.frameCCD.setLineWidth(5)
            self.ui.frameCCD.setPalette(self.ui.pWarn)
            print '***'
            return

    def doResetPLCServer(self):
        self.logger.info('Sending restore connection to PLC Server')
        self.con.plcc.reset_connection()


    def getGratingName(self,gratingID):
        self.settings['grName']=''
#        if (str(gratingID)=='00000'):
        if (gratingID==0): # Empty!***
#            print '00000'
            # use backup name:
            grName = str(self.ui.lineEditManGratingName.text())
        else:        
            # translate ID to name*** # not tested yet
            grName = 'gr%i'%gratingID #'TODO'
        self.settings['grName'] = grName

    
    def validateExpTime(self):
        if(self.dbglvl>2):
            print 'changed!'
        validator = QtGui.QDoubleValidator(bottom=0)
#        validator.setDecimals(0)
        state = validator.validate(self.ui.displExpTime.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            self.ui.displExpTime.setPalette(self.ui.lValid)
            self.ui.lastValidExpTime = self.ui.displExpTime.text()
#            self.ui.displExpTime.setText(str(int(self.ui.displExpTime.text())))
        else:
            self.ui.displExpTime.setPalette(self.ui.lInvalid)
            try:
                self.ui.displExpTime.setText(self.ui.lastValidExpTime) # this might be a better way to deal actually - never allow invalid
            except:
                self.ui.displExpTime.setText("10") # if hasn't been set initially
#        self.ui.displExpTime.setText('111')
#        pass



    #*** need to also validate frame number ***

    def validateNexp(self):
        if(self.dbglvl>2):
            print 'changed!'
        validator = QtGui.QIntValidator(bottom=1)
        state = validator.validate(self.ui.displNexp.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            self.ui.displNexp.setPalette(self.ui.lValid)
            self.ui.lastValidNexp = self.ui.displNexp.text()
        else:
            self.ui.displNexp.setPalette(self.ui.lInvalid)
            try:
                self.ui.displNexp.setText(self.ui.lastValidNexp) # this might be a better way to deal actually - never allow invalid
            except:
                self.ui.displNexp.setText("1")
    
    def validateGratingAngle(self):
        if(self.dbglvl>2):
            print 'changed!'
        validator = QtGui.QDoubleValidator(bottom=-13.0,top=25.0)
        state = validator.validate(self.ui.lineEditGratingAngleDeg.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            self.ui.lineEditGratingAngleDeg.setPalette(self.ui.lValid)

        else:
            self.ui.lineEditGratingAngleDeg.setPalette(self.ui.lInvalid)

    def validateFocus(self):
        if(self.dbglvl>2):
            print 'changed!'
        validator = QtGui.QDoubleValidator(bottom=1.0,top=8.0, decimals=3)
        state = validator.validate(self.ui.lineEditFocus.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            self.ui.lineEditFocus.setPalette(self.ui.lValid)

        else:
            self.ui.lineEditFocus.setPalette(self.ui.lInvalid)

   
   


    def closeEvent(self):

        quit_msg = "Are you sure you want to quit?"
        #reply = QtGui.QMessageBox.question(self, 'Message', quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        reply = QtGui.QMessageBox.question(None, 'Removal', quit_msg, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if (reply == QtGui.QMessageBox.Yes):
#            QtGui.QApplication.quit() # Nope!
            QtGui.qApp.closeAllWindows()
        else:
            pass


def main():
    import sys
    app = QtGui.QApplication(['SpUpNIC Control'])
    myapp = CSGui()
    myapp.show()
    myapp.setWindowTitle('SpUpNIC Control')
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
