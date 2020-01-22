#/usr/bin/env python

# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'QLgui/mainwindow.ui'
#
# Created: Tue Jan 22 16:18:30 2013
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from QtVariant import QtCore, QtGui
from QxtSpanSlider import QxtSpanSlider

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from pylab import imshow,close
import numpy as np
from scipy.optimize import curve_fit

#import pyfits,time
import time

import os
from pyinotify import WatchManager, Notifier, ThreadedNotifier, EventsCodes, ProcessEvent

from zscale import zscale
#from peakdet import peakdet,gaussian
from specutils import *


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


idebugGUI=False 



from qlGUI import Ui_MainWindow


# http://stackoverflow.com/questions/8356336/how-to-capture-output-of-pythons-interpreter-and-show-in-a-text-widget
class EmittingStream(QtCore.QObject):

    textWritten = QtCore.pyqtSignal(str)

    def write(self, text):
        self.textWritten.emit(str(text))




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
 
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        #--
#        self.mpl_toolbar = NavigationToolbar(self.canvas, None) 
#        self.vbl.addWidget(self.mpl_toolbar)
        #--
        self.setLayout(self.vbl)
# ----
class matplotlibWidgetwBar(QtGui.QWidget):
 
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        #--
        self.mpl_toolbar = NavigationToolbar(self.canvas, None) 
        self.vbl.addWidget(self.mpl_toolbar)
        #--
        self.setLayout(self.vbl)

        
        
        
        
        
        
        
        
        
        
        self.extractWin_x0=self.extractWinSlider.lowerValue
        self.extractWin_x1=self.extractWinSlider.upperValue
        self.SB1_x0=self.SB1Slider.lowerValue
        self.SB1_x1=self.SB1Slider.upperValue
        self.SB2_x0=self.SB2Slider.lowerValue
        self.SB2_x1=self.SB2Slider.upperValue
 

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        self.autoExtraction=True
        self.lamp=None
        self.grating=None

#        self.raw2DGraph = QtGui.QGraphicsView(self.centralWidget)
#        self.raw2DGraph = QtGui.QWidget(self.centralWidget)
        self.raw2DGraph = matplotlibWidget(self.centralWidget)
        self.raw2DGraph.setGeometry(QtCore.QRect(20, 40, 1061, 81))
        self.raw2DGraph.setObjectName(_fromUtf8("raw2DGraph"))

        
        self.skysub2DView = matplotlibWidget(self.centralWidget)
        self.skysub2DView.setGeometry(QtCore.QRect(20, 160, 1061, 81))
        self.skysub2DView.setObjectName(_fromUtf8("skysub2DView"))

        self.useSB1Check = QtGui.QCheckBox(self.centralWidget)
        self.useSB1Check.setGeometry(QtCore.QRect(1110, 560, 85, 21))
        self.useSB1Check.setChecked(True)
        self.useSB1Check.setObjectName(_fromUtf8("useSB1Check"))
        self.useSB2Check = QtGui.QCheckBox(self.centralWidget)
        self.useSB2Check.setGeometry(QtCore.QRect(1260, 560, 85, 21))
        self.useSB2Check.setChecked(True)
        self.useSB2Check.setObjectName(_fromUtf8("useSB2Check"))
        
        
        # ----
        
        # ---- 1D extracted spec
        
        self.spec1DView = matplotlibWidgetwBar(self.centralWidget)
        self.spec1DView.setGeometry(QtCore.QRect(20, 290, 1061, 381))
        self.spec1DView.setObjectName(_fromUtf8("spec1DView"))

        # ----

        # ---- 2D spec related stuff
        self.z0Fill = QtGui.QLineEdit(self.centralWidget)
        self.z0Fill.setGeometry(QtCore.QRect(1140, 50, 113, 23))
        self.z0Fill.setObjectName(_fromUtf8("z0Fill"))
        self.z1Fill = QtGui.QLineEdit(self.centralWidget)
        self.z1Fill.setGeometry(QtCore.QRect(1290, 50, 113, 23))
        self.z1Fill.setObjectName(_fromUtf8("z1Fill"))
        self.ScalingCombo = QtGui.QComboBox(self.centralWidget)
        self.ScalingCombo.setGeometry(QtCore.QRect(1160, 80, 78, 24))
        self.ScalingCombo.setObjectName(_fromUtf8("ScalingCombo"))
        self.ScalingCombo.addItem(_fromUtf8(""))
        self.ScalingCombo.addItem(_fromUtf8(""))
        self.ScalingCombo.addItem(_fromUtf8(""))
        self.zScaleUpdateBut = QtGui.QPushButton(self.centralWidget)
        self.zScaleUpdateBut.setGeometry(QtCore.QRect(1280, 80, 99, 24))
        self.zScaleUpdateBut.setObjectName(_fromUtf8("zScaleUpdateBut"))

        # ----



        # ---- log window and file loaders

#        self.logView = QtGui.QListView(self.centralWidget)
        self.logView = QtGui.QTextEdit(self.centralWidget)
        self.logView.setReadOnly(True)
#        self.logView.setGeometry(QtCore.QRect(20, 710, 531, 192))
        self.logView.setGeometry(QtCore.QRect(20, 710, 731, 192))
        self.logView.setObjectName(_fromUtf8("logView"))
        self.CurrFrameBut = QtGui.QPushButton(self.centralWidget)
        #self.CurrFrameBut.setGeometry(QtCore.QRect(740, 720, 99, 24))
        self.CurrFrameBut.setGeometry(QtCore.QRect(940, 720, 99, 24))
        self.CurrFrameBut.setObjectName(_fromUtf8("CurrFrameBut"))
        self.currFrameFill = QtGui.QLineEdit(self.centralWidget)
        #self.currFrameFill.setGeometry(QtCore.QRect(600, 720, 131, 23))
        self.currFrameFill.setGeometry(QtCore.QRect(800, 720, 131, 23))
        self.currFrameFill.setObjectName(_fromUtf8("currFrameFill"))
        self.curArcBut = QtGui.QPushButton(self.centralWidget)
        self.curArcBut.setEnabled(False)
        self.curArcBut.setGeometry(QtCore.QRect(940, 770, 99, 24))
        self.curArcBut.setObjectName(_fromUtf8("curArcBut"))
        self.currArcFill = QtGui.QLineEdit(self.centralWidget)
        self.currArcFill.setGeometry(QtCore.QRect(800, 770, 131, 23))
        self.currArcFill.setObjectName(_fromUtf8("currArcFill"))
        self.curFlatBut = QtGui.QPushButton(self.centralWidget)
        self.curFlatBut.setEnabled(False)
        self.curFlatBut.setGeometry(QtCore.QRect(940, 820, 99, 24))
        self.curFlatBut.setObjectName(_fromUtf8("curFlatBut"))
        self.currFlatFill = QtGui.QLineEdit(self.centralWidget)
        self.currFlatFill.setGeometry(QtCore.QRect(800, 820, 131, 23))
        self.currFlatFill.setObjectName(_fromUtf8("currFlatFill"))
        self.curFluxCalBut = QtGui.QPushButton(self.centralWidget)
        self.curFluxCalBut.setEnabled(False)
        self.curFluxCalBut.setGeometry(QtCore.QRect(940, 870, 99, 24))
        self.curFluxCalBut.setObjectName(_fromUtf8("curFluxCalBut"))
        self.currFluxCalFill = QtGui.QLineEdit(self.centralWidget)
        self.currFluxCalFill.setGeometry(QtCore.QRect(800, 870, 131, 23))
        self.currFluxCalFill.setObjectName(_fromUtf8("currFluxCalFill"))
        
        
        # -- temp stuff fold old spectrograph only (no grating keyword in header)
        self.gratingCombo = QtGui.QComboBox(self.centralWidget)
        self.gratingCombo.setGeometry(QtCore.QRect(1150, 790, 171, 24))
        self.gratingCombo.setObjectName(_fromUtf8("gratingCombo"))
        self.gratingCombo.addItem(_fromUtf8(""))
        self.gratingCombo.addItem(_fromUtf8(""))
        self.gratingCombo.addItem(_fromUtf8(""))
        self.label_15 = QtGui.QLabel(self.centralWidget)
        self.label_15.setGeometry(QtCore.QRect(1130, 760, 141, 16))
        self.label_15.setObjectName(_fromUtf8("label_15"))
        # --
        
        
        
        
        self.parentDataDir=None #'/home/gilbank/Proj/CassSpect/data/Deatrick/'
        self.imagetype=""
        self.lastWLCFile=None
        self.wlc=False
        self.latestFrame=None ##self.parentDataDir+'a3000142.fits'
        #self.ReadOutDelay = 2.0 # seconds # not used now
        
        # ----
        
        
        self.widget = QtGui.QWidget(self.centralWidget)

        """
        self.widget.setGeometry(QtCore.QRect(740, 130, 333, 23))
        self.widget.setObjectName(_fromUtf8("widget"))
        """
        
        self.horizontalLayout = QtGui.QHBoxLayout(self.widget)
        self.horizontalLayout.setMargin(0)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.showSBCheck = QtGui.QCheckBox(self.widget)
        self.showSBCheck.setChecked(True)
        self.showSBCheck.setObjectName(_fromUtf8("showSBCheck"))
        self.horizontalLayout.addWidget(self.showSBCheck)
        self.showEWCheck = QtGui.QCheckBox(self.widget)
        self.showEWCheck.setChecked(True)
        self.showEWCheck.setObjectName(_fromUtf8("showEWCheck"))
        self.horizontalLayout.addWidget(self.showEWCheck)
        MainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QtGui.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 1437, 21))
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuFile = QtGui.QMenu(self.menuBar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtGui.QToolBar(MainWindow)
        self.mainToolBar.setObjectName(_fromUtf8("mainToolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtGui.QStatusBar(MainWindow)
        self.statusBar.setObjectName(_fromUtf8("statusBar"))
        MainWindow.setStatusBar(self.statusBar)
        self.actionSet_data_dir = QtGui.QAction(MainWindow)
        self.actionSet_data_dir.setObjectName(_fromUtf8("actionSet_data_dir"))
        self.actionQuit = QtGui.QAction(MainWindow)
        self.actionQuit.setObjectName(_fromUtf8("actionQuit"))
        self.actionLoadData = QtGui.QAction(MainWindow)
        self.actionLoadData.setObjectName(_fromUtf8("actionLoadData"))
        self.menuFile.addAction(self.actionSet_data_dir)
        self.menuFile.addAction(self.actionLoadData)
        self.menuFile.addAction(self.actionQuit)
        self.menuBar.addAction(self.menuFile.menuAction())

        self.scaledImage=None

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.showEWCheck, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.displayRaw2D)
        QtCore.QObject.connect(self.showSBCheck, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.displayRaw2D)
        QtCore.QObject.connect(self.zScaleUpdateBut, QtCore.SIGNAL(_fromUtf8("clicked()")), self.updatezScaling)
        QtCore.QObject.connect(self.autoExtractBut, QtCore.SIGNAL(_fromUtf8("clicked()")), self.updateAutoExtraction)
        #QtCore.QObject.connect(self.showSBCheck, QtCore.SIGNAL(_fromUtf8("clicked()")), MainWindow.show)
        QtCore.QObject.connect(self.SB1Slider, QtCore.SIGNAL(_fromUtf8("sliderReleased()")), self.readExtractionSliders)
        QtCore.QObject.connect(self.SB2Slider, QtCore.SIGNAL(_fromUtf8("sliderReleased()")), self.readExtractionSliders)
        QtCore.QObject.connect(self.extractWinSlider, QtCore.SIGNAL(_fromUtf8("sliderReleased()")), self.readExtractionSliders)
        QtCore.QObject.connect(self.useSB1Check, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.updateAll)
        QtCore.QObject.connect(self.useSB2Check, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.updateAll)
        
        QtCore.QObject.connect(self.actionLoadData, QtCore.SIGNAL(_fromUtf8("activated()")), self.showDialog)
#        QtCore.QObject.connect(self.CurrFrameBut, QtCore.SIGNAL(_fromUtf8("clicked()")), self.showDialog)
        
        QtCore.QObject.connect(self.actionSet_data_dir, QtCore.SIGNAL(_fromUtf8("activated()")), self.showDialogSetDatadir)        
        QtCore.QObject.connect(self.actionQuit, QtCore.SIGNAL(_fromUtf8("activated()")), self.closeEvent)        
        
        QtCore.QMetaObject.connectSlotsByName(MainWindow)        


        if(self.parentDataDir==None):
            self.showDialogSetDatadir()

        # -- pyinotify code for monitoring dir for new files:
        wm = WatchManager()
#        notifier = ThreadedNotifier(wm, PTmp())
        notifier = ThreadedNotifier(wm, self.process_IN_CREATE)
#        notifier = ThreadedNotifier(wm, self.processMonitorData)
#        notifier = Notifier(wm, PTmp())
        notifier.start()
        mask = EventsCodes.ALL_FLAGS['IN_CREATE']
        wdd = wm.add_watch(self.parentDataDir, mask, rec=False)
        


    
    def updateAll(self):
#        self.findLatestFrame()
#        self.loadLatestFrame()

#        self.determinezScaling()
        self.updatezScaling()

        self.displayRaw2D()
        #print self.imagetype
        if(idebugGUI): print "updating All"

        if(self.imagetype=='science'):   
            self.profile1D()
            self.plotProfile()
            self.updateExtractionSliders()
        
            self.doSkySubtraction()
            self.displaySS2D()
        else:
            self.data2DSkySub=self.data2D
#            self.data2DSkySub=np.copy(self.data2D)*0.0
            #self.plotProfile()
            self.medsky1d=0.0
            self.medskylev=0.0
            self.displaySS2D()
#            self.skysub2DView.canvas.clear()
        
        self.extract1D()
        self.display1D()
    

       
        #***
    """
    def updateAll(self):
        self.updatezScaling()
        self.displayRaw2D()
##        self.updateAll()
        if self.imagetype=='science' : self.doScience()
        if self.imagetype=='arc' : self.doArc()
        if self.imagetype=='other': self.doOther()
    """
    
    # -- This is kinda dumb but works!:
    def monitorFileSize(self):
        # watch size of newly created file and return when it stops increasing
        lastSize=0.0
        currSize = os.path.getsize(self.latestFrame)
#        print currSize
        print "monitoring read out..."
        while(currSize>lastSize):
#            print currSize
            lastSize=currSize
            time.sleep(1.0)
            currSize = os.path.getsize(self.latestFrame)
            #pass
#        print self.latestFrame+" finished reading out."
        print "finished reading out."
        currSize=0.0 # reset
   
#    class PTmp(ProcessEvent):
    def process_IN_CREATE(self, event):
    #def processMonitorData(self, event):
        newfile=os.path.join(event.path, event.name)
        # mkae sure .fit* file
        ext=newfile.split('.')[-1]
        try:
            if (ext[0:3] <>'fit'):
                return
        except:
            return
        
#        print "%s found: %s" % (time.ctime(), os.path.join(event.path, event.name))
        print "%s found: %s" % (time.ctime(), newfile)
        
        #print "waiting %.1f sec"%(self.ReadOutDelay)
        #time.sleep(self.ReadOutDelay) # wait for copy to finish??
        

        # this is replacement for findLatestFrame:
        self.latestFrame=self.parentDataDir+event.name
#        print time.ctime()

        self.monitorFileSize()
        
        self.loadLatestFrame()
#        self.updateAll()
        

        self.updatezScaling()
        self.displayRaw2D()
##        self.updateAll()
        if self.imagetype=='science' : self.doScience()
        if self.imagetype=='arc' : self.doArc()
        if self.imagetype=='other': self.doOther()
        
        #--
#        self.extract1D()
#        self.display1D()
    """
    def process_IN_CLOSE_NOWRITE(self, event):
    #def processMonitorData(self, event):
        print "%s closed: %s" % (time.ctime(), os.path.join(event.path, event.name))
        #print "waiting %.1f sec"%(self.ReadOutDelay)
        #time.sleep(self.ReadOutDelay) # wait for copy to finish??
        
        # this is replacement for findLatestFrame:
        self.latestFrame=self.parentDataDir+event.name
#        print time.ctime()

        self.loadLatestFrame()
#        self.updateAll()
        

        self.updatezScaling()
        self.displayRaw2D()
##        self.updateAll()
        if self.imagetype=='science' : self.doScience()
        if self.imagetype=='arc' : self.doArc()
        if self.imagetype=='other': self.doOther()
   """


    def loadLatestFrame(self):
        try:
            tmpdata2D,hdr = pyfits.getdata(self.latestFrame,header=True)
        except:
            time.sleep(1) # seems too quick and file is not finished writing??***
            tmpdata2D,hdr = pyfits.getdata(self.latestFrame,header=True)
        self.data2D=tmpdata2D[:,::-1] # flip so blue on left as God intended
        self.currFrameFill.setText(self.latestFrame)
        
        # -- determine object type
        self.imagetype='other'
        if(hdr['IMAGETYP']=='object'): self.imagetype='science'
        if(hdr['IMAGETYP']=='COMPARISON'): self.imagetype='arc'

        # if arc, get lamp name from header
        if(self.imagetype=='arc'): 
            self.lamp=hdr['ARC-LAMP']
            # unset WLC flag
            self.wlc=False
            self.wavelen=np.arange(self.data2D.shape[1])
        
        # get grating name (eventually from header, for now from combo box)
        self.grating=self.gratingCombo.currentText()
        print "[%s] [%s] [%ss]    [pipeline recipe: %s]" % (hdr['imagetyp'], self.grating, hdr['exptime'], self.imagetype)
        if(self.imagetype=='arc'): print self.lamp

#        print self.imagetype
#        if(self.imagetype<>'object' ): self.data2DSkySub=self.data2D



    def updatezScaling(self):
        if(str(self.ScalingCombo.currentText())=='ZScale'): 
            self.z0,self.z1=zscale(self.data2D)
        if(str(self.ScalingCombo.currentText())=='Min/Max'): 
            self.z0=np.min(self.data2D)
            self.z1=np.max(self.data2D)
        if(str(self.ScalingCombo.currentText())=='Manual'): 
            self.z0=float(self.z0Fill.text())
            self.z1=float(self.z1Fill.text())
        
        if (idebugGUI): print "updating scaling"
        self.z0Fill.setText("%.3f"%(self.z0))
        self.z1Fill.setText("%.3f"%(self.z1))
        self.displayRaw2D()
        try: 
            self.displaySS2D()
            self.plotProfile()
        except:
            pass
        
 
# http://matplotlib.1069221.n5.nabble.com/key-press-events-in-matplotlib-embedded-in-pyqt4-td27958.html
    def displayRaw2D(self):
        #foo=pyfits.getdata('/home/gilbank/Proj/CassSpect/data/Deatrick/a3000142.fits')
        self.raw2DGraph.canvas.ax.clear()
#        self.raw2DGraph.canvas.ax.imshow(self.data2D, interpolation='nearest',vmin=278,vmax=356)
        #self.raw2DGraph.canvas.ax.imshow(self.data2D, interpolation='nearest',vmin=self.z0,vmax=self.z1)
        #self.raw2DGraph.canvas.ax.imshow(self.data2D, interpolation='nearest',vmin=self.z0,vmax=self.z1,origin='lower',\
        self.raw2DGraph.canvas.ax.imshow(self.data2D, interpolation='nearest',vmin=float(self.z0Fill.text()),vmax=float(self.z1Fill.text()),origin='lower',\
                                         cmap='gist_gray')
                                         #cmap=matplotlib.cm("gist_gray"))
   #     self.raw2DGraph.canvas.ax.axhline(self.extractWin_x0,color='y',alpha=0.3)
        #self.profileGraph.canvas.ax.axis('tight')
        #self.profileGraph.canvas.ax.set_aspect('equal', 'datalim')
        if (self.showEWCheck.isChecked()):
#        if ( (self.showEWCheck.isChecked()) & (self.imagetype=='science') ):
            self.raw2DGraph.canvas.ax.axhspan(self.extractWin_x0,self.extractWin_x1,color='g',alpha=0.3)
        if ( (self.showSBCheck.isChecked()) & (self.imagetype=='science') ):
            if (self.useSB1Check.isChecked()):
                self.raw2DGraph.canvas.ax.axhspan(self.SB1_x0,self.SB1_x1,color='r',alpha=0.2)
            else:
                self.raw2DGraph.canvas.ax.axhspan(self.SB1_x0,self.SB1_x1,color='r',alpha=0.)
            if (self.useSB2Check.isChecked()):
                self.raw2DGraph.canvas.ax.axhspan(self.SB2_x0,self.SB2_x1,color='r',alpha=0.2)
            else:
                self.raw2DGraph.canvas.ax.axhspan(self.SB2_x0,self.SB2_x1,color='r',alpha=0.)

        # draw original plot again to fix axhspan messing with plot limits!
        self.raw2DGraph.canvas.ax.imshow(self.data2D, interpolation='nearest',vmin=float(self.z0Fill.text()),vmax=float(self.z1Fill.text()),origin='lower',\
                                         cmap='gist_gray')
        
        #self.raw2DGraph.canvas.draw()
#        self.raw2DGraph.canvas.ax.axhspan(40,80,color='r',alpha=0.3)
        #self.raw2DGraph.canvas.ax.
        self.raw2DGraph.canvas.draw()



    def doScience(self):
        #self.loadLatestFrame()
        #self.updatezScaling()
        #self.displayRaw2D()
        self.updateAll()

        self.profile1D()
        self.plotProfile()
        self.updateExtractionSliders()
    
        self.doSkySubtraction()
        self.displaySS2D()


    def doArc(self):
        #self.loadLatestFrame()
        #self.updatezScaling()
        #self.displayRaw2D()
        self.plotProfile()        
        #*** clear unnec. plots
        #*** set extraction pars
        
        self.doWLC()
        
        
    def doOther(self): 
        # bias, flats, ...
        pass
        


    def doWLC(self):
        arcspec=self.spec1d
        # library directory for ref file lists. Must be writable for temp files:
        
        libdir='/home/gilbank/Proj/CassSpect/data/libs/'
        #refArcFile=libdir+'deatrick_arc_intv.txt'
        refArcFile=libdir+'ref_%s_%s.txt'%(self.grating,self.lamp)
#        refArcFile='/home/gilbank/Proj/CassSpect/data/Deatrick/gr4_4_CuAr_tweakdc.txt'
        try:
            reflam,refspec = np.loadtxt(refArcFile,unpack='true')
        except:
            print "%s not found!"%(refArcFile)
            return
        
        
        linelist=libdir+'CuAr.lis'
        xp,lam0,lam1=np.loadtxt(linelist,usecols=(0,1,2),unpack='t')
        lines=lam0

        #arclis()
        #*** need to make _tmparc.lis
        print 'Attempting wavelength calibration (WLC)'
        try:
            #**** need to set up some options for pars with diff gratings (/lamps?)    
            if(self.grating=='gr4'):
                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=800.,maxFineShift=50.,iplot=False,idebug=False)
#                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=800.,maxFineShift=50.,iplot=False,idebug=True)##***
            if(self.grating=='gr7'):
                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=500.,maxFineShift=10.,iplot=False,idebug=False)
            
            
#            foo=interactiveWLC(lamcal,arcspec,lines)
        except:
            print 'WLC failed!!!'
            self.wlc=False
            self.wavlen=np.arange(self.data2D.shape[0])
            return
        # need to check if lamcal successful
        
        # -- update solution
        self.wavelen=lamcal
        self.lastWLCFile=str(self.currFrameFill.text())
        self.wlc=True
        self.currArcFill.setText(self.lastWLCFile)
        print 'WLC updated'
        self.display1D()
        
        # overplot, ref arc:
        close('all') #*** close all silly plot wins
        self.spec1DView.canvas.ax.plot(reflam,refspec,'g-',alpha=0.3)
        self.spec1DView.canvas.ax.set_xlabel(r'$\lambda (\AA)$')

        for line in lines:
            if ( (line > np.min(self.wavelen)) & (line < np.max(self.wavelen)) ):
                self.spec1DView.canvas.ax.axvline(line,color='k',alpha=0.2)
            else: pass

        self.spec1DView.canvas.draw()
        

    def plotProfile(self):
        self.profileGraph.canvas.ax.clear()
        xx=np.arange(self.data2D.shape[0])
        self.profileGraph.canvas.ax.plot(xx,np.median(self.data2D,1),'b-')
        if self.imagetype=='science':
            popt=self.objectGaussCoeffs
            self.profileGraph.canvas.ax.plot(xx,gaussian(xx,popt[0],popt[1],popt[2],popt[3]),'k--')
        self.profileGraph.canvas.ax.axvspan(self.extractWin_x0,self.extractWin_x1,color='g',alpha=0.4)
        #self.profileGraph.canvas.ax.axvline(self.extractWin_x0)
        #self.profileGraph.canvas.ax.axvline(self.extractWin_x1)
        if ( (self.useSB1Check.isChecked()) & (self.imagetype=='science') ):
            self.profileGraph.canvas.ax.axvspan(self.SB1_x0,self.SB1_x1,color='r',alpha=0.3)
        else:
            self.profileGraph.canvas.ax.axvspan(self.SB1_x0,self.SB1_x1,color='r',alpha=0.)
        if ( (self.useSB2Check.isChecked()) & (self.imagetype=='science') ):
            self.profileGraph.canvas.ax.axvspan(self.SB2_x0,self.SB2_x1,color='r',alpha=0.3)
        else:
            self.profileGraph.canvas.ax.axvspan(self.SB2_x0,self.SB2_x1,color='r',alpha=0.)
        #self.profileGraph.canvas.ax.axvline(self.SB1_x0)
        #self.profileGraph.canvas.ax.axvline(self.SB1_x1)
        #self.profileGraph.canvas.ax.axvline(self.SB2_x0)
        #self.profileGraph.canvas.ax.axvline(self.SB2_x1)
        self.profileGraph.canvas.ax.set_ylim([self.z0,self.z1])
        self.profileGraph.canvas.draw()
        
        self.updateExtractionSliders()



    def profile1D(self):
        xsec=np.median(self.data2D[:,:],1)
        pkmax,pkmin=peakdet(xsec,100)
        npks=np.shape(pkmax)[0]
        #*** only works for 1 peak, so far


        # fit Gaussian to position
        width=3.0 #***
        p0=[pkmax[0][1],pkmax[0][0],width,np.min(xsec)]
        x=np.arange(len(xsec))
        y=xsec

        e=np.ones(len(xsec))
        popt, pcov = curve_fit(gaussian, x, y, sigma=e)
        
        self.objectPkPos=pkmax[0][1]
        self.objectGaussCoeffs=popt
        
        #-- calculate extract,sideband windows. Need some error checking!
        extract_halfwidth=3.0*popt[2]
#        extract_halfwidth=np.min([extract_halfwidth,8.0])
        # Need proper error-checking on these extraction pars! ***
 
        
        SBoffset=10.0
        SBwidth=30.0
        
        if (self.autoExtraction):
            self.extractWin_x0=popt[1]-extract_halfwidth
            self.extractWin_x1=popt[1]+extract_halfwidth
            self.SB1_x1=self.extractWin_x0-SBoffset
            self.SB1_x0=self.SB1_x1-SBwidth
            self.SB2_x0=self.extractWin_x1+SBoffset
            self.SB2_x1=self.SB2_x0+SBwidth
        
        # -- update 2D:
        self.displayRaw2D()
        
        
    def updateAutoExtraction(self):
        self.autoExtraction=True
        self.updateAll()
        
    def readExtractionSliders(self): 
#        if (self.useSB1Check.isChecked()):
        self.SB1_x0=self.SB1Slider.lowerValue
        self.SB1_x1=self.SB1Slider.upperValue
#        else:
#            self.SB1_x0=None
#            self.SB1_x1=None            
        self.SB2_x0=self.SB2Slider.lowerValue
        self.SB2_x1=self.SB2Slider.upperValue
        self.extractWin_x0=self.extractWinSlider.lowerValue
        self.extractWin_x1=self.extractWinSlider.upperValue
        self.autoExtraction=False
        self.plotProfile()

    def updateExtractionSliders(self):
        
        self.extractWinSlider.setRange(0, 132)
        self.extractWinSlider.setSpan(self.extractWin_x0,self.extractWin_x1)
        self.SB1Slider.setRange(0, 132)
        self.SB1Slider.setSpan(self.SB1_x0,self.SB1_x1)
        self.SB2Slider.setRange(0, 132)
        self.SB2Slider.setSpan(self.SB2_x0,self.SB2_x1)
        self.displayRaw2D()
        self.doSkySubtraction()
        self.displaySS2D()
        self.extract1D()
        self.display1D()
        #self.updateAll()

    def doSkySubtraction(self):
        if(self.imagetype<>'science'): return
        # let's just do a simple 0th order fit to SBs for now:
        spec2d=self.data2D
        b0=self.SB1_x0
        b1=self.SB1_x1
        b2=self.SB2_x0
        b3=self.SB2_x1
        xx=np.arange(spec2d.shape[0])
        sky01=spec2d[b0:b1,:]
        sky23=spec2d[b2:b3,:]
        if (self.useSB1Check.isChecked() & self.useSB2Check.isChecked() ):
            sky03=np.append(sky01,sky23,axis=0)
        else:
            if (self.useSB1Check.isChecked()):
                sky03=sky01
            else:
                if (self.useSB2Check.isChecked()):
                    sky03=sky23
                else:
                    sky03=np.copy(sky01)*0.0
        medsky1d=np.median(sky03,0)
        medsky2d=np.transpose(np.repeat(np.reshape(medsky1d,(-1,1)),(spec2d.shape[0]),axis=1))
    
        #popt=self.objectGaussCoeffs
        #profile1d=gaussian(xx,popt[0],popt[1],popt[2],popt[3])
        #profile2d=np.repeat(np.reshape(profile1d,(-1,1)),(spec2d.shape[1]),axis=1)

        self.data2DSkySub = spec2d - medsky2d
        self.medsky1d=medsky1d
        self.medskylev=np.median(medsky1d)


    def displaySS2D(self):
        if (self.imagetype=='science'):
            #foo=pyfits.getdata('/home/gilbank/Proj/CassSpect/data/Deatrick/a3000142.fits')
            self.skysub2DView.canvas.ax.clear()
            self.skysub2DView.canvas.ax.imshow(self.data2DSkySub, interpolation='nearest',\
                                               vmin=self.z0-self.medskylev,vmax=self.z1-self.medskylev,origin='lower')
            self.skysub2DView.canvas.draw()
        else:
            self.skysub2DView.canvas.ax.clear()
            self.skysub2DView.canvas.draw()
       
    def extract1D(self):
        # simple box-sum extraction:
        if (self.imagetype=='science'):
            self.spec1d = np.sum(self.data2DSkySub[self.extractWin_x0:self.extractWin_x1,:],axis=0)
        else:
            self.spec1d = np.sum(self.data2D[self.extractWin_x0:self.extractWin_x1,:],axis=0)


    def display1D(self):
        self.spec1DView.canvas.ax.clear()
        xx=np.arange(self.spec1d.size)
#        if(self.lastWLCFile):
#            xx=self.wavelen

        if (self.wlc):
            xx=self.wavelen        
        #self.spec1DView.canvas.ax.plot(xx,self.spec1d,'b-')
        self.spec1DView.canvas.ax.plot(xx,self.spec1d,'b-')
        #medskyUnderAp = self.medsky1d*(self.extractWin_x1-self.extractWin_x0)
        if (self.imagetype=='science'):
            medskyUnderAp = ( self.medsky1d  - self.medskylev ) *(self.extractWin_x1-self.extractWin_x0)
            self.spec1DView.canvas.ax.plot(xx,medskyUnderAp,'r-')
        #self.spec1DView.canvas.ax.plot(xx,medskyUnderAp,'r-')
#        print np.shape(xx)
#        print np.shape(medskyUnderAp)
        self.spec1DView.canvas.ax.set_title(self.latestFrame)
        y0,y1=zscale(self.spec1d)
#        self.spec1DView.canvas.ax.set_ylim([y0,y1])
        if(self.imagetype=='arc'): y1=np.max(self.spec1d)
        if(y0<100.0): y0=0.0
        if(y1>100.0):
            self.spec1DView.canvas.ax.set_ylim([y0,y1])
        if (self.wlc):
            self.spec1DView.canvas.ax.set_xlabel(r'$\lambda (\AA)$')
        else:
            self.spec1DView.canvas.ax.set_xlabel('x (pix)')
        self.spec1DView.canvas.draw()

    """
    def OLDdisplayRaw2D(self):
        self.dpi = 100
#        self.fig = Figure( (11,1), dpi=self.dpi)
        self.fig = Figure(  dpi=self.dpi)
        self.fig.clear() 
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.raw2DGraph)
        #-- tmp hardwired stuff:
        foo=pyfits.getdata('/home/gilbank/Proj/CassSpect/data/Deatrick/a3000142.fits')
        self.axes = self.fig.add_subplot(111)
        self.axes.imshow(foo, interpolation='nearest',vmin=278,vmax=356)
        self.canvas.draw()
    """


    def showDialog(self):
#        print "foo"
        fname = QtGui.QFileDialog.getOpenFileName(None,'Open file', self.parentDataDir)
        
        f = open(fname, 'r')
        
        with f:        
            data = f.read()
            self.currFrameFill.setText(fname) 
            self.latestFrame=str(fname)
            self.loadLatestFrame()
            print "loaded "+self.latestFrame+' [%s] [%s]'%(self.imagetype,self.grating)
            self.updateAll()

    def showDialogSetDatadir(self):
#        print "foo"
        #fname = QtGui.QFileDialog.getOpenFileName(None,'Open file', self.parentDataDir)
        fname = str(QtGui.QFileDialog.getExistingDirectory(None, "Set directory where the telescope is writing data"))
        self.parentDataDir=fname+'/' 
        print "data directory is now: "+str(self.parentDataDir)


    def normalOutputWritten(self, text):
        """Append text to the QTextEdit."""
        # Maybe QTextEdit.append() works as well, but this is how I do it:
        cursor = self.logView.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.logView.setTextCursor(cursor)
        self.logView.ensureCursorVisible()
        
    def closeEvent(self):

        quit_msg = "Are you sure you want to exit the program?"
        #reply = QtGui.QMessageBox.question(self, 'Message', quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        reply = QtGui.QMessageBox.question(None, 'Removal', quit_msg, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if (reply == QtGui.QMessageBox.Yes):
#            QtGui.QApplication.quit() # Nope!
            QtGui.qApp.closeAllWindows()
        else:
            pass
        
