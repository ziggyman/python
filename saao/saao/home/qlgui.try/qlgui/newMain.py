#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time,os,sys,json,glob
from PyQt4 import QtCore, QtGui

#import logging as log
import logging, StringIO, time
import numpy as np
from wx import FrameNameStr
try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from astropy.io import fits
from zscale import zscale

from qlGUI import Ui_MainWindow

from profPane import ProfPane
from im2DPane import Im2DPane
from im1DPane import Im1DPane
from logPane import logPane,logBuffer

from mlabWidgets import *
from QxtSpanSlider import QxtSpanSlider

from specutils import wlc_arc,xcor2


from pyinotify import WatchManager, Notifier, ThreadedNotifier, EventsCodes, ProcessEvent



dbglvl=2



class QLGui(QtGui.QMainWindow, ProfPane, Im2DPane, Im1DPane, logPane):
    def __init__(self, parent=None, infile=None):
        QtGui.QWidget.__init__(self, parent)
        
        #set up the main UI
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.dbglvl = dbglvl

#        self.libdir='/home/gilbank/Proj/spup/GUI/libs/'
#        self.libdir='libs/'
        self.libdir = '/usr/local/lib/qlgui/'

        self.spup = True # else work with old spectrograph data


        self.ui.logView.setReadOnly(True)
        # Install the custom output stream
        sys.stdout = EmittingStream(textWritten=self.normalOutputWritten)


        # set some default size values. use to try to set sliders, etc. to sensible ranges:
        if(self.spup):
            self.sx = 2148
            self.sy = 256
        else:
            self.sx = 2098
            self.sy = 512

        self.xbin=1
        self.ybin=1



        # -- update sliders to double sliders:
        geom = self.ui.SB1Slider.geometry()
        self.ui.SB1Slider = QxtSpanSlider(self.ui.centralWidget)
        self.ui.SB1Slider.setGeometry(geom)
        #self.ui.SB1Slider.setOrientation(QtCore.Qt.Horizontal)
        """
        self.SB1Slider = QxtSpanSlider(self.centralWidget)
        self.SB1Slider.setGeometry(QtCore.QRect(1100, 520, 160, 23))
        self.SB1Slider.setOrientation(QtCore.Qt.Horizontal)
        self.SB1Slider.setObjectName(_fromUtf8("SB1Slider"))
        """

        if(self.spup):
            self.ui.SB1Slider.setRange(0,132)#***
            self.ui.SB1Slider.setSpan(80,120)#***
        else:
            self.ui.SB1Slider.setRange(0,132)
            self.ui.SB1Slider.setSpan(80,120)
        # --

        geom = self.ui.SB2Slider.geometry()
        self.ui.SB2Slider = QxtSpanSlider(self.ui.centralWidget)
        self.ui.SB2Slider.setGeometry(geom)
        #self.ui.SB1Slider.setOrientation(QtCore.Qt.Horizontal)
        """
        self.SB1Slider = QxtSpanSlider(self.centralWidget)
        self.SB1Slider.setGeometry(QtCore.QRect(1100, 520, 160, 23))
        self.SB1Slider.setOrientation(QtCore.Qt.Horizontal)
        self.SB1Slider.setObjectName(_fromUtf8("SB1Slider"))
        """
        if(self.spup):
            self.ui.SB2Slider.setRange(0,132) #***
            self.ui.SB2Slider.setSpan(50,10) #***
        else:
            self.ui.SB2Slider.setRange(0,132)
            self.ui.SB2Slider.setSpan(50,10)
        
        geom = self.ui.extractWinSlider.geometry()
        self.ui.extractWinSlider = QxtSpanSlider(self.ui.centralWidget)
        self.ui.extractWinSlider.setGeometry(geom)

        
        self.extractWin_x0=self.ui.extractWinSlider.lowerValue
        self.extractWin_x1=self.ui.extractWinSlider.upperValue
        self.SB1_x0=self.ui.SB1Slider.lowerValue
        self.SB1_x1=self.ui.SB1Slider.upperValue
        self.SB2_x0=self.ui.SB2Slider.lowerValue
        self.SB2_x1=self.ui.SB2Slider.upperValue



        self.ui.lineEditz0.setText('0.')
        self.ui.lineEditz1.setText('10000.')


        # header keyword:
        if(self.spup):
            self.imtype='EXPTYPE'
            self.sciKW = 'SCIENCE' # keyword for science
#            self.arcKW='arc'
            self.arcKW='ARC'
        else:
            self.imtype='IMAGETYP'
            self.sciKW = 'object' # keyword for science
            self.arcKW='COMPARISON'

        pars = {}
        pars['smoothPix'] = 0
        pars['SkySub'] = False
        pars['useSB1'] = True
        pars['useSB2'] = True
        pars['LockDisplays'] = False
#        pars['showExtWin'] = True
#        pars['showSideBands'] = True
#        pars['dataDir'] = '/home/gilbank/Proj/CassSpect/data/Deatrick/'
#        pars['dataDir'] = '/home/ccd/cassspectr/cassspectr/'
        f=open('/home/ccd/.spup/datadir.json','r')
        pars['dataDir'] = json.load(f)
        f.close()
#        pars['outDataDir'] = '/home/ccd/cassspectr/cassspectr/'
        pars['outDataDir'] = pars['dataDir']
        pars['currFrame'] = 'none'
        pars['obsType'] = 'none'
        pars['currArc'] = 'none'
        pars['currFlat'] = 'none'
        pars['currFluxCal'] = 'none'
        pars['currBias'] = 'none'
        pars['scaling'] = 'ZScale'
        pars['z0'] = np.NAN
        pars['z1'] = np.NAN
        if(self.spup):
            pars['extractWin_x0'] = 120.0 #**
            pars['extractWin_x1'] = 140.0 #**
        else:
            pars['extractWin_x0'] = 50.0
            pars['extractWin_x1'] = 70.0
        pars['SBoffset'] = 5
        pars['SBwidth'] = 15
        
        self.extractWin_x0 = pars['extractWin_x0'] ; self.extractWin_x1 = pars['extractWin_x1'] 
        
        
        pars['autoExtraction'] = True

#        pars['isScience']=True

        self.pars = pars


        # Decide extraction type:




        self.ui.lineEditDataDir.setText(self.pars['dataDir'])
        self.ui.lineEditOutDataDir.setText(self.pars['outDataDir'])


#        self.imagetype='science'
        self.imagetype='none'
        self.wlc=False
        
        # -- Initialise values:
 
        self.autoExtraction=True
        self.lamp=None
        self.grating=None

        self.parentDataDir=None #'/home/gilbank/Proj/CassSpect/data/Deatrick/'
        self.imagetype=""
        self.lastWLCFile=None
        self.wlc=False
        self.latestFrame=None ##self.parentDataDir+'a3000142.fits'

        self.scaledImage=None
        
        
        
        # disable unimplemented options:
        self.ui.useSB1Check.setEnabled(False)
        self.ui.useSB2Check.setEnabled(False)
        self.ui.checkBoxLockDisplays.setEnabled(False)
        self.ui.checkBoxWrite1DFITS.setEnabled(False)
        self.ui.pushButCurrArc.setEnabled(False)

        
        
        # -- Connect signals and slots
        QtCore.QObject.connect(self.ui.actionQuit, QtCore.SIGNAL(_fromUtf8("activated()")), self.closeEvent)
        QtCore.QObject.connect(self.ui.pushButPan, QtCore.SIGNAL(_fromUtf8("pressed()")), self.doPanAll)        
        QtCore.QObject.connect(self.ui.pushButZoom, QtCore.SIGNAL(_fromUtf8("pressed()")), self.doZoomAll)
        QtCore.QObject.connect(self.ui.pushButResetView, QtCore.SIGNAL(_fromUtf8("pressed()")), self.resetZoomAll)
        QtCore.QObject.connect(self.ui.pushButSmoothObj, QtCore.SIGNAL(_fromUtf8("pressed()")), self.display1D)
        QtCore.QObject.connect(self.ui.pushButNoSmooth, QtCore.SIGNAL(_fromUtf8("pressed()")), self.noSmoothing)
        QtCore.QObject.connect(self.ui.pushButCurrFrame, QtCore.SIGNAL(_fromUtf8("pressed()")), self.loadFrame)
#        QtCore.QObject.connect(self.ui.pushButCurrFrame, QtCore.SIGNAL(_fromUtf8("pressed()")), self.loadAndUpdate)
        QtCore.QObject.connect(self.ui.checkBoxShowExtractWin, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.display2D)
        QtCore.QObject.connect(self.ui.checkBoxShowSB, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.display2D)
        QtCore.QObject.connect(self.ui.zScaleUpdateBut, QtCore.SIGNAL(_fromUtf8("clicked()")), self.updatezScaling)
        QtCore.QObject.connect(self.ui.pushButAutoExtract, QtCore.SIGNAL(_fromUtf8("clicked()")), self.updateAutoExtraction)
        QtCore.QObject.connect(self.ui.useSB1Check, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.updateAll)
        QtCore.QObject.connect(self.ui.useSB2Check, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.updateAll)
        QtCore.QObject.connect(self.ui.SB1Slider, QtCore.SIGNAL(_fromUtf8("sliderReleased()")), self.readExtractionSliders)
        QtCore.QObject.connect(self.ui.SB2Slider, QtCore.SIGNAL(_fromUtf8("sliderReleased()")), self.readExtractionSliders)
        QtCore.QObject.connect(self.ui.extractWinSlider, QtCore.SIGNAL(_fromUtf8("sliderReleased()")), self.readExtractionSliders)
        QtCore.QObject.connect(self.ui.comboBoxColorMap, QtCore.SIGNAL(_fromUtf8("currentIndexChanged(int)")), self.display2D)
        QtCore.QObject.connect(self.ui.pushButChooseDataDir, QtCore.SIGNAL(_fromUtf8("clicked()")), self.showDialogSetDatadir)
        QtCore.QObject.connect(self.ui.pushButChooseOutDataDir, QtCore.SIGNAL(_fromUtf8("pressed()")), self.showDialogSetOutDatadir)
        QtCore.QObject.connect(self.ui.pushButCurrFrameBrowse, QtCore.SIGNAL(_fromUtf8("pressed()")), self.showDialogSetCurrFrame)

        QtCore.QObject.connect(self.ui.pushButIgnoreWLC, QtCore.SIGNAL(_fromUtf8("clicked()")), self.ignoreWLC)
#        QtCore.QObject.connect(self.ui.pushButIgnoreWLC, QtCore.SIGNAL(_fromUtf8("clicked()")), self.ignoreWLC)
        QtCore.QObject.connect(self.ui.pushButIgnoreWLClibArc, QtCore.SIGNAL(_fromUtf8("clicked()")), self.useLibArc)

        QtCore.QObject.connect(self.ui.useSB1Check, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.updateAll)
        QtCore.QObject.connect(self.ui.useSB2Check, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.updateAll)
        

        
        # ---- Set up logging
        # see: https://docs.python.org/2/howto/logging-cookbook.html
        
        logfile = 'QL%s.log'%(time.strftime("%Y%m%d", time.gmtime()))
        logging.basicConfig(filename=logfile,level=logging.DEBUG,format='%(asctime)s %(message)s')
#        logFormatter = logging.Formatter('%(asctime)s   %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M:%S')
        logFormatter = logging.Formatter('%(asctime)s   %(levelname)-8s %(message)s', datefmt='%H:%M:%S')
        self.logBuffer = logBuffer()
        self.logBuffer.bufferMessage.connect(self.on_logBuffer_bufferMessage)

###        logFormatter = logging.Formatter('%(levelname)s: %(message)s')

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
        self.makeDummy()


        #-- Monitor data directory for new images        
        wm = WatchManager()
        notifier = ThreadedNotifier(wm, self.process_IN_CREATE)
        notifier.start()
        self.notifier = notifier
        mask = EventsCodes.ALL_FLAGS['IN_CREATE']
        wdd = wm.add_watch(self.pars['dataDir'], mask, rec=False)
        #--




        # read values from GUI:
        self.pars['showSideBands'] = self.ui.checkBoxShowSB.isChecked()
        self.pars['smoothPix'] = self.ui.spinBoxSmoothWidth.value()
                
        self.logger.debug(self.ui.checkBoxShowSB.isChecked())
        
        
        
        
        
        
        
            
        
#        self.profileGraph = matplotlibWidgetwBar(self.centralWidget)
#        self.profPlot = matplotlibWidget(self.ui.centralWidget)
###        sc = matplotlibWidgetwBar(self.ui.profileGraph)
        # **** clunky, but effective!
        #profWin = matplotlibWidgetwBar(self.ui.centralWidget)
        profWin = matplotlibWidget(self.ui.centralWidget)
        #self.profPlot.setGeometry(QtCore.QRect(1100, 140, 330, 271))
        self.profWin = profWin
        
#        imname = 'a3000134.fits'
#        self.logger.info(imname)
        
#        self.image = self.pars['dataDir']+imname
#        data2D,hdr = fits.getdata(self.image, header=True)
#        self.data2D = data2D
#        self.pars['currFrame'] = imname
#        self.pars['obsType'] = hdr[self.imtype]
#        print self.pars['obsType']


#        self.ui.lineEditCurrFrame.setText(pars['currFrame'])


        self.plotProfile()
        
        #self.profWin.mpl_toolbar.zoom()
        


        #graph2D = matplotlibWidgetwBar(self.ui.centralWidget)
        graph2D = matplotlibWidget(self.ui.centralWidget)
        self.graph2D = graph2D
        
        self.display2D()

        specWin = matplotlibWidget(self.ui.centralWidget)
        self.specWin = specWin
        
#        self.logger.debug('**',self.sky1d)
        self.display1D()


        """
        self.drawIndicators()


        print self.ax0
        print self.ax1
        print self.ax2
        """




    def doScience(self):
        #self.loadLatestFrame()
        #self.updatezScaling()
        #self.displayRaw2D()
#        self.subZero()
        self.updateAll()

#        self.profile1D()
        """
        self.fitSky()
        self.extract1D()
        self.plotProfile()
        self.updateExtractionSliders()
    
#        self.doSkySubtraction()
        self.display2D()
        """

    def findNearestRefArc(self):
###        print self.lamp, self.grating, self.grangle


        try:
            reffiles = glob.glob(self.libdir+self.lamp+'/'+self.grating+'_*_dc.txt')
            self.reffile=''
        except: return
        bestgrang = -99.0
        for f in reffiles:
#            print f
#            print f.split('_')
#            print f.split('_')[-2]

            tgrang = f.split('_')[-2]
#            print tgrang
            tgrang = float(tgrang)
            if ( np.abs(tgrang-float(self.grangle))< np.abs(bestgrang-float(self.grangle)) ):
                bestgrang = tgrang

        reffile = self.libdir+self.lamp+'/'+self.grating+'_%.2f'%(bestgrang)+'_dc.txt'
        print reffile
        self.reffile = reffile


    def showRefArc(self):
        # display closest arc from reference library -
        # same lamp, same grating, closest gr-angle:
        #
        self.findNearestRefArc()
        try:
            x,y = np.loadtxt(self.reffile,unpack=True)
        except:
            print 'no available reference arc for this setup'
            return

        if(self.wlc==False):
            self.logger.info('WLC failed - using nearest available reference arc\n%s'%self.reffile)
            wlcapprox=True
            self.wlc=True
            self.wavelen = x
        else:
            wlcapprox=False
        self.display1D()

#        print x,y
        # plot:
        specWin = self.specWin

        if (wlcapprox):
#            print 'APPROX WLC!'
            specWin.canvas.ax.set_xlabel(r'$\lambda(\AA)$ APPROX! - BLUE SPECTRUM IS ACCURATE')


        geom = self.ui.spec1DView.geometry()
        specWin.setGeometry(geom)
        specWin.canvas.ax.plot(x,y,'b-',label=self.reffile,alpha=0.7)
        specWin.canvas.ax.legend(bbox_to_anchor=(1.12, 1.0))
        specWin.canvas.draw()
        specWin.show()


    def useLibArc(self):
        # ignore WLC (might be wrong) and use nearest library arc instead:
        self.wlc=False
        self.showRefArc()


    def doArc(self):
        #self.loadLatestFrame()
        self.updatezScaling()
        self.display2D()
        #self.updateAll()
        self.extract1D()
        self.plotProfile()  
        #self.updateExtractionSliders()
        #*** clear unnec. plots
        #*** set extraction pars
        self.display1D()
#        if not self.spup:
#            self.grating = self.ui.comboBoxGrating.currentText()
#        else:
        self.grating = self.hdr['grating']
        self.grangle = self.hdr['GR-ANGLE']
        self.lamp = self.hdr['arc-lamp']
        self.doWLC()

        #**
#        self.wlc=False
        #**

        if self.wlc==False:
            self.showRefArc()

        self.checkOutputReq()
        self.writeResults()

        
    def doOther(self): 
        # bias, flats, ...
        self.updatezScaling()
        self.display2D()
        #self.updateAll()
        self.extract1D()
        self.plotProfile()
        #self.updateExtractionSliders()
        #*** clear unnec. plots
        #*** set extraction pars
        self.display1D()

#        pass


    def doHartmann(self):
        self.updatezScaling()
        self.display2D()
        self.extract1D()
        self.plotProfile()
        self.display1D()

        # determine where we are in the Hartmann sequence:
        shutterPos = self.hdr['HARTPOS']
        hartSeqNumString = self.hdr['HARTSEQ']
        # Exposure 1/3
        hartExpNum = int(hartSeqNumString.split('/')[0])
        hartNExps = int(hartSeqNumString.split('/')[1])
        focusPos = self.hdr['FOCUSPOS']

        i = hartExpNum-1 # assuming 1-based in ui***
        print('Hartmann seq: %i/%i; %s; Focus=%.3f'%(hartExpNum,hartNExps,shutterPos,focusPos))
#        print self.pars['currFrame']

        # keep track in array/file of:
        # focuspos, A file number, B file number, <FWHM> and cross-corelation between A/B frames

        # Reset for new run:
        # Need to be careful, do we always start with shutterPos B at exp=1??
        if ( (hartExpNum==1) & (shutterPos=='B') ):
            # New sequence started:
            self.hartFocusPos = np.zeros(hartNExps)
            self.hartAfile = list(np.repeat('',hartNExps))
            self.hartBfile = list(np.repeat('',hartNExps))
#            self.hartAspec1d = list(np.repeat('',hartNExps))
#            self.hartBspec1d = list(np.repeat('',hartNExps))
            self.hartFWHM_A = np.zeros(hartNExps)
            self.hartFWHM_B = np.zeros(hartNExps)
            self.hartXcor = np.zeros(hartNExps)
            self.hartFocusPos = np.zeros(hartNExps)
            self.hartXcorShift = np.zeros(hartNExps)
        self.hartFocusPos[i] = focusPos
        if (shutterPos=='A'):
            self.hartAfile[i] = str(self.pars['currFrame'])
            self.hartAspec1d = self.obj1d
            self.hartFWHM_A[i] = self.getFWHMframe()
        elif (shutterPos=='B'):
            self.hartBfile[i] = str(self.pars['currFrame'])
            self.hartBspec1d = self.obj1d
            self.hartFWHM_A[i] = self.getFWHMframe()
        else:
            print 'Not in Hartmann Seq.?!?!'
            return
        print self.hartAfile
        print self.hartBfile

        if( (self.hartAfile[i]<>'') & (self.hartBfile[i]<>'') ):
            print 'Have both A & B images for this position...\nDo Xcor!'
            # check if we have both A & B exposures at this position, if so, run XCor:***

            # xcor A with B, we only care about size of shift
            maxshift = 20 # pixels

            # Caution: don't want to include bias level. Make sure we have real data. come in by 100 pix each side?
            #*** binning?
            inspec = self.hartAspec1d[100:-100]
            refspec = self.hartBspec1d[100:-100]
            shft,CC = xcor2(inspec,refspec,maxshift,dx=1.,quiet='false')
            #print 'measured shift = %.2f; CC = %.2f'%(shft,CC)
            print 'measured shift = %.2f'%(shft)
            self.thisShiftA = shft
            self.hartXcorShift[i] = shft

            self.display1DHart()

            if (hartExpNum==hartNExps):
                print 'final image in sequence\nDo analysis!'
                self.doAnalyseXcor()


        else:
            print 'waiting for both A,B images at this position'


    def doAnalyseXcor(self):
        print 'focus  best-fit shift'
        print '====================='
        for i in range(len(self.hartFocusPos)):
            print self.hartFocusPos[i], self.hartXcorShift[i]

#        pass


    def getFWHMframe(self):
        return -1.0

        
    def doWLC(self):
        print 'doWLC()'
        arcspec=self.obj1d
        # library directory for ref file lists. Must be writable for temp files:
        
#        libdir='/home/gilbank/Proj/CassSpect/data/libs/'
        #libdir='./libs/'
        #refArcFile=libdir+'deatrick_arc_intv.txt'
#        refArcFile=libdir+'ref_%s_%s.txt'%(self.grating,self.lamp)
        self.findNearestRefArc()
        refArcFile = self.reffile
#        refArcFile='/home/gilbank/Proj/CassSpect/data/Deatrick/gr4_4_CuAr_tweakdc.txt'
        try:
            reflam,refspec = np.loadtxt(refArcFile,unpack='true')
        except:
            print "%s not found!"%(refArcFile)
            return
        
        
#        linelist=libdir+'CuAr.lis'
        linelist=self.libdir+'%s.lis'%self.lamp
        try:
            xp,lam0,lam1=np.loadtxt(linelist,usecols=(0,1,2),unpack='t')
            lines=lam0
        except:
            lines = np.loadtxt(linelist) #*** just list of wavelengths of lines

#        print '**',lines

        #arclis()
        #*** need to make _tmparc.lis
        print 'Attempting wavelength calibration (WLC)'
        try:
            # mute output:
            sys.stdout = sys.__stdout__ # send back to terminal
            
            #**** need to set up some options for pars with diff gratings (/lamps?)    
            if(self.grating=='gr4'):
                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=800.,maxFineShift=50.,iplot=False,idebug=False)
#                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=800.,maxFineShift=50.,iplot=False,idebug=True)##***
            elif(self.grating=='gr7'):
##                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=500.,maxFineShift=10.,iplot=False,idebug=False)
                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=1200.,maxFineShift=10.,iplot=False,idebug=False)
#                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=500.,maxFineShift=10.,iplot=True,idebug=True)
#                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=1200.,maxFineShift=50.,iplot=False,idebug=False)

            elif(self.grating=='gr5'):
##                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=500.,maxFineShift=10.,iplot=False,idebug=False)
                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=1200.,maxFineShift=50.,iplot=False,idebug=False)
            else:
                lamcal = wlc_arc(arcspec,reflam,refspec,linelist=linelist,maxCoarseShift=1200.,maxFineShift=10.,iplot=False,idebug=False)

            # set back to log window
            sys.stdout = EmittingStream(textWritten=self.normalOutputWritten)
            
#            foo=interactiveWLC(lamcal,arcspec,lines)
        except:
            print 'WLC failed!!!'
            self.wlc=False
            self.wavlen=np.arange(self.data2D.shape[0])
            return
        # need to check if lamcal successful
        
        # -- update solution
        self.wavelen=lamcal
        self.lastWLCFile=str(self.ui.lineEditCurrFrame.text())
        self.wlc=True
        self.ui.lineEditCurrArc.setText(self.lastWLCFile)
        print 'WLC updated'
        self.display1D()
        
        # overplot, ref arc:
        #close('all') #*** close all silly plot wins
        # scale ref spec for diplay:
        srefspec = refspec - arcspec.min()
        srefspec = srefspec*arcspec.max()/srefspec.max()

        self.specWin.canvas.ax.plot(reflam,srefspec,'g-',alpha=0.3)
        self.specWin.canvas.ax.set_xlabel(r'$\lambda (\AA)$')

        for line in lines:
            if ( (line > np.min(self.wavelen)) & (line < np.max(self.wavelen)) ):
                self.specWin.canvas.ax.axvline(line,color='k',alpha=0.2)
                #print line
            else: pass

        self.specWin.canvas.draw()
        
    def ignoreWLC(self):
        self.wlc = False
        self.display1D()
        #**
        # clear "arc used" box:
        self.ui.lineEditCurrArc.setText('')

    def makeDummy(self):
        # make a dummy image and plots when none is specified:
        self.data2D = np.zeros((self.sy,self.sx))
        self.sky1d = np.zeros(self.sx)
        self.obj1d = np.zeros(self.sx)
        self.imagetype='none'
        self.currFrame='none'
        

    def drawIndicators(self):
        
        pars = self.pars
        
        #self.ui.lineEditCurrFrame.setText(pars['currFrame'])







    def cfRange(self):
        # compare ranges on various plots, if "lock displays" is True
        
        # coords, x,y,flux [where x may be pixels or lambda!]
        
        # make some assumptions:
        # flux 2D plot always overrides flux profWin
        
        # profile:
        ax0y = np.array(self.ax0[2,3])
        ax0flux = np.array(self.ax0[0,1])
        # 1D:
        ax1x = np.array(self.ax1[0,1])
        ax1flux = np.array(self.ax1[2,3])
        #2D: 
        ax2x = np.array(self.ax2[0,1])
        ax2y = np.array(self.ax2[2,3])
        
        
        


    def noSmoothing(self):
        # undo smoothing:
        self.ui.spinBoxSmoothWidth.setValue(0.0)
        self.pars['smoothPix'] = 0
        self.display1D()


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



    def doPanAll(self):
        self.profWin.mpl_toolbar.pan()
        self.graph2D.mpl_toolbar.pan()
        self.specWin.mpl_toolbar.pan()

    def doZoomAll(self):
        self.profWin.mpl_toolbar.zoom()
        self.graph2D.mpl_toolbar.zoom()
        self.specWin.mpl_toolbar.zoom()

    def resetZoomAll(self):
        self.profWin.mpl_toolbar.home()
        self.graph2D.mpl_toolbar.home()
        self.specWin.mpl_toolbar.home()



    def normalOutputWritten(self, text):
        """Append text to the QTextEdit."""
        # Maybe QTextEdit.append() works as well, but this is how I do it:
        cursor = self.ui.logView.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.ui.logView.setTextCursor(cursor)
        self.ui.logView.ensureCursorVisible()


    def determineObsType(self):
        # 
#        if (self.pars['obsType']==self.sciKW):
#            self.pars['isScience']=True
#        else:
#            self.pars['isScience']=False
            # -- determine object type
        self.imagetype='other'
        #*** watch case-sensitivity later!***
        if(self.hdr[self.imtype]==self.sciKW):
            self.imagetype='science'
        elif(self.hdr[self.imtype]==self.arcKW):
            self.imagetype='arc'
        elif(self.hdr[self.imtype]=='HARTMANN'):
            self.imagetype='hartmann'

        self.ui.lineEditExpType.setText(self.imagetype)

        # if arc, get lamp name from header
        if(self.imagetype=='arc'): 
            self.lamp = self.hdr['ARC-LAMP']
            # unset WLC flag
            self.wlc=False
            self.wavelen=np.arange(self.data2D.shape[1])
            self.ui.lineEditArcLamp.setText(self.lamp)
        else:
            self.ui.lineEditArcLamp.setText('')


        # get grating name (eventually from header, for now from combo box)
        #self.grating=self.ui.comboBoxGrating.currentText()
        self.grating = self.hdr['grating']
        self.grangle = self.hdr['GR-ANGLE']
        self.ui.lineEditGrating.setText(self.grating)
        self.ui.lineEditGratingAngle.setText(self.grangle)
        print "[%s] [%s][%s] [%ss]\n[pipeline recipe: %s]" % (self.hdr[self.imtype], self.grating, self.grangle, self.hdr['exptime'], self.imagetype)
#        if(self.imagetype=='arc'): print self.lamp

    """
    def loadAndUpdate(self):
        self.loadFrame()
        ###self.updateAll()
    """

    def loadFrame(self):
        # first see if file is valid
        frame = self.ui.lineEditCurrFrame.text()
        filename = '%s%s'%(self.pars['dataDir'],frame)
        self.logger.debug(filename)
        try:
            data2D,hdr = fits.getdata(filename, header=True)
            self.logger.debug('success!')
        except:
            try:
                data2D,hdr = fits.getdata(filename+'.fits', header=True)
                frame = frame+'.fits'
            except:
                self.logger.info('%s could not be loaded'%filename)
                return
    
        # flip so blue is on left, as God intended!
        data2D = data2D[:,::-1]
        self.hdr = hdr
        self.pars['currFrame']=frame
        #*** refresh all
        imname = frame
        self.logger.info(imname)
        
        self.image = self.pars['dataDir']+imname
        #data2D,hdr = fits.getdata(self.image, header=True)
        self.data2D = data2D
        self.sy,self.sx = np.shape(data2D)
        self.pars['currFrame'] = imname
        self.pars['obsType'] = hdr[self.imtype]
        print self.pars['obsType']
        self.determineObsType()
        
        
        if self.imagetype=='science' : self.doScience()
        if self.imagetype=='arc' : self.doArc()
        if self.imagetype=='hartmann': self.doHartmann()
        if self.imagetype=='other': self.doOther()

        """
        self.plotProfile()
        self.updatezScaling()
        self.display2D()
        self.display1D()
        """
    """
    def loadLatestFrame(self):
        try:
            tmpdata2D,hdr = pyfits.getdata(self.latestFrame,header=True)
        except:
#            time.sleep(1) # seems too quick and file is not finished writing??***
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
    """
    
    
    
    def updatezScaling(self):
        if(str(self.ui.comboBoxScaling.currentText())=='ZScale'): 
            self.z0,self.z1=zscale(self.data2D)
        if(str(self.ui.comboBoxScaling.currentText())=='Min/Max'): 
            self.z0=np.min(self.data2D)
            self.z1=np.max(self.data2D)
        if(str(self.ui.comboBoxScaling.currentText())=='Manual'): 
            self.z0=float(self.ui.lineEditz0.text())
            self.z1=float(self.ui.lineEditz1.text())
        
        self.logger.debug("updating scaling")
        self.ui.lineEditz0.setText("%.3f"%(self.z0))
        self.ui.lineEditz1.setText("%.3f"%(self.z1))
#        self.displayRaw2D()
        try: 
            self.display2D()
            self.plotProfile()
        except:
            pass
        
 
        
        
        
        
    def updateAll(self):
#        self.findLatestFrame()
#        self.loadLatestFrame()

#        self.determinezScaling()
        self.updatezScaling()

        self.display2D()
        #print self.imagetype
        self.logger.debug("updating All")

        if(self.imagetype=='science'):   
#            self.profile1D()
            self.fitSky()
            self.extract1D()
            self.plotProfile()
            self.updateExtractionSliders()
        
#            self.doSkySubtraction()
            self.display2D()
            self.display1D()
        else:
            self.data2DSkySub=self.data2D
#            self.data2DSkySub=np.copy(self.data2D)*0.0
            #self.plotProfile()
            self.medsky1d=0.0
            self.medskylev=0.0
            self.display2D()
            self.display1D()
#            self.skysub2DView.canvas.clear()
        
        self.extract1D()
        self.display1D()
        self.checkOutputReq()
        self.writeResults()
        
           
    def checkOutputReq(self):    
       if (self.ui.checkBoxWrite1Dascii.isChecked):
           self.writeAscii=True
       else:
           self.writeAscii=False
       if (self.ui.checkBoxWrite1DFITS.isChecked):
           self.writeFits=True
       else:
           self.writeFits=False
           
           
       
    def writeResults(self):
        
        if(self.writeAscii):
            outfile1D = str( self.pars['outDataDir']+self.pars['currFrame']).replace('.fits','_1D.txt') 
            #            outfile1D = 'text.txt'
#            print outfile1D
            x = np.arange(len(self.obj1d))
            if(self.wlc):
                xxx=self.wavelen
            else:
                xxx = np.arange(len(self.obj1d))
#            np.savetxt(outfile1D, np.transpose( (xxx,self.obj1d,x) ) )
            if(self.imagetype=='science'):
                np.savetxt(outfile1D, np.transpose( (xxx,self.obj1d,self.sky1d,x) ), header='wavelen objCounts skyAndBiasCounts pixnum' )
            else:
                np.savetxt(outfile1D, np.transpose( (xxx,self.obj1d,x*0.,x) ), header='wavelen objCounts skyAndBiasCounts pixnum' )
            self.logger.debug('written %s'%outfile1D)
#            self.logger.info('written %s'%outfile1D)
        if(self.writeFits):
            pass
        
        
    def updateAutoExtraction(self):
        self.autoExtraction=True
        self.updateAll()
        

        
    def showDialogSetDatadir(self):
#        print "foo"
        #fname = QtGui.QFileDialog.getOpenFileName(None,'Open file', self.parentDataDir)
        fname = str(QtGui.QFileDialog.getExistingDirectory(None, "Set directory for input data"))
        self.pars['dataDir']=fname+'/' 
        self.notifier.stop()
        wm = WatchManager()
        notifier = ThreadedNotifier(wm, self.process_IN_CREATE)
        self.notifier = notifier

        print "input data directory is now: "+str(self.pars['dataDir'])
        self.ui.lineEditDataDir.setText(self.pars['dataDir'])
        self.notifier.start()
        mask = EventsCodes.ALL_FLAGS['IN_CREATE']
        wdd = wm.add_watch(self.pars['dataDir'], mask, rec=False)
#        sel

    def showDialogSetOutDatadir(self):
#        print "foo"
        #fname = QtGui.QFileDialog.getOpenFileName(None,'Open file', self.parentDataDir)
        fname = str(QtGui.QFileDialog.getExistingDirectory(None, "Set directory for writing output data"))
        self.pars['outDataDir']=fname+'/' 
#        print "data directory is now: "+str(self.parentDataDir)
        print "output data directory is now: "+str(self.pars['outDataDir'])
        self.ui.lineEditOutDataDir.setText(self.pars['outDataDir'])

    def showDialogSetCurrFrame(self):
#        print "foo"
        fname = QtGui.QFileDialog.getOpenFileName(None,'Open file', self.pars['dataDir'])
        #fname = str(QtGui.QFileDialog.getExistingDirectory(None, "choose frame"))
#        self.pars['dataDir']=fname+'/' 
#        print "input data directory is now: "+str(self.pars['dataDir'])
#        self.ui.lineEditDataDir.setText(self.pars['dataDir'])
        frame = fname.split('/')[-1]
        self.ui.lineEditCurrFrame.setText(frame)
#        filename = '%s%s'%(self.pars['dataDir'],frame)
        self.pars['currFrame'] = frame

        # update GUI
#***        self.
        self.loadFrame()
        self.updateAll()


    def closeEvent(self):

        quit_msg = "Are you sure you want to quit?"
        #reply = QtGui.QMessageBox.question(self, 'Message', quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        reply = QtGui.QMessageBox.question(None, 'Removal', quit_msg, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if (reply == QtGui.QMessageBox.Yes):
#            QtGui.QApplication.quit() # Nope!
            self.notifier.stop()
            QtGui.qApp.closeAllWindows()
        
        else:
            pass

    # -- This is kinda dumb but works!:
    def monitorFileSize(self):
        # watch size of newly created file and return when it stops increasing
        lastSize=0.0
        currSize = os.path.getsize(self.latestFrame)
#        print currSize
        print "monitoring read out..."
        while( (currSize>lastSize) | (currSize==0)):
#            print currSize
            lastSize=currSize
            #time.sleep(1.0)
            time.sleep(0.1)
            currSize = os.path.getsize(self.latestFrame)
#            print currSize
            #pass
#        print self.latestFrame+" finished reading out."
        print "finished reading out."
        currSize=0.0 # reset


    def process_IN_CREATE(self, event):
        newfile=os.path.join(event.path, event.name)
        # make sure .fit* file
        ext=newfile.split('.')[-1]
        try:
            if (ext[0:3] <>'fit'):
                return
        except:
            return
        
#        print "%s found: %s" % (time.ctime(), os.path.join(event.path, event.name))
        print "%s found: %s" % (time.ctime(), newfile)

        
        
        frameName = (newfile.split('/')[-1]).split('.')[0]
        self.latestFrame = newfile
        self.monitorFileSize()
        
        self.ui.lineEditCurrFrame.setText(frameName)
        self.loadFrame()
###        self.updateAll()
        
        #print "waiting %.1f sec"%(self.ReadOutDelay)
        #time.sleep(self.ReadOutDelay) # wait for copy to finish??
        

        """        # this is replacement for findLatestFrame:
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
        """
def run():
    import sys
#    app = QtGui.QApplication(sys.argv)
#    MainWindow = QtGui.QMainWindow()
#    ui = Ui_MainWindow()
#    ui.setupUi(MainWindow)
#    MainWindow.show()
#    sys.exit(app.exec_())
    #setup_logging('_tmp.log')
    #log.info("Starting {}.".format(__file__))
    app = QtGui.QApplication([])
    myapp = QLGui()
    myapp.setWindowTitle('SpUpNIC Quick Look')
    myapp.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    run()
