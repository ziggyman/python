
import time
from PyQt4 import QtCore, QtGui

#import logging as log
import logging, StringIO, time
import numpy as np
try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from scipy.optimize import curve_fit

from astropy.modeling import models, fitting

from specutils import peakdet, gaussian

from astropy.io import fits
from zscale import zscale

class ProfPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl
        
    def plotProfile(self):
        #sc = matplotlibWidgetwBar(self.ui.centralWidget)
        #self.profPlot.setGeometry(QtCore.QRect(1100, 140, 330, 271))
        sc = self.profWin
        
        image = self.data2D
        xsec = np.median(image,1)
        
        self.profMax = np.max(xsec)
        self.whereProfMax = np.argmax(xsec)

        geom = self.ui.profileGraph.geometry()
        sc.setGeometry(geom)

        sc.canvas.ax.clear()
        xx = np.arange(len(xsec))
        #sc.canvas.ax.set_aspect("auto")
        #sc.canvas.ax.plot(xx,xsec,'b-')
        sc.canvas.ax.plot(xsec,xx,'k-',lw=2)
        sc.canvas.ax.set_ylim((xx.min(),xx.max()))
        sc.canvas.draw()
        ##
        sc.show()
        ##
        
        self.ax0 = sc.canvas.ax.axis()

        """
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
        """
#        self.ui.profWidget.canvas.draw()
        
#        def getPeak(self):
        if(self.imagetype=='science'): # only do on science frame
            pkmax,pkmin=peakdet(xsec,20)
            npks=np.shape(pkmax)[0]
            self.logger.debug('npks %s'%npks)
            #*** only works for 1 peak, so far

            self.logger.debug(pkmax)
            pkpos = np.array(pkmax[0][0])
            sc.canvas.ax.axhline(pkpos, color='b',alpha=0.2)
            sc.canvas.ax.set_title('peak: y = %.0f'%pkpos)
            sc.canvas.draw()
            ##
            sc.show()
            ##

            # try to fit gaussian to peak:
                # fit Gaussian to position
            width=5.0 # initial guess at object width
            
            """    
            p0=[pkmax[0][1],pkmax[0][0],width,np.min(xsec)]

#            g_init = models.Gaussian1D(amplitude=pkmax[0][1], mean=pkmax[0][0], stddev=width, n_models=1)
            g_init = models.Gaussian1D(amplitude=[pkmax[0][1],200], mean=[pkmax[0][0],len(xx)/2.], stddev=[width,100], n_models=2)
            #g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
            fit_g = fitting.LevMarLSQFitter()
            x=np.arange(len(xsec))
            y=xsec
            g = fit_g(g_init, x, y)
            print g
            
            sc.canvas.ax.plot(g(xx),xx,'g-')
            """
            
            x=np.arange(len(xsec))
            y=xsec
 
                        
            e=np.ones(len(xsec))
            popt, pcov = curve_fit(gaussian, x, y, sigma=e)
            self.logger.debug(popt)
            fitpos=popt[1]
            fitsigma=np.min([3.0*np.abs(popt[2]),2.0]) # ** careful, sigma can be -ve.!
        
            self.logger.debug(fitpos,fitsigma)
            profile1d=gaussian(xx,popt[0],popt[1],popt[2],popt[3])
            
            sc.canvas.ax.plot(profile1d,xx,'g-')
            
            #-- calculate extract,sideband windows. Need some error checking!
            extract_halfwidth=3.0*popt[2]
            extract_halfwidth=np.min([extract_halfwidth,8.0])
            # Need proper error-checking on these extraction pars! ***
            SBoffset = self.pars['SBoffset'] #5.0
            SBwidth = self.pars['SBwidth'] #15.0

            if (self.autoExtraction):
                self.extractWin_x0=popt[1]-extract_halfwidth
                self.extractWin_x1=popt[1]+extract_halfwidth
                self.SB1_x1=self.extractWin_x0-SBoffset
                self.SB1_x0=self.SB1_x1-SBwidth
                self.SB2_x0=self.extractWin_x1+SBoffset
                self.SB2_x1=self.SB2_x0+SBwidth
        
                self.pars['extractWin_x0'] = self.extractWin_x0
                self.pars['extractWin_x1'] = self.extractWin_x1
                self.pars['SB1_x0'] =self.SB1_x0
                self.pars['SB1_x1'] =self.SB1_x1
                self.pars['SB2_x0'] =self.SB1_x0
                self.pars['SB2_x1'] =self.SB1_x1
                
        
                if(self.dbglvl>1):
                    self.logger.debug(self.SB1_x1,self.SB1_x0,self.SB2_x0,self.SB2_x1)
        
        self.showBands1D()
        """
            self.fitSky()
            #print self.fitSky()
            self.extract1D()
            
            # -- update 2D:
            #self.displayRaw2D()
        else:
            # just extract object (full width of frame)?
            self.extract1D()
            #**
#            self.sky1d = np.zeros(self.data2D.shape[1])    
            #**
        """
        
    def showBands1D(self):
        sc = self.profWin
        
#        if(self.pars['isScience']):
        if(self.imagetype=='science'):
            if(self.pars['useSB1']):
                sc.canvas.ax.axhline(self.SB1_x0,color='r',alpha=0.3)
                sc.canvas.ax.axhline(self.SB1_x1,color='r',alpha=0.3)
            if(self.pars['useSB2']):
                sc.canvas.ax.axhline(self.SB2_x0,color='r',alpha=0.3)
                sc.canvas.ax.axhline(self.SB2_x1,color='r',alpha=0.3)
            #sc.canvas.ax.axis('tight')
            sc.canvas.ax.fill_between(np.arange(0,self.profMax),self.SB1_x0,self.SB1_x1,color='r',alpha=0.3)
#            sc.canvas.ax.axis('tight')
            sc.canvas.ax.fill_between(np.arange(0,self.profMax),self.SB2_x0,self.SB2_x1,color='r',alpha=0.3)
            
            
        if(1):
            sc.canvas.ax.axhline(self.extractWin_x0,color='g',alpha=0.5)
            sc.canvas.ax.axhline(self.extractWin_x1,color='g',alpha=0.5)
#            sc.canvas.ax.axis('tight')
            sc.canvas.ax.fill_between(np.arange(0,self.profMax),self.extractWin_x0,self.extractWin_x1,color='g',alpha=0.3)
            
        
        sc.canvas.draw()
        ##
        sc.show()
        ##

        
    def fitSky(self):
        if not self.ui.checkBoxdoSkySub.isChecked():
            self.sky1d = np.zeros(self.data2D.shape[1])
            return
        #self.useSB1 = True
        #self.useSB2 = True
        #print self.imagetype
    
        b0 = np.floor(self.SB1_x0) ; b1 = np.ceil(self.SB1_x1)
        b2 = np.floor(self.SB2_x0) ; b3 = np.ceil(self.SB2_x1)
    
        sky01 = self.data2D[b0:b1,:]
        sky23 = self.data2D[b2:b3,:]

        if (self.pars['useSB1'] == False):
            sky01[:,:] = np.nan#sky01*0.0
        if (self.pars['useSB2'] == False):
            sky23[:,:] = np.nan#sky23*0.0
        
        sky03 = np.append(sky01,sky23,axis=0)        
        
#        print 'here'
        sky1d = np.median(sky03,0)
        
        self.sky1d = sky1d 
#        print self.sky1d 

        """
        if( self.useSB1 & self.useSB2):
            sky03=np.append(sky01,sky23,axis=0)        
        elif( self.useSB1==False ):
            sky03 = sky23
        elif( self.useSB2==False ):
            sky03 = sky01
        """ 
        
        
    def extract1D(self):
        
        # Extraction type:
        extrType = self.ui.comboBoxExtractionType.currentIndex()

        b0 = np.floor(self.extractWin_x0) ; b1 = np.ceil(self.extractWin_x1)
        #print b0,b1,'**'
        self.logger.debug("extraction window: %.2f--%.2f"%(b0,b1))
        if(self.imagetype=='science'):

            if extrType ==0:
                # Peak row:
                # just peak:
                self.obj1d = self.data2D[self.whereProfMax,:] - self.sky1d
                self.extractionUnits = 'peak row in extraction window'
            elif extrType ==1:
                self.obj1d = np.sum(self.data2D[b0:b1,:],0) - self.sky1d
                self.extractionUnits = 'summed counts within extraction window'
            elif extrType ==2:
                self.obj1d = np.median(self.data2D[b0:b1,:],0)*float(b1-b0) - self.sky1d
                self.extractionUnits = 'median counts within extraction window x window width'
            elif extrType ==3:
                self.obj1d = np.median(self.data2D[b0:b1,:],0) - self.sky1d
                self.extractionUnits = 'median counts per pixel within extraction window'
        else:
            self.obj1d = np.median(self.data2D[b0:b1,:],0) #- self.sky1d
        
        
        
        
    def subZero(self):
        # subtract zero level, either from overscan or bias, if available:
        # 0-21, 103-end
        biassec1=np.arange(0,20)+1
        biaslev = np.median(self.data2D[:,biassec1])
        print biaslev
        self.data2D = self.data2D - biaslev
        
    def readExtractionSliders(self): 
#        if (self.useSB1Check.isChecked()):
        self.SB1_x0=self.ui.SB1Slider.lowerValue
        self.SB1_x1=self.ui.SB1Slider.upperValue
#        else:
#            self.SB1_x0=None
#            self.SB1_x1=None

        """
        if(self.ui.useSB1Check.isChecked()):
            self.pars['useSB1']=True
            self.ui.useSB1Check.setEnabled(False)
        else:
            self.pars['useSB1']=False
            self.ui.useSB1Check.setEnabled(False)
        
        if(self.ui.useSB2Check.isChecked()):
            self.pars['useSB2']=True
            self.ui.useSB2Check.setEnabled(False)
        else:
            self.pars['useSB2']=False
            self.ui.useSB1Check.setEnabled(False)
        """
        
            
        self.SB2_x0=self.ui.SB2Slider.lowerValue
        self.SB2_x1=self.ui.SB2Slider.upperValue
        
        self.extractWin_x0=self.ui.extractWinSlider.lowerValue
        self.extractWin_x1=self.ui.extractWinSlider.upperValue
        self.autoExtraction=False
        self.plotProfile()
        if(self.imagetype=='science'):
            self.fitSky()
        self.extract1D()        
        self.display2D()
        self.updateExtractionSliders()
        self.display1D()
        

    def updateExtractionSliders(self):
#        self.ui.extractWinSlider.setRange(0, self.sy)
        self.ui.extractWinSlider.setSpan(self.extractWin_x0,self.extractWin_x1)
#        self.ui.SB1Slider.setRange(0, self.sy)

        # reset range so it can't overlap with extraction window:
        self.extractWin_x0=self.ui.extractWinSlider.lowerValue
        self.extractWin_x1=self.ui.extractWinSlider.upperValue
        self.ui.SB1Slider.setRange(0,self.extractWin_x0)
        if(self.SB1_x1>self.extractWin_x0):
            self.SB1_x1=self.extractWin_x0
        if(self.SB2_x0<self.extractWin_x1):
            self.SB2_x0=self.extractWin_x1
            
        self.ui.SB2Slider.setRange(self.extractWin_x1,self.sy)
        self.ui.SB1Slider.setSpan(self.SB1_x0,self.SB1_x1)
#        self.ui.SB2Slider.setRange(0, self.sy)
        self.ui.SB2Slider.setSpan(self.SB2_x0,self.SB2_x1)

        self.display2D()
        #self.doSkySubtraction()
        #self.displaySS2D()
        if(self.imagetype=='science'):
            self.fitSky()

        self.extract1D()
        self.display1D()
        self.display2D()
        #self.updatezScaling()
        #self.updateAll()

        