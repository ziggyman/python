import time
from PyQt4 import QtCore, QtGui

#import logging as log
import logging, StringIO, time
import numpy as np

from scipy.stats import nanmedian

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


    def sideBandToggled(self):
        self.pars['useSB1'] = self.ui.useSB1Check.isChecked()
        self.pars['useSB2'] = self.ui.useSB2Check.isChecked()
#        print '<<<',self.pars['useSB1']
        self.updateAll()

    def findObject(self):
        # find the best candidate peak for object (i.e. nearest centre of slit).
        # deal with multiple peaks
        """
        print 'findObject called!'
        print self.pars['extractWin_x0']
        print self.extractWin_x0
        """
        image = self.data2D
        xsec = np.median(image,1)
        xx=np.arange(self.data2D.shape[0])
        ny = self.data2D.shape[0]

        pkmax,pkmin=peakdet(xsec,20)
        npks=np.shape(pkmax)[0]
        #print npks
        self.logger.debug('npks %s'%npks)
        #*** only works for 1 peak, so far
        if npks > 1:
            #print 'multiple peaks'
            cendis = 9999.
            cenx = ny/2
            for i in range(npks):
                tdis = np.abs( np.array(pkmax[i][0])-cenx )
                if tdis < cenx:
                    besti = i
                    cendis = tdis
        else:
            besti = 0

        self.logger.debug(pkmax)
        try:
            pkpos = np.array(pkmax)[besti][0]
        except:
            pkpos,pkmax,npks = ny/2,-1,0
        #print ':::::',npks,pkmax,pkpos
        self.logger.debug(':::::',npks,pkmax,pkpos)

        self.whereProfMax = pkpos
        self.profMax = xsec[pkpos]

        # try to fit gaussian to peak:
        # fit Gaussian to position
        width=5.0 # initial guess at object width

        if(0): # Gaussian fitting
            x=np.arange(len(xsec))
            y=xsec
            e=np.ones(len(xsec))
            popt, pcov = curve_fit(gaussian, x, y, sigma=e)
            self.logger.debug(popt)
            fitpos=popt[1]
            fitsigma=np.min([3.0*np.abs(popt[2]),2.0]) # ** careful, sigma can be -ve.!

            self.logger.debug(fitpos,fitsigma)
            profile1d=gaussian(xx,popt[0],popt[1],popt[2],popt[3])
            #sc.canvas.ax.set_title('Gauss. peak: y = %.0f'%fitpos)
            # but need better error trapping!
            print '>>>',popt


        # if bad: return
        popt = np.array([100.,pkpos,width])


        return popt




    def plotProfile(self):
        #sc = matplotlibWidgetwBar(self.ui.centralWidget)
        #self.profPlot.setGeometry(QtCore.QRect(1100, 140, 330, 271))
        sc = self.profWin
        
        image = self.data2D
        xsec = np.median(image,1)

        popt = self.findObject()

        """
        self.profMax = np.max(xsec)
        self.whereProfMax = np.argmax(xsec)
        """

#        self.whereProfMax = self.sy/2
#        self.profMax = xsec[self.whereProfMax]




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
            """
            pkmax,pkmin=peakdet(xsec,20)
            npks=np.shape(pkmax)[0]
            self.logger.debug('npks %s'%npks)
            #*** only works for 1 peak, so far

            self.logger.debug(pkmax)
            pkpos = np.array(pkmax[0][0])
            print ':::::',npks,pkmax
            sc.canvas.ax.axhline(pkpos, color='b',alpha=0.2)
            sc.canvas.ax.set_title('peak: y = %.0f'%pkpos)
            sc.canvas.draw()
            ##
            sc.show()
            ##
            """

#            pkpos = self.findObject()
            popt = self.findObject()

            """
            # try to fit gaussian to peak:
                # fit Gaussian to position
            width=5.0 # initial guess at object width


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
 
            """
            e=np.ones(len(xsec))
            popt, pcov = curve_fit(gaussian, x, y, sigma=e)
            self.logger.debug(popt)
            fitpos=popt[1]
            fitsigma=np.min([3.0*np.abs(popt[2]),2.0]) # ** careful, sigma can be -ve.!

            self.logger.debug(fitpos,fitsigma)
            profile1d=gaussian(xx,popt[0],popt[1],popt[2],popt[3])
            sc.canvas.ax.set_title('Gauss. peak: y = %.0f'%fitpos)
            # but need better error trapping!
            print '>>>',popt
            """

##            sc.canvas.ax.plot(profile1d,xx,'g-')
            
            #-- calculate extract,sideband windows. Need some error checking!
            extract_halfwidth=3.0*np.abs(popt[2]) # fitted width can be -ve(!) which really screws things up!
            extract_halfwidth=np.min([extract_halfwidth,8.0])
            # Need proper error-checking on these extraction pars! ***
            SBoffset = self.pars['SBoffset'] #5.0
            SBwidth = self.pars['SBwidth'] #15.0

#            print self.pars['useSB1']

            if (self.autoExtraction):
#                print '**',popt[1],extract_halfwidth
                self.extractWin_x0=popt[1]-extract_halfwidth
                self.extractWin_x1=popt[1]+extract_halfwidth
                self.SB1_x1=self.extractWin_x0-SBoffset
                self.SB1_x0=self.SB1_x1-SBwidth
                self.SB2_x0=self.extractWin_x1+SBoffset
                self.SB2_x1=self.SB2_x0+SBwidth

#                print self.SB1_x0,self.SB1_x1

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
                sc.canvas.ax.fill_between(np.arange(0,self.profMax),self.SB1_x0,self.SB1_x1,color='r',alpha=0.3)


            if(self.pars['useSB2']):
                sc.canvas.ax.axhline(self.SB2_x0,color='r',alpha=0.3)
                sc.canvas.ax.axhline(self.SB2_x1,color='r',alpha=0.3)
                sc.canvas.ax.fill_between(np.arange(0,self.profMax),self.SB2_x0,self.SB2_x1,color='r',alpha=0.3)
            #sc.canvas.ax.axis('tight')
#            sc.canvas.ax.axis('tight')

            
            
        if(1):
            sc.canvas.ax.axhline(self.extractWin_x0,color='g',alpha=0.5)
            sc.canvas.ax.axhline(self.extractWin_x1,color='g',alpha=0.5)
#            sc.canvas.ax.axis('tight')
            sc.canvas.ax.fill_between(np.arange(0,self.profMax),self.extractWin_x0,self.extractWin_x1,color='g',alpha=0.3)
            
        # suppress printing of x,y coords:
        sc.canvas.ax.format_coord = lambda x, y: ''

        sc.canvas.draw()
        ##
        sc.show()
        ##

        
    def fitSky(self):
#        if not self.ui.checkBoxdoSkySub.isChecked():
#            self.sky1d = np.zeros(self.data2D.shape[1])
#            return
        #self.useSB1 = True
        #self.useSB2 = True
        #print self.imagetype
    
        b0 = np.floor(self.SB1_x0) ; b1 = np.ceil(self.SB1_x1)
        b2 = np.floor(self.SB2_x0) ; b3 = np.ceil(self.SB2_x1)
    
        sky01 = np.copy(self.data2D[b0:b1,:])
        sky23 = np.copy(self.data2D[b2:b3,:])

        if (self.pars['useSB1'] == False):
            sky01[:,:] = np.nan#sky01*0.0
        if (self.pars['useSB2'] == False):
            sky23[:,:] = np.nan#sky23*0.0
        
        sky03 = np.append(sky01,sky23,axis=0)        
        
#        print 'here'
#        print np.shape(sky03)
#        sky1d = np.median(sky03,0)
#        notNaN = np.where(np.isfinite(sky03))
#        sky1d = np.median(sky03[notNaN[0],notNaN[1]],0)
        sky1d = nanmedian(sky03,0)

        if np.sum(np.isnan(sky1d))==np.count_nonzero(sky1d) : # in case both sky bands are turned off, replace with zeros
            #print 'sky all nans'
            sky1d[:]=0.

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
        self.pars['extrType'] = str(extrType)

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
                # for summed counts, we need to multiply the sky [per pixel] by the number of pixels in the extraction window:
                self.sky1d = self.sky1d * float(b1-b0)
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

        