
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

from specutils import peakdet

from astropy.io import fits
#from astropy.convolution import Gaussian1DKernel as gkern
#from astropy.convolution import convolve as convl
from zscale import zscale


class Im1DPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl
        #self.sky1d = self.ui.sky1d
        
        
###    def display1D(self,resetView=False):
    def display1D(self,resetView=True):
#        print 'display1d()'
        specWin = self.specWin
        geom = self.ui.spec1DView.geometry()
        specWin.setGeometry(geom)

        if np.nanmax(self.data2D)==0.:
            resetView = True # no real data loaded yet

        #** get and set previous limits
        xmin, xmax = specWin.canvas.ax.get_xlim()
        # let's just keep it x range for now in case user wants to drastically change exptime.
        #xmin, xmax = specWin.canvas.ax.get_xlim()

        specWin.canvas.ax.clear()
###        specWin.canvas.figure.add_subplot(111,sharex=self.ax2)

        self.pars['smoothPix'] = self.ui.spinBoxSmoothWidth.value()
        self.logger.debug(self.pars['smoothPix'])
        if(self.wlc):
            xxx=self.wavelen
            specWin.canvas.ax.set_xlabel(r'$\lambda(\AA)$')
        else:
            specWin.canvas.ax.set_xlabel('pix')
            xxx = np.arange(len(self.obj1d))
#        specWin.canvas.ax.set_ylabel('median counts/pix')
        specWin.canvas.ax.set_ylabel(self.extractionUnits)
        if(self.pars['smoothPix']<=1):
            specWin.canvas.ax.plot(xxx,self.obj1d,'k-',label='object')
#            specWin.canvas.ax.plot(xxx,self.obj1d,'k-',label='object')
            #print self.obj1d
            #print 'here'
        else:
#            specWin.canvas.ax.clear()
#            kern = gkern(self.pars['smoothPix'])
#            smobj1d = convl(self.obj1d,kern)

            nn = self.pars['smoothPix']
            kern=np.ones(nn)/float(nn)
            smobj1d=np.convolve(self.obj1d,kern,mode='same')
            specWin.canvas.ax.plot(xxx,smobj1d,'k-',label='object')
        if(self.imagetype=='science'):
            specWin.canvas.ax.plot(xxx,self.sky1d,'r-',label='sky\n(+bias\n   level)')
        
        
        
##        specWin.canvas.ax.set_xlim([xxx.min()*0.99,xxx.max()*1.01])
        
#        specWin.canvas.ax.subplots_adjust(bottom=0.2)
        #specWin.canvas.ax.set_ymargin=0.2
        specWin.canvas.ax.legend(bbox_to_anchor=(1.12, 1.0))
        specWin.canvas.draw()
        ##
        specWin.show()
        ##

        # suppress printing of x,y coords:
        specWin.canvas.ax.format_coord = lambda x, y: ''

        #** set to previous zoom
        if not resetView:
            specWin.canvas.ax.set_xlim([xmin,xmax])
        specWin.canvas.draw()
        ##
        specWin.show()


#        self.ax1 = specWin.canvas.ax#.axis()
###        specWin.canvas.figure.add_subplot(111,sharex=self.ax2)
        
        

    def write1Dascii(self):
        if(self.wlc):
            xxx=self.wavelen
        else:
            xxx = np.arange(len(self.obj1d))
        outfile = '%s%s_1D.txt'%(self.pars['outDataDir'],self.pars['frame'])
        self.logger.info('written %s'%outfile)
        np.savetxt(xxx,self.obj1d)

    
    def display1DHart(self):
#        print 'display1d()'
        specWin = self.specWin
        geom = self.ui.spec1DView.geometry()
        specWin.setGeometry(geom)

        specWin.canvas.ax.clear()

        if(0):
            xxx=self.wavelen
        else:
            specWin.canvas.ax.set_xlabel('pix')
            xxx = np.arange(len(self.obj1d))
        specWin.canvas.ax.set_ylabel('median counts/pix')
        if(1):
            #specWin.canvas.ax.plot(xxx,self.obj1d,'k-',label='object')
            specWin.canvas.ax.plot(xxx,self.hartAspec1d,'r-',label='A')
            specWin.canvas.ax.plot(xxx,self.hartBspec1d,'b-',label='B')
            specWin.canvas.ax.plot(xxx-self.thisShiftA,self.hartAspec1d,'r--',label='best-fit shift A')
            #print self.obj1d
            #print 'here'



##        specWin.canvas.ax.set_xlim([xxx.min()*0.99,xxx.max()*1.01])

#        specWin.canvas.ax.subplots_adjust(bottom=0.2)
        #specWin.canvas.ax.set_ymargin=0.2
        specWin.canvas.ax.legend(bbox_to_anchor=(1.12, 1.0))
        specWin.canvas.draw()
        ##
        specWin.show()
        ##

        self.ax1 = specWin.canvas.ax.axis()

