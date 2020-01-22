
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
from zscale import zscale

class Im2DPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl

        
        
    def display2D(self):
        graph2D = self.graph2D
        
        geom = self.ui.raw2DGraph.geometry()
        
        graph2D.setGeometry(geom)
        graph2D.canvas.ax.clear()
        
#        z0,z1 = zscale(self.data2D)
        z0 = float(self.ui.lineEditz0.text())
        z1 = float(self.ui.lineEditz1.text())
        
        # aspect=auto is crucial to filling window!
        graph2D.canvas.ax.imshow(self.data2D,vmin=z0,vmax=z1, extent=(0., self.data2D.shape[1], 0., self.data2D.shape[0]), aspect='auto', origin='lower',\
                                 interpolation='nearest', cmap=str(self.ui.comboBoxColorMap.currentText()))
        #  get_cmap('')
        # to list options
        
        #foo = graph2D.gca()
        # remove axis labels:
        graph2D.canvas.ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
        graph2D.canvas.ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
        #graph2D.canvas.ax.format_coord = lambda x, y: ''
        graph2D.canvas.draw()
        ##
        graph2D.show()
        ##
        self.ax2 = graph2D.canvas.ax.axis()
        
        #        if (self.showSideBands):    

#        self.showSideBands=True
#        self.showExtractionWindow=True

        nx = self.data2D.shape[1]
        if (self.ui.checkBoxShowSB.isChecked() & (self.imagetype=='science')):
            if(self.pars['useSB1']):
                graph2D.canvas.ax.axhline(self.SB1_x0,color='r',alpha=0.3)
                graph2D.canvas.ax.axhline(self.SB1_x1,color='r',alpha=0.3)
            if(self.pars['useSB2']):
                graph2D.canvas.ax.axhline(self.SB2_x0,color='r',alpha=0.3)
                graph2D.canvas.ax.axhline(self.SB2_x1,color='r',alpha=0.3)
#            graph2D.canvas.ax.axis('tight')
            graph2D.canvas.ax.fill_between(np.arange(0,nx),self.SB1_x0,self.SB1_x1,color='r',alpha=0.3)
#            graph2D.canvas.ax.axis('tight')
            graph2D.canvas.ax.fill_between(np.arange(0,nx),self.SB2_x0,self.SB2_x1,color='r',alpha=0.3)

        if (self.ui.checkBoxShowExtractWin.isChecked()  & (self.imagetype=='science')):
            graph2D.canvas.ax.axhline(self.extractWin_x0,color='g',alpha=0.5)
            graph2D.canvas.ax.axhline(self.extractWin_x1,color='g',alpha=0.5)
#            graph2D.canvas.ax.axis('tight')
            graph2D.canvas.ax.fill_between(np.arange(0,nx),self.extractWin_x0,self.extractWin_x1,color='g',alpha=0.3)
            
        graph2D.canvas.draw()
        ##
        graph2D.show()
        ##
