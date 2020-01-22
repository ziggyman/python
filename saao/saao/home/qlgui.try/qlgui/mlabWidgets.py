#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


dbglvl=2


# ----
# http://matplotlib.1069221.n5.nabble.com/key-press-events-in-matplotlib-embedded-in-pyqt4-td27958.html
class MplCanvas(FigureCanvas):
 
    def __init__(self):
#        self.fig = Figure((self.width(),self.height()))
        self.fig = Figure()
#        self.fig = Figure((1,0.3),dpi=100)
#        self.fig = Figure((100,30))
        #print self.size

        self.ax = self.fig.add_subplot(111)
 
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
 
class matplotlibWidget(QtGui.QWidget):
 
    def __init__(self, parent = None):
        QtGui.QWidget.__init__(self, parent)


        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()
        
#        print self.width(), self.height()
#        self.canvas.setMinimumSize(self.canvas.size())
        
        self.vbl.addWidget(self.canvas)
        #--
        self.mpl_toolbar = NavigationToolbar(self.canvas, self) 
        self.mpl_toolbar.hide()
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




# http://stackoverflow.com/questions/8356336/how-to-capture-output-of-pythons-interpreter-and-show-in-a-text-widget
class EmittingStream(QtCore.QObject):

    textWritten = QtCore.pyqtSignal(str)

    def write(self, text):
        self.textWritten.emit(str(text))





