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

from im2DPane import Im2DPane


dbglvl=2


from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

# http://stackoverflow.com/questions/12695678/how-to-modify-the-navigation-toolbar-easily-in-a-matplotlib-figure-window
class PanOnlyToolbar(NavigationToolbar):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar.toolitems if
#    t[0] in ("Home", )]
                 t[0] in ("Pan", "Zoom", "Home")]

    def __init__(self, *args, **kwargs):
        super(PanOnlyToolbar, self).__init__(*args, **kwargs)
#        self.layout().takeAt(1)  #or more than 1 if you have more buttons
        self.layout().takeAt(3)  #or more than 1 if you have more buttons
#        self.layout().takeAt(0)  #or more than 1 if you have more buttons




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


    """
    # http://stackoverflow.com/questions/30237175/matplotlib-gui-imshow-coordinates
    def mouse_moved(self, mouse_event):
        if mouse_event.inaxes:
            x, y = mouse_event.xdata, mouse_event.ydata
            self.emit(QtCore.SIGNAL("MOUSE_MOVED"),x,y)
#            self.xcoo,self.ycoo = x,y

#            print x,y

#            mouse_moved_signal.emit(x,y)
    """


class matplotlibWidget(QtGui.QWidget):
 
    def __init__(self, parent = None):
#    def __init__(self, parent =):
        QtGui.QWidget.__init__(self, parent)


        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()


###        self.canvas.mpl_connect('motion_notify_event', self.canvas.mouse_moved)

#        self.connect(Im2DPane.display2D.graph2D,QtCore.SIGNAL("MOUSE_MOVED"),Im2DPane.display2Dlabel)

#        print self.width(), self.height()
#        self.canvas.setMinimumSize(self.canvas.size())
        
        self.vbl.addWidget(self.canvas)
        #--
#        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        self.mpl_toolbar = PanOnlyToolbar(self.canvas, self)
        self.mpl_toolbar.hide()
        self.vbl.addWidget(self.mpl_toolbar)
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
#        self.mpl_toolbar = NavigationToolbar(self.canvas, None)
        self.mpl_toolbar = PanOnlyToolbar(self.canvas, self)
        self.vbl.addWidget(self.mpl_toolbar)
        #--
        self.setLayout(self.vbl)




# http://stackoverflow.com/questions/8356336/how-to-capture-output-of-pythons-interpreter-and-show-in-a-text-widget
class EmittingStream(QtCore.QObject):

    textWritten = QtCore.pyqtSignal(str)

    def write(self, text):
        self.textWritten.emit(str(text))





