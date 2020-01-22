import time
from PyQt4 import QtCore, QtGui

#import logging as log
import logging, StringIO, time
import numpy as np

from scipy.stats import nanmedian

from astropy.io import fits
from glob import glob

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

class ObsLogPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl



    def mklog(self):
        files = glob('%s/*fits'%self.datadir)
        files = np.sort(files)
        self.ui.obsLogScroll.setReadOnly(True)
        self.ui.obsLogScroll.setText(files)