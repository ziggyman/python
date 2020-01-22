import time
from PyQt4 import QtCore, QtGui

#import logging as log
import logging, StringIO, time

from cassspectr.controller.controller import Controller

plc_host = 'localhost'
plc_port = 9090
detector_host = 'localhost'
detector_port = 9091

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class logBuffer(QtCore.QObject, StringIO.StringIO):
    bufferMessage = QtCore.pyqtSignal(str)

    def __init__(self, *args, **kwargs):
        QtCore.QObject.__init__(self)
        StringIO.StringIO.__init__(self, *args, **kwargs)

    def write(self, message):
        if message:
            self.bufferMessage.emit(unicode(message))

        StringIO.StringIO.write(self, message)

class logPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl





