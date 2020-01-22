
# https://wiki.python.org/moin/PyQt/Making%20non-clickable%20widgets%20clickable


import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *

def clickable(widget):

    class Filter(QObject):
    
        clicked = pyqtSignal()
        
        def eventFilter(self, obj, event):
        
            if obj == widget:
                if event.type() == QEvent.MouseButtonRelease:
                    if obj.rect().contains(event.pos()):
                        self.clicked.emit()
                        # The developer can opt for .emit(obj) to get the object within the slot.
                        return True
            
            return False
    
    filter = Filter(widget)
    widget.installEventFilter(filter)
    return filter.clicked

