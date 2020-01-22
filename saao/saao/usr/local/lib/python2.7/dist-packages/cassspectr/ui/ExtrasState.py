"""Store state of everything (in dictionary) which is not already dealt with by PLC/Detector"""

# this should probably be made part of Controller later***

import logging as log
import numpy as np

class ExtrasState:
    def __init__(self):
        
        headInfo = {}
        # header Pane:
        headInfo['OBSERVER'] = ''#(self.ui.lineEditObserver.text(), 'OBSERVER NAME')
        headInfo['PROP-ID'] = ''#(self.ui.lineEditPropID.text(), 'PROPOSAL ID')
        headInfo['GRATING'] = ''#(self.ui.comboBoxGrating.currentText(), 'GRATING NAME')
        headInfo['COMMENT'] = ''#(self.ui.lineEditComment.text(), '')
        # Instrument Pane:
        headInfo['GR-ANGLE'] = 0.0#(self.ui.labGratingAngleCurrent.text(), 'GRATING ANGLE')
#        headInfo['SLITWID'] = ("%.4f"%float(self.ui.labSlitWidthCurrent.text()) ,'SLIT WIDTH (arcsec)')
        headInfo['SLITWID'] = 0.0#(self.ui.labSlitWidthCurrent.text() ,'SLIT WIDTH (arcsec)')
        headInfo['FILTER'] = ''#(self.ui.comboBoxFilter.currentText(), 'FILTER NAME')
        # CCD Pane:
        headInfo['FRAME'] = ''#(self.ui.displFrameNum.text(), 'FRAME NUMBER')
        #actualExptime = self.ui.displExpTime.text()
        #print '*** Need to take actual time exposed, not requested here ***' # *** TODO
        headInfo['EXPTIME'] = 0.0#(actualExptime, 'EXPOSURE TIME (s)')
        headInfo['EXPTYPE'] = ''#(self.ui.comboExpType.currentText(), 'EXPOSURE TYPE')
        headInfo['BINNING'] = ''#(self.ui.comboCCDbinning.currentText(), 'CCD binning')
#        headInfo['CCDSUM']
        # don't call GAIN or RDNOISE, to avoid possible confusion with measured values later
        headInfo['SGAIN'] = ''#(self.ui.comboCCDmodeGain.currentText(), 'CCD Mode')
        headInfo['SRDNOISE'] = ''#(self.ui.comboCCDmodeNoise.currentText(), 'CCD Mode')

        headInfo['CCDTEMP'] = 0.0#(self.ui.labCCDTemp.text(), 'CCD temperature')
        headInfo['CF-TEMP'] = 0.0#(self.ui.labColdPlateTemp.text(), 'Cold Plate Temperature')
        # Target Pane:
#        headInfo['TARGNAME'] = 
#        headInfo['TARG-RA'] = 
#        headInfo['TARG-DEC'] = 
#        headInfo['EPOCH'] = 
        # ZZZ
