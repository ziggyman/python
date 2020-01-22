import time
from PyQt4 import QtCore, QtGui

import logging as log

import numpy as np

from cassspectr.controller.controller import Controller

from astropy.io import fits as pyfits

from collections import OrderedDict

plc_host = 'localhost'
plc_port = 9090
detector_host = 'localhost'
detector_port = 9091

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


class HeaderPane:
    def __init__(self, ui):
        self.ui = ui
        self.dbglvl = self.ui.dbglvl

    #*** DESIGN CHANGE TO IMPLEMENT***:
    """Keep all variables which are not already stored in Carel's status
    dictionaries (like the PLC controller) in another python
    dictionary. This means that code such as that below is then not
    dependent on the GUI to read critical values.

    e.g. dict['GRATING NAME']=self.ui.comboBoxGrating.currentText()
    
    then the code below can search that dictionary directly..

    NOTE: most of the changes will only be needed in the header and
    target panes, but there will be some in CCD pane (such as exp
    type)

    """

    def createStartHeader(self):
        # create a brief header of values at start of exposure:
###        d = {}
        d = OrderedDict({})
        self.getTime()
        d['UT-START'] = (self.ut, 'UT time (start)')
        d['LC-START'] = (self.sast, 'local time (start)')
        d['DATE-OBS'] = (self.utdate, 'UT date (start)')
        d['HJD-OBS'] = (self.hjd, 'Heliocentric Julian Date (start)')
        # ****

        """
        DATE-OBS= '2008-05-02'                 / UT date (start)
        UT-DATE = '2008-05-02'                 / UT date (start)
        UT-TIME = '02:36:50'                   / UT time (start) 9410 seconds
        LC-TIME = '22:36:50'                   / local time (start) 81410 seconds
        AIRMASS =                1.096         / airmass (start)
        TEL-ELEV=                65.77         / telescope elevation
        ST      =              45339.8         / sidereal time: 12:35:39 (start)
        """

        return d

    # The way doneReadout() is currently set up, this header is
    # written at the END of the exposure. So, we could calculate
    # actual values here (like time of midpoint of exposure,...):
    def createHeader(self):
#    def gatherKeywords(self):
        # create a dictionary of relevant keywords:
###        d = {}
        d = OrderedDict({})
        # header Pane:
        d['OBSERVER'] = (self.settings['observer'], 'OBSERVER NAME')
        d['RUN-NO'] = (self.settings['propID'], 'RUN NUMBER')
        d['GRATING'] = (str(self.ui.labGratingName.text()), 'GRATING NAME')
        d['COMMENT'] = (str(self.ui.lineEditComment.text()), '')
        # Instrument Pane:
        d['GR-ANGLE'] = (str(self.ui.labGratingAngleCurrent.text()), 'GRATING ANGLE')
        slitString = str(self.ui.labSlitWidthCurrent.text())
        #print slitString
        slitwidArcsec = slitString.split(' ')[-1].replace('"','')
        slitwidPos = slitString.split(' ')[0].replace(':','')
        d['SLITPOS'] = (str(slitwidPos) ,'SLIT WIDTH (index)')
        d['SLITWID'] = (str(slitwidArcsec) ,'SLIT WIDTH (arcsec)')
        d['FILTER'] = (str(self.ui.comboBoxFilter.currentText()), 'FILTER NAME')
        # CCD Pane:
        d['FRAME'] = (self.settings['frameNum'], 'FRAME NUMBER')
        d['OBJECT'] = (str(self.ui.displExpName.text()), 'Target name')
        actualExptime = str(self.ui.displExpTime.text())
        print '*** Need to take actual time exposed, not requested here ***' # *** TODO
        d['EXPTIME'] = (actualExptime, 'ACTUAL EXPOSURE TIME (s)')
        d['EXPTYPE'] = (self.settings['expType'], 'EXPOSURE TYPE')
        d['BINNING'] = (self.settings['ccdBinning'], 'CCD binning')
        d['CCDSUM'] = (self.settings['ccdBinning'].replace('x',' '), 'CCD binning')
        d['SGAIN'] = (str(self.ui.comboCCDmodeGain.currentText()), 'CCD Mode GAIN')
        d['SRDNOISE'] = (str(self.ui.comboCCDmodeNoise.currentText()), 'CCD Mode RDNOISE')
        d['CCDTEMP'] = (str(self.ui.labCCDTemp.text()), 'CCD temperature')
        d['CF-TEMP'] = (str(self.ui.labColdFingerTemp.text()), 'Cold Finger Temperature')
        d['ARC-LAMP'] = (self.settings['ARC'], 'arc lamp ID')
        d['FOCUSPOS'] = (self.settings['focusPosition'], 'Focus position (mm)')
        d['HARTPOS'] = (self.settings['HartmannState'], 'Position of Hartmann Shutter')
        if self.targetsLoaded:
            d['TARG-RA'] = self.settings['TARG-RA']
            d['TARG-DEC'] = self.settings['TARG-DEC']
#            d['TARGNAME'] = self.settings['TARG-NAME']


        try:
            d['HARTSEQ'] = (self.settings['HartmannSeq'], 'Position in Hartmann sequence')
        except:
            d['HARTSEQ'] = ('0/0', 'Position in Hartmann sequence')
#            pass # may not be defined yet if first Hartmann sequence hasn't been started
        return d

        """
        other keywords to use:
        [may not be assigned here]

        arc1name
        arc2name
        arc1on
        arc2on
        grating angle
        slitwidth
        exp start time ut
        exp end time
        exp midpoint
        airmass start
        airmass end
        airmass mid

from LDSS3 (edited):

TELESCOP= 'Clay_Mag_2'                 / telescope
SITENAME= 'LCO'
SITEALT =                 2400         / meters
SITELAT =            -29.01423
SITELONG=            -70.69242
TIMEZONE=                    4
DATE-OBS= '2008-05-02'                 / UT date (start)
UT-DATE = '2008-05-02'                 / UT date (start)
UT-TIME = '02:36:50'                   / UT time (start) 9410 seconds
UT-END  = '02:36:52'                   / UT time (end) 9412 seconds
LC-TIME = '22:36:50'                   / local time (start) 81410 seconds
NIGHT   = '01May2008'                  / local night
INSTRUME= 'LDSS3-Two'                  / instrument name
SCALE   =                0.188         / arcsec/pixel
EGAIN   =                 0.75         / electrons/DU (est.)
ENOISE  =                 4.07         / electrons/read (est.)
DISPAXIS=                    1         / dispersion axis
RA      = ' 13:48:03.5'                / right ascension
RA-D    =          207.0144583         / degrees
DEC     = '-11:40:05.0'                / declination
DEC-D   =          -11.6680556         / degrees
EQUINOX =           2000.00000         / of above coordinates
ASECS   =                  0.0         / arcseconds
DSECS   =                  0.0         / arcseconds
EPOCH   =           2008.33356         / epoch (start)
AIRMASS =                1.096         / airmass (start)
TEL-ELEV=                65.77         / telescope elevation
ST      =              45339.8         / sidereal time: 12:35:39 (start)
ROTANGLE=               207.85         / rotator offset angle
ROTATORE=               -13.89         / rotator encoder angle
FILENAME= 'ccd4048c2'
EXPTYPE = 'Object'                     / exposure type
EXPTIME =                2.000         / exposure time
NLOOPS  =                    1         / # of loops per sequence
LOOP    =                    1         / # within this sequence
BINNING = '1x1'                        / binning
SPEED   = 'Slow'                       / readout speed
NOVERSCN=                    0 / overscan
NBIASLNS=                    0 / biaslines
BIASSEC = '[2033:2160,750:3500]'       / NOAO: bias section
APERTURE= '1347_11a'                   / slit-mask
FILTER  = 'VR'                         / filter
GRISM   = 'Old_Med_Red'                / grism
FOCUS   =                  801         / LDSS3-focus
COMMENT =                              / comment
SUBRASTR= '1:2032,1:4064'              / subraster x1:x2,y1:y2
TEMPCCD =               -110.0         / CCD temperature [C]
SOFTWARE= 'Version 1.06 (Aug 26 2007, 13:27:39)'
FITSVERS= '1.02'                       / FITS header version

        """


