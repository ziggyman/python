
import logging as log
import time

from thrift import Thrift
from thrift.transport import TSocket
from thrift.transport import TTransport
from thrift.protocol import TBinaryProtocol

from cassspectr.interface.spectrograph_interface.PLC import Client

class PLCController:
    def __init__(self, tcp_host='localhost', tcp_port=9090):
        self.tcp_host = tcp_host
        self.tcp_port = tcp_port
        self.open_connection()

    def open_connection(self):
        try:
            self.socket = TSocket.TSocket(self.tcp_host, self.tcp_port)
            self.transport = TTransport.TBufferedTransport(self.socket)
            protocol = TBinaryProtocol.TBinaryProtocol(self.transport)
            self.client = Client(protocol)
            print("Opening transport")
            self.transport.open()
            print("Transport opened")
        except Thrift.TException, tx:
            print '%s' % (tx.message)
            exit(1)

    def reset_connection(self):
        self.transport.close()
        self.socket.close()
        self.open_connection()

    def __exit__(self):
        self.transport.close()

    def reconnect_serial(self):
        self.client.reconnect_serial()

    def add_plc_info(self, dict):
        pass
        # TODO: add the PLC info we need to the dict. The format is 
        # {key1: (value1, comment1), ...}
        # So: dict[key1] = (value1, comment1)...

        # status = self.get_status()
        # dict['FOCUSPOS'] = (status['FocusPosition'], 'Focus position (mm)')
        # dict['GR-ANGLE'] = (status['GratingAngle'], 'GRATING ANGLE')
        # slitString = status['SlitWidthPosition']
        # slitwidArcsec = slitString.split(' ')[-1].replace('"','')
        # slitwidPos = slitString.split(' ')[0].replace(':','')
        # dict['SLITPOS'] = (str(slitwidPos) ,'SLIT WIDTH (index)')
        # dict['SLITWID'] = (str(slitwidArcsec) ,'SLIT WIDTH (arcsec)')
        # dict['FILTER'] = (str(self.ui.comboBoxFilter.currentText()), 'FILTER NAME')
        # # CCD Pane:
        # dict['FRAME'] = (self.settings['frameNum'], 'FRAME NUMBER')
        # dict['OBJECT'] = (str(self.ui.displExpName.text()), 'Target name')
        # actualExptime = str(self.ui.displExpTime.text())
        # print '*** Need to take actual time exposed, not requested here ***' # *** TODO
        # dict['EXPTIME'] = (actualExptime, 'EXPOSURE TIME (s)')
        # dict['EXPTYPE'] = (self.settings['expType'], 'EXPOSURE TYPE')
        # dict['BINNING'] = (self.settings['ccdBinning'], 'CCD binning')
        # dict['CCDSUM'] = (self.settings['ccdBinning'].replace('x',' '), 'CCD binning')
        # dict['SGAIN'] = (str(self.ui.comboCCDmodeGain.currentText()), 'CCD Mode GAIN')
        # dict['SRDNOISE'] = (str(self.ui.comboCCDmodeNoise.currentText()), 'CCD Mode RDNOISE')
        # dict['CCDTEMP'] = (str(self.ui.labCCDTemp.text()), 'CCD temperature')
        # dict['CF-TEMP'] = (str(self.ui.labColdFingerTemp.text()), 'Cold Finger Temperature')
        # dict['ARC-LAMP'] = (self.settings['ARC'], 'arc lamp ID')
        # dict['FOCUSPOS'] = (self.settings['focusPosition'], 'Focus position (mm)')
        # dict['HARTPOS'] = (self.settings['HartmannState'], 'Position of Hartmann Shutter')

    def get_status(self):
        """Update and display the current status of the filter wheel."""
        log.debug("Fetching status...")
        status = self.client.get_status()
        status = self.clean_status(status)
        return status
        
    def get_errors(self):
        """Return a list of failed plc controls"""
        return self.client.get_errors()

    def set_gm_centred(self):
        self.client.set_gm_centred()

    def set_gm_inbeam(self):
        self.client.set_gm_inbeam()

    def set_fw_init(self):
        self.client.set_fw_init()

    def set_fw_reset(self):
        self.client.set_fw_reset()

    def set_fw_move(self, pos):
        self.client.set_fw_move(pos)

    def set_arc_mirror(self, in_beam):
        self.client.set_arc_mirror(in_beam)

    def set_arc_lamp_1(self, on):
        self.client.set_arc_lamp_1(on)

    def set_arc_lamp_2(self, on):
        self.client.set_arc_lamp_2(on)

    def set_slit_init(self):
        self.client.set_slit_init()

    def set_slit_reset(self):
        self.client.set_slit_reset()

    def set_slit_width(self, width):
        self.client.set_slit_width(width)

    def set_slit_shutter(self, shutterstatus):
        self.client.set_slit_shutter(shutterstatus)

    def set_slit_illumination_status(self, status):
        if status == 0:
            status = False
        else:
            status = True
        self.client.set_slit_illumination_status(status)

    def set_slit_illumination_value(self, value):
        self.client.set_slit_illumination_value(value)

    def set_rear_of_slit_mirror(self, status):
        if status == 0:
            status = False
        else:
            status = True
        self.client.set_rear_of_slit_mirror(status)
    
    def set_grating_angle_init(self):
        self.client.set_grating_angle_init()

    def set_grating_angle_reset(self):
        self.client.set_grating_angle_reset()

    def grating_angle_move_abs(self, value):
        self.client.grating_angle_move_abs(value)

    def grating_angle_move_rel(self, value):
        self.client.grating_angle_move_rel(value)

    # def set_grating_angle_rel(self):
    #     self.client.set_grating_angle_rel()

    # def set_grating_angle(self, angle):
    #     self.client.set_grating_angle(angle)

    # def set_grating_angle_steps(self, steps):
    #     self.client.set_grating_angle_steps(steps)

    def set_camera_focus_init(self):
        self.client.set_camera_focus_init()

    def set_camera_focus_reset(self):
        self.client.set_camera_focus_reset()

    # def set_camera_focus_abs(self):
    #     self.client.set_camera_focus_abs()

    # def set_camera_focus_rel(self):
    #     self.client.set_camera_focus_rel()

    # def set_camera_focus_value(self, val):
    #     self.client.set_camera_focus_value(val)

    def camera_focus_move_abs(self, pos):
        self.client.camera_focus_move_abs(pos)

    def camera_focus_move_rel(self, val):
        self.client.camera_focus_move_rel(val)

    def set_hartman_a(self, in_beam):
        self.client.set_hartman_a(in_beam)

    def set_hartman_b(self, in_beam):
        self.client.set_hartman_b(in_beam)

    def clean_status(self, status):
        # Thrift returns the status as a map of strings to floats. But
        # some of these are really ints, and this function converts
        # them to ints in the status dictionary
        status_ints = ["GMCentred",
                       "GMInbeam",
                       "GMMoving",
                       "GMFailure",
                       "FilterInit",
                       "FilterCentred",
                       "FilterMoving",
                       "FilterFailure",
                       "ARCMirror",
                       "ARC1",
                       "ARC2",
                       "SlitShutter",
                       "SlitIllumination",
                       "SlitIlluminationValue",
                       "RearOfSlitMirror",
                       "HartmanA",
                       "HartmanB",

                       "SlitWidthInitPos",
                       "SlitWidthAtLimits",
                       "SlitWidthMoving",
                       "SlitWidthFailure",
                       "GratingAngleInit",
                       "GratingAngleLimit1",
                       "GratingAngleLimit2",
                       "GratingAngleMoving",
                       "GratingAngleFailure",
                       "CameraFocusInit",
                       "CameraFocusLimit1",
                       "CameraFocusLimit2",
                       "CameraFocusMoving",
                       "CameraFocusFailure",
                       "SlitWidthInitReq",
                       "AngleInitReq",

                       "FilterwheelPosition",
                       "GratingID",
                       "GratingInserted",
                       "GratingHatchClosed",
                       "SlitShutterFailure",
                       "ARCMirrorFailure",
                       "RoSMirrorFailure",
                       "HartmanFailure",

                       "GratingAngleSteps",

                       "SlitWidthPosition",

                       "TopCrateInterlock",
                       "FilterInterlock",
                       "PneumaticsInterlock",
                       "ARCInterlock",
                       "SlitWidthInterlock",
                       "BottomSignalInterlock",
                       "BottomDriveInterlock",
                       "SlitWidthErrorCodes",
                       "GratingAngleErrorCodes",
                       "FocusError",
                       "FocusMotorOn",
                       "FocusMoving",
                       "FocusReferencing",
                       "FocusAtPosition",
                       "FocusNegativeLimit",
                       "FocusReference",
                       "FocusPositiveLimit",
                       "FocusIn1",
                       "FocusIn2",
                       "FocusIn3",
                       "FocusIn4"]

        for key in status_ints:
            status[key] = int(status[key])
        return status

if __name__ == '__main__':
    plcc = PLCController()
    print("Status = {}".format(plcc.get_status()))
    print("")
    print("Setting gm_centred")
    plcc.client.set_gm_centred()
    status = plcc.get_status()
    print("GM_Centred = [{}]".format(status["GMCentred"]))
    print("GM_Inbeam = [{}]".format(status["GMInbeam"]))
    print("")
    print("Setting gm_inbeam")
    plcc.client.set_gm_inbeam()
    status = plcc.get_status()
    print("GM_Centred = [{}]".format(status["GMCentred"]))
    print("GM_Inbeam = [{}]".format(status["GMInbeam"]))
    plcc.client.set_slit_init()
    plcc.client.set_slit_reset()
    print("Sleeping")
    time.sleep(2)
    print("Done sleeping. Setting slit width")
    plcc.client.set_slit_width(2000)
    print("Done")
