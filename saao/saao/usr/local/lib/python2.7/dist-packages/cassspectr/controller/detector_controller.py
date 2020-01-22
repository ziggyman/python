import logging as log
import time

from thrift import Thrift
from thrift.transport import TSocket
from thrift.transport import TTransport
from thrift.protocol import TBinaryProtocol

from cassspectr.interface.spectrograph_interface.DetectorService import Client
from numpy import array, reshape, uint16
from astropy.io import fits

class DetectorController:
    def __init__(self, tcp_host='localhost', tcp_port=9091):
        self.tcp_host = tcp_host
        self.tcp_port = tcp_port
        self.open_connection()
        self.is_initialized = False

    def open_connection(self):
        try:
            self.socket = TSocket.TSocket(self.tcp_host, self.tcp_port)
            self.transport = TTransport.TBufferedTransport(self.socket)
            protocol = TBinaryProtocol.TBinaryProtocol(self.transport)
            self.client = Client(protocol)
            print("Opening transport")
            self.transport.open()
            print("Transport opened")
        except Thrift.TException as tx:
            print("Instantiation of DetectorController failed: {}".format(tx.message))
            exit(1)

    def reset_connection(self):
        self.transport.close()
        self.socket.close()
        self.open_connection()
        self.is_initialized = False

    def __exit__(self):
        self.transport.close()

    def initialize(self, numrows, numcols):
        """Initialize the detector.
        numrows and numcols are the actual number of rows and columsn on
        the detector, not just those in the visible area.
        """
        self.client.initialise(numrows, numcols)
        self.is_initialized = True

    def get_ccd_temperature(self):
        """Return the temperature of the CCD."""
        return self.client.getCCDTemperature()

    def get_cold_finger_temperature(self):
        """Return the cold plate temperature."""
        return self.client.getColdFingerTemperature()

    def get_utility_board_temperature(self):
        """Return the utility board temperature."""
        return self.client.getUtilityBoardTemperature()

    def get_desired_temperature(self):
        """Return the desired temperature."""
        return self.client.getDesiredTemperature()

    def get_control_temperature(self):
        """Return the control temperature."""
        return self.client.getControlTemperature()

    def set_control_temperature(self, temp):
        """Set the control temperature"""
        print("Calling driver's setControlTemperature({})".format(temp))
        self.client.setControlTemperature(temp)

    def start_exposure(self, exp_time, rows, cols, rowbin, colbin, rowcen, colcen):
        """Start an exposure.

        The exposure time, detector dimensions, binning and subframe
        are specified.  This call is non-blocking. To halt the exposure,
        call finish_exposure (retains data) or abort_exposure (discards
        data). To monitor the exposure status, use get_exposure_status(),
        get_remaining_time() or get_pixel_count().
        """
        self.nRows = rows
        self.nCols = cols
        self.expTime = exp_time
        self.rowBin = rowbin
        self.colBin = colbin
        print("Calling startExposure with {} {} {} {} {} {} {}".
              format(exp_time, 
                     rows, cols, 
                     rowbin, colbin, 
                     rowcen, colcen))
        self.client.startExposure(exp_time, 
                                  rows, cols, 
                                  rowbin, colbin, 
                                  rowcen, colcen)

    def start_engineering_exposure(self, exp_time):
        self.client.startEngineeringExposure(exp_time)

    def finish_exposure(self):
        self.client.finishExposure()

    def stop_exposure(self):
        self.client.stopExposure()

    def abort_exposure(self):
        print("Controller calling abortExposure")
        self.client.abortExposure()

    def get_exposure_status(self):
        return self.client.getExposureStatus()

    def get_remaining_time(self):
        return self.client.getRemainingTime()

    def get_pixel_count(self):
        return self.client.getPixelCount();

    def reset_device(self):
        self.client.resetDevice()

    def reset_controller(self):
        self.client.resetController()

    def test_data_link(self, id, msg):
	print("Controiller got TDL request")
        self.client.testDataLink(id, msg)

    def ioctl_command(self, ioctl_cmd):
        return self.client.ioctlCommand(ioctl_cmd)

    def read_memory_locations(self, board_id, memtype, memaddress, num_locations):
        return self.client.readMemoryLocations(board_id, memtype, memaddress, num_locations)

    def get_gains(self):
        gains =  self.client.getGains()
        return gains

    def set_gain_index(self, gain_id):
        self.client.setGainIndex(gain_id)

    def set_readout_speed_index(self, readout_speed_id):
        self.client.setReadoutSpeedIndex(readout_speed_id)

    def get_current_gain_index(self):
        return self.client.getCurrentGainIndex()

    def get_current_readout_speed_index(self):
        return self.client.getCurrentReadoutSpeedIndex()

    def get_image_data(self):
        """Return an array of the data.

        The array is shaped to the number of rows and columns. The array
        includes only the pixels read out so far, the remainder of
        elements are set to zero.
        """
        data = self.client.getImageData()
        lendata = len(data)
        print("Length of data: {}  nrows = {}  ncols = {}".format(lendata, self.nRows, self.nCols))
        print("Binning: rowBin = {}  colBin = {}".format(self.rowBin, self.colBin))
        data += [0] * (self.nRows/self.rowBin * self.nCols/self.colBin - len(data))
        arr = array(data, dtype=uint16)
        arr = arr.reshape(self.nRows/self.rowBin, self.nCols/self.colBin)
        return arr

    def add_detector_info(self, dict):
        """Return a dictionary of exposure parameters."""
        dict["EXPOSURE"] = (self.expTime, "Total exposure time")
        dict["HBIN"] = (self.rowBin, "Horizontal binning")
        dict["VBIN"] = (self.colBin, "Vertical binning")
        # TODO: Add subrect
        # Add gains info:
        gn =  self.get_current_gain_index()
        ro = self.get_current_readout_speed_index()
        gains_table = self.get_gains()
        dict["GAINVAL"] = (gains_table[ro][gn]["gain"], "Gain value") 
        dict["NOISEADU"] = (gains_table[ro][gn]["noise_adu"], "Noise (ADU)") 
        dict["NOISEEL"] = (gains_table[ro][gn]["noise_e"], "Noise (electrons)") 
        

    def get_errors(self):
        return self.client.getErrors()

if __name__ == '__main__':
    try:
        nrows = 16
        ncols = 1024
        rfact = 1
        cfact = 1
        shutter_open = True
        temperature = -8
        exposure_time = 10
        
        dc = DetectorController()
        dc.initialize(nrows, ncols)
#        dc.set_temperature(temperature)
#        print("The detector temperature is {}".format(dc.get_temperature()))
        print("Starting exposure...")
        dc.start_exposure(exposure_time, nrows, ncols, rfact, cfact, shutter_open)
        print("Done starting exposure")
        remaining_time = dc.get_remaining_time()
        print("Remaining time = {}".format(remaining_time));
        print("Pixel count  = {}".format(dc.get_pixel_count()));
        while remaining_time < exposure_time:
            time.sleep(1)
            remaining_time = dc.get_remaining_time()
            print("Remaining time = {}".format(remaining_time));
            print("Pixel count  = {}".format(dc.get_pixel_count()));
    except Exception as e:
        print("Uh oh. Caught an exception: {}".format(str(e)))
    print("Done")
