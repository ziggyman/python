#!/usr/bin/python

""" This class provides a command line interface to the detector driver code.
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
from future.builtins import *

import logging as log
import time
import argparse
import sys
import datetime

from cassspectr.controller.controller import Controller
from cassspectr.controller.detector_controller import DetectorController
from cassspectr.interface.spectrograph_interface.ttypes import DetectorException

try:
    from pyds9 import ds9
except:
    from ds9 import ds9
    
from astropy.io import fits

welcome_str = """
==============================================
Welcome to the Detector Driver Remote CLI.
==============================================
"""

help_str = """
Options:
--------

  1  -- Initialize the detector
  1a -- Reset the connection to the driver
  2  -- Get the CCD temperature
  3  -- Get the cold finger temperature
  3a -- Get the utility board temperature
  3b -- Get the desired temperature
  3c -- Get the control temperature
  3d -- Set the control temperature
  4  -- Start an exposure
  5  -- Finish the exposure
  6  -- Abort the exposure
  7  -- Get the current exposure status
  8  -- Get the remaining time for the exposure 
  9  -- Get the number of pixels read out
  10 -- Get the image so far
  11 -- Reset the device
  12 -- Reset the controller
  13 -- Display the available gains
  14 -- Set the gain
  15 -- Set the readout speed
  16 -- Get the current gain index
  17 -- Get the current readout speed index
  18 -- Get a description of the acquisition parameters
  19 -- Get an engineering image (full CCD, no binning, no bias)
  20 -- Take an exposure and monitor its progress
  21 -- Test the data link
  22 -- Send ioctl command
  23 -- Read memory locations
  24 -- Get errors
  h  -- Show this help menu
  q  -- Quit or exit the program
"""

# The defaults for initializing the detector
nRows_init = 512    # (but lower half masked)
nCols_init = 2148   # 2048 + 50 prescan + 50 overscan

# The defaults for taking exposures
nRows = 256    # 512 rows, but lower half is masked
nCols = 2148
nRowBin = 1    # No binning
nColBin = 1    
nRowCen = 384
nColCen = 1074

def process_commandline_args():
    """Process the command line arguments."""
    parser = argparse.ArgumentParser(description="Driver program for the detector controller.")
    parser.add_argument('--log', 
                        help="Specify the logging level.", default='info', 
                        choices=['debug', 'info', 'warn', 'error', 'critical'])
    parser.add_argument('--logfile', 
                        help="Specify the log file.")
    parser.add_argument('--detector_host', dest='detector_host', help="Specify the detector tcp host.", default='localhost')
    parser.add_argument('--detector_port', dest='detector_port', help="Specify the detector tcp port.", default=9091)
    args = parser.parse_args()
    return args

def setup_logging(args):
    # Set up logging
    logfile = __file__.split('.py')[0] + ".log"
    if "logfile" in args and args.logfile != None:
        logfile = args.logfile
    print("Writing logs to {}".format(logfile))
    loglevel = args.log

    log.basicConfig(filename = logfile, 
                    format='%(asctime)s: %(message)s', 
                    level=getattr(log, loglevel.upper()))

def initialize(dc):
        """Initialize the detector"""
        if dc.is_initialized:
            print("Detector is already initialized")
        else:
            rows = input("Please enter the number of rows [{}]: ".format(nRows_init))
            if (len(rows) == 0):
                rows = nRows_init
            rows = int(rows)
            cols = input("Please enter the number of columns [{}]: ".format(nCols_init))
            if (len(cols) == 0):
                cols = nCols_init
            cols = int(cols)
            print("Initializing detector with {} rows and {} cols".format(rows, cols))
            dc.initialize(rows, cols)

def reset_connection(dc):
    """Reset the connection to the driver"""
    dc.reset_connection()

def get_ccd_temperature(dc):
    print("CCD temperature is {}".format(dc.get_ccd_temperature()))
    
def get_cold_finger_temperature(dc):
    print("Cold finger temperature is {}".format(dc.get_cold_finger_temperature()))

def get_utility_board_temperature(dc):
    print("Utility board temperature is {}".format(dc.get_utility_board_temperature()))

def get_desired_temperature(dc):
    print("Desired temperature is {}".format(dc.get_desired_temperature()))

def get_control_temperature(dc):
    print("Control temperature is {}".format(dc.get_control_temperature()))

def set_control_temperature(dc):
        control_temp = input("Please enter the control temperature in ADU: ")
        if (len(control_temp) == 0):
            print("Invalid control temperature")
            return
        control_temp = int(control_temp)
        dc.set_control_temperature(control_temp)

def start_exposure(controller):
    print("Starting exposure...")
    exp_time = input("Please enter the exposure time: ")
    if (len(exp_time) == 0):
        print("Invalid exposure time")
        return
    exp_time = float(exp_time)
    rows = input("Please enter the number of rows [{}]: ".format(nRows))
    if (len(rows) == 0):
        rows = nRows
    rows = int(rows)
    cols = input("Please enter the number of columns [{}]: ".format(nCols))
    if (len(cols) == 0):
        cols = nCols
    cols = int(cols)
    rowbin = input("Please enter the number of rows to bin [{}]: ".format(nRowBin))
    if (len(rowbin) == 0):
        rowbin = nRowBin
    rowbin = int(rowbin)
    colbin = input("Please enter the number of columns to bin [{}]: ".format(nColBin))
    if (len(colbin) == 0):
        colbin = nColBin
    colbin = int(colbin)
    rowcen = input("Please enter the centre row [{}]: ".format(nRowCen))
    if (len(rowcen) == 0):
        rowcen = nRowCen
    rowcen = int(rowcen)
    colcen = input("Please enter the centre column [{}]: ".format(nColCen))
    if (len(colcen) == 0):
        colcen = nColCen
    colcen = int(colcen)
    print("Starting exposure with time = {}, and {} {} {} {} {} {}".
          format(exp_time, rows, cols, rowbin, colbin, rowcen, colcen))
    supp_fits_info = {"OBSERVER": ("Carel", "Observer name"),
                      "OBJECT": ("Vega", "The object name")}
    controller.start_exposure(exp_time, 
                              rows, cols, 
                              rowbin, colbin, 
                              rowcen, colcen,
                              supp_fits_info)

def abort_exposure(dc):
    print("Aborting exposure...")
    dc.abort_exposure()
    print("Exposure aborted")

def finish_exposure(dc):
    print("Finishing exposure...")
    dc.finish_exposure()
    print("Exposure finished")

def get_exposure_status(dc):
    print("Getting acquisition status")
    status = dc.get_exposure_status()
    if status == 0:
        print("The acquisition status is: idle")
    elif status == 1:
        print("The acquisition status is: exposing")
    elif status == 2:
        print("The acquisition status is: reading out")
    else:
        print("Unknown status: {}".format(status))

def get_remaining_time(dc):
    time = dc.get_remaining_time()
    print("The remaining time is: {}".format(time))

def get_pixel_count(dc):
    count = dc.get_pixel_count()
    print("The number of pixels read is: {}".format(count))

def show_image(controller):
    filename = "tmp.fits"
    print("Generating fits...")
    print("Getting data...")
    arr = controller.dc.get_image_data()
    print("Generating fits with keys...")
    controller.gen_fits_with_keys(arr, filename)
    print("Displaying {}".format(filename))
    d = ds9()
    d.set("file {}".format(filename))
    
def reset_device(dc):
    print("Resetting device...")
    dc.reset_device()

def reset_controller(dc):
    print("Resetting controller...")
    dc.reset_controller()

def get_gains(dc):
    print("Gains:")
    gains = dc.get_gains()
    for readout_speed, val in gains.items():
        print("")
        print(" Readout speed index = {}".format(readout_speed))
        for gain_index, dict in val.items():
            print("    Gain index = {}".format(gain_index))
            for k, v in dict.items():
                print("           {} = {}".format(k, v))

def set_gain_index(dc):
    gain_id = input("Please enter the gain id [0|1]: ")
    dc.set_gain_index(int(gain_id))

def set_readout_speed_index(dc):
    readout_speed_id = input("Please enter the readout speed id [0,1]: ")
    dc.set_readout_speed_index(int(readout_speed_id))
    dc.set_gain_index(dc.get_current_gain_index())

def get_current_gain_index(dc):
    idx = dc.get_current_gain_index()
    gain_descrip = "Faint" if idx else "Bright"
    print("The current gain index is {}. That's {}".format(idx, gain_descrip))

def get_current_readout_speed_index(dc):
    idx = dc.get_current_readout_speed_index()
    readout_descrip = "Fast" if idx else "Slow"
    print("The current readout speed index is {}. That's {}".format(idx, readout_descrip))

def get_acquisition_description(dc):
    gn =  dc.get_current_gain_index()
    ro = dc.get_current_readout_speed_index()
    gain_descrip = "Faint" if gn else "Bright"
    readout_descrip = "Fast" if ro else "Slow"
    gains_table = dc.get_gains()
    print("The acquisition is {} and {}".format(gain_descrip, readout_descrip))
    print("The software gain is {}".format(gains_table[ro][gn]["sw_gain"]))
    print("The software speed is {}".format(gains_table[ro][gn]["sw_speed"]))
    print("The offset is {}".format(gains_table[ro][gn]["offset"]))
    print("The gain is {}".format(gains_table[ro][gn]["gain"]))
    print("The noise in ADU is {}".format(gains_table[ro][gn]["noise_adu"]))
    print("The noise in electrons is {}".format(gains_table[ro][gn]["noise_e"]))

def get_engineering_exposure(dc):
    exp_time = input("Please enter the exposure length: ")
    dc.start_engineering_exposure(float(exp_time))


def expose_and_monitor(controller, dc):
    status_map = {0: "Idle",
                  1: "Exposing",
                  2: "Readout"}
    start_exposure(controller)
    start_time = datetime.datetime.now()
    exp_end_time = start_time
    status = dc.get_exposure_status()
    oldstatus = status
    while status != 0 or oldstatus == 0:
        print("Time: {}  Status index = {} Status = {}".format(datetime.datetime.now(), 
                                                               status,
                                                               status_map[status]))
        if status != oldstatus:
            if status == 2:
                exp_end_time = datetime.datetime.now()
            oldstatus = status
        time.sleep(0.1)
        status = dc.get_exposure_status()
    read_end_time = datetime.datetime.now()
    print("Time: {}  Status index = {} Status = {}".format(datetime.datetime.now(), 
                                                           status,
                                                           status_map[status]))
    print("Exposure time was {}".format((exp_end_time - start_time).total_seconds()))
    print("Readout time was {}".format((read_end_time - exp_end_time).total_seconds()))

def test_data_link(dc):
    id = input("Please enter the board id: ")
    if id == None or id == '' or int(id) not in [1,2,3]:
        print("Incorrect board id entered. Must be 1, 2, or 3")
        return
    msg = input("Please enter the msg to send: ")
    if msg == None or msg == '':
        print("Incorrect msg entered. Must be an integer.")
        return
    try:
        msg = int(msg)
    except ValueError:
        print("Incorrect msg entered. Must be an integer.")
        return
    retval = dc.test_data_link(int(id), msg)
    if retval == 0:
        print("Call of test_data_link to board {} with message {} FAILED".format(id, msg))
    else:
        print("Call of test_data_link to board {} with message {} SUCCEEDED".format(id, msg))

ioctl = {1: "ASTROPCI_GET_HCTR",
         4: "ASTROPCI_GET_HSTR"}

def ioctl_commands():
    ioctl_str = "Please choose an ioctl command:\n"
    for k, v in ioctl.items():
        ioctl_str += "{}: {}\n".format(k,v)
    return ioctl_str

def ioctl_command(dc):
    ioctl_cmd = input(ioctl_commands())
    if ioctl_cmd and int(ioctl_cmd) in [1, 2, 3]:
        print("You chose option {}. Sending command {}".format(ioctl_cmd, ioctl[int(ioctl_cmd)]))
        result = dc.ioctl_command(ioctl_cmd)

def read_memory_locations(dc):
    id = input("Please enter the board id: ")
    if id == None or id == '' or int(id) not in [1,2,3]:
        print("Incorrect board id entered. Must be 1, 2, or 3")
        return
    memtype = input("Please enter the memory type [X|Y|P|R]: ")
    if memtype is None or memtype == '' or memtype not in ['X', 'Y', 'P', 'R']:
        print("Incorrect memory type entered. Must be one of X, Y, P or R")
        return
    memaddress = input("Please enter the memory address: ")
    try:
        memaddress = int(memaddress)
    except ValueError:
        print("Memory address must be an integer")
        return
    num_locations = input("Please enter the number of memory locations to read: ")
    try:
        num_locations = int(num_locations)
    except ValueError:
        print("Number of locations must be an integer")
        return
    memlocs = dc.read_memory_locations(id, memtype, memaddress, num_locations)
    i = 0
    for memloc in memlocs:
        print("At position [{}] read [{0X}] ([{}])".format(i, memloc, memloc))
        i = i+1

def get_errors(dc):
    errors = dc.get_errors()
    print("The following errors were reported by the dector:")
    for err in errors:
        print(err)

def run_event_loop(controller, dc):
    """Prompt for and process user input."""
    while True:
        try:
            print(help_str)
            action = input("Please enter an option: ")
            log.debug("Selected action = {}".format(action))
            if len(action) == 0:
                print("No response entered. Try again!")
                continue
            act = action.upper()
            if act == 'Q':
                break
            elif act == 'H':
                print(help_str.format(1,8))
            elif act == '1':
                initialize(dc)
            elif act == '1A':
                reset_connection(dc)
            elif act == '2':
                get_ccd_temperature(dc)
            elif act == '3':
                get_cold_finger_temperature(dc)
            elif act == '3A':
                get_utility_board_temperature(dc)
            elif act == '3B':
                get_desired_temperature(dc)
            elif act == '3C':
                get_control_temperature(dc)
            elif act == '3D':
                set_control_temperature(dc)
            elif act == '4':
                start_exposure(controller)
            elif act == '5':
                finish_exposure(dc)
            elif act == '6':
                abort_exposure(dc)
            elif act == '7':
                get_exposure_status(dc)
            elif act == '8':
                get_remaining_time(dc)
            elif act == '9':
                get_pixel_count(dc)
            elif act == '10':
                show_image(controller)
            elif act == '11':
                reset_device(dc)
            elif act == '12':
                reset_controller(dc)
            elif act == '13':
                get_gains(dc)
            elif act == '14':
                set_gain_index(dc)
            elif act == '15':
                set_readout_speed_index(dc)
            elif act == '16':
                get_current_gain_index(dc)
            elif act == '17':
                get_current_readout_speed_index(dc)
            elif act == '18':
                get_acquisition_description(dc)
            elif act == '19':
                get_engineering_exposure(dc)
            elif act == '20':
                expose_and_monitor(controller, dc)
            elif act == '21':
                test_data_link(dc)
            elif act == '22':
                ioctl_command(dc)
            elif act == '23':
                read_memory_locations(dc)
            elif act == '24':
                get_errors(dc)
            else:
                # Invalid action selected
                print("Invalid option entered. Try again!")
        except DetectorException as e:
            msg = e.message
            print("Caught DetectorException with message {}".format(msg))
        except Exception as e:
            type, value, traceback = sys.exc_info()
            print('Error type, str %s' % (type))
            print("Uh oh. Caught an exception: {}".format(e))

def main():
    print(welcome_str)
    # Basics
    args = process_commandline_args()
    setup_logging(args)
    log.info("Starting {}.".format(__file__))

    # Create the detector controller class
    log.debug("Instantiating detector controller.")
    try:
        print("Host: [{}]".format(args.detector_host))
        print("Port: [{}]".format(args.detector_port))
#        dc = DetectorController(detector_host=args.detector_host, detector_port=args.detector_port)
        controller = Controller(plc_host = "Empty", plc_port = 0, 
                                detector_host = args.detector_host, 
                                detector_port = args.detector_port)
        dc = controller.dc
    except Exception as e:
        print("Error initializing detector controller. Reason: " + str(e))
        sys.exit(1)

    # Run the main event loop. Control stays here until the user quits.
    run_event_loop(controller, dc)

    log.info("Exiting {}.".format(__file__))
    print("We're done... Goodbye.\n")

if __name__ == '__main__':
    main()
