#!/usr/bin/python

""" This class provides a command line interface to the plc driver code.
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
from future.builtins import *

import logging as log
import time
import argparse
import sys
from math import floor

from cassspectr.controller.controller import Controller
from cassspectr.controller.plc_controller import PLCController
from thrift.Thrift import TType, TMessageType, TException, TApplicationException
from thrift.transport.TTransport import TTransportException

fw_min = 1
fw_max = 8

welcome_str = """
==============================================
Welcome to the PLC Driver Remote CLI.
==============================================
"""

help_str = """
Options:
--------

  1  -- Get and display status                 1a -- Get failed states                    1b -- Reset the connection to the driver

  2  -- Set guide mirror centred               3  -- Set guide mirror inbeam              4  -- Set rear-of_slit mirror in beam/out of beam
  5  -- Move ARC mirror in beam/out of beam    6  -- Switch ARC lamp 1 on/off             7  -- Switch ARC lamp 2 on/off

  8  -- Reset filterwheel                      9  -- Initialize filterwheel               10 -- Move filterwheel

  11 -- Reset slit                             12 -- Initialize slit                      13 -- Set slit width
  14 -- Open/close slit shutter                15 -- Set slit illumination on/off         16 -- Set slit illumination value

  17 -- Set grating angle abs                  18 -- Set grating angle rel
  19 -- Reset grating angle                    20 -- Initialize grating angle             

  22 -- Set camera focus abs                   23 -- Set camera focus rel
  24 -- Reset camera focus                     25 -- Initialize camera focus              

  27 -- Set Hartman A in beam/out of beam      28 -- Set Hartman B in beam/out of beam

  29 -- Reconnect the serial port

  h  -- Show this help menu                    q  -- Quit or exit the program
"""

def get_status(plcc):
    status = plcc.get_status()
    display_status(status)

def display_status(status):
    """Display the status of the PLC controller."""

    status_string = """
======================================================  STATUS =================================================================

FilterwheelPosition   : {FilterwheelPosition}        FilterInit            : {FilterInit}       FilterCentred         : {FilterCentred}
FilterMoving          : {FilterMoving}        FilterFailure         : {FilterFailure}

SlitWidthPosition     : {SlitWidthPosition}        SlitIllumination      : {SlitIllumination}       SlitIlluminationValue : {SlitIlluminationValue}
SlitWidthInitPos      : {SlitWidthInitPos}        SlitWidthInitReq      : {SlitWidthInitReq}       SlitWidthMoving       : {SlitWidthMoving}        SlitWidthFailure      : {SlitWidthFailure}
SlitShutter           : {SlitShutter}        SlitShutterFailure    : {SlitShutterFailure}

GratingAngleSteps     : {GratingAngleSteps:06d}   GratingAngle          : {GratingAngle:03.2f}   GratingAngleInit      : {GratingAngleInit}
AngleInitReq          : {AngleInitReq}        GratingAngleMoving    : {GratingAngleMoving}       GratingAngleLimit1    : {GratingAngleLimit1}        GratingAngleLimit2    : {GratingAngleLimit2}
GratingID             : {GratingID:05d}    GratingInserted       : {GratingInserted}       GratingHatchClosed    : {GratingHatchClosed}        GratingAngleFailure   : {GratingAngleFailure}

CameraFocusInit       : {CameraFocusInit}        CameraFocusMoving     : {CameraFocusMoving}       CameraFocusLimit1     : {CameraFocusLimit1}        CameraFocusLimit2     : {CameraFocusLimit2}
FocusAtPosition       : {FocusAtPosition}        FocusPosition         : {FocusPosition}     FocusLVDT             : {FocusPositionPot}

GMCentred             : {GMCentred}        GMInbeam              : {GMInbeam}       GMMoving              : {GMMoving}        GMFailure             : {GMFailure}

RearOfSlitMirror      : {RearOfSlitMirror}        RoSMirrorFailure      : {RoSMirrorFailure}

ARCMirror             : {ARCMirror}        ARC1                  : {ARC1}       ARC2                  : {ARC2}        ARCMirrorFailure      : {ARCMirrorFailure}

HartmanA              : {HartmanA}        HartmanB              : {HartmanB}       HartmanFailure        : {HartmanFailure}

TopCrateInterlock     : {TopCrateInterlock}        FilterInterlock       : {FilterInterlock}       PneumaticsInterlock   : {PneumaticsInterlock}         ARCInterlock          : {ARCInterlock}
SlitWidthInterlock    : {SlitWidthInterlock}        BottomSignalInterlock : {BottomSignalInterlock}       BottomDriveInterlock  : {BottomDriveInterlock}

================================================================================================================================"""
    
    print(status_string.format(**status))
    print("")

def get_errors(plcc):
    failed_states = plcc.get_errors()
    print("Failed states: ")
    for failed in failed_states:
        print("  {}".format(failed))
    print("")

def reset_connection(plcc):
    plcc.reset_connection()

def set_gm_centred(plcc):
    plcc.set_gm_centred()

def set_gm_inbeam(plcc):
    plcc.set_gm_inbeam()

def set_fw_init(plcc):
    plcc.set_fw_init()

def set_fw_reset(plcc):
    plcc.set_fw_reset()

def get_response(text, response_type, response_range=None, default=None):
    response = input(text)
    if len(response) == 0:
        if default:
            return default
        else:
            return None
    try:
        if response_type == 'i' or response_type == 'b':
            response = int(response)
        elif response_type == 'f':
            response = float(response)
    except ValueError:
        print("The value entered was not of the expected type")
        return None
    if response_range is not None:
        try:
            assert(response in response_range)
        except AssertionError:
            print("The value entered was not in the expected range {}".format(response_range))
            return None
    if response_type == 'b':
        if int(response) == 1:
            response = True
        elif int(response) == 0:
            response = False
    return response
            

def set_fw_move(plcc):
    pos = get_response("Please enter the required filterwheel position: ", 'i', range(fw_min, fw_max+1))
    if pos is not None:
        plcc.set_fw_move(pos)

def set_arc_mirror(plcc):
    ib = get_response("Move the arc mirror in beam (1) or out of beam (0): ", 'b', [0,1])
    print("Arc Mirror IB = {}".format(ib))
    if ib is not None:
        plcc.set_arc_mirror(ib)

def set_arc_lamp_1(plcc):
    on = get_response("Switch arc lamp 1 on (1) or off (0): ", 'b', [0,1])
    if on is not None:
        plcc.set_arc_lamp_1(on)

def set_arc_lamp_2(plcc):
    on = get_response("Switch arc lamp 2 on (1) or off (0): ", 'b', [0,1])
    if on is not None:
        plcc.set_arc_lamp_2(on)

def set_slit_init(plcc):
    plcc.set_slit_init()

def set_slit_reset(plcc):
    plcc.set_slit_reset()

def set_slit_width(plcc):
    width = get_response("Please enter the required slit width: ", "f")
    if width is not None:
        plcc.set_slit_width(width)

def set_slit_shutter(plcc):
    shutterstatus = get_response("Open (1) or close (0) the slit shutter: ", 'b', [0,1])
    if shutterstatus is not None:
        plcc.set_slit_shutter(shutterstatus)

def set_slit_illumination_status(plcc):
    status = get_response("Set slit illumination status on (1) or off (0): ", 'b', [0,1])
    if status is not None:
        plcc.set_slit_illumination_status(status)

def set_slit_illumination_value(plcc):
    value = get_response("Set slit illumination value: ", 'f')
    if value is not None:
        plcc.set_slit_illumination_value(value)

def set_rear_of_slit_mirror(plcc):
    status = get_response("Set rear of slit mirror in-beam (1) or out-of-beam (0): ", 'b', [0,1])
    if status is not None:
        plcc.set_rear_of_slit_mirror(status)
    
def set_grating_angle_init(plcc):
    plcc.set_grating_angle_init()

def set_grating_angle_reset(plcc):
    plcc.set_grating_angle_reset()

def grating_angle_move_abs(plcc):
    value = get_response("Please enter the required grating angle as a decimal: ", 'f')
    if value is not None:
        plcc.grating_angle_move_abs(value)

def grating_angle_move_rel(plcc):
    value = get_response("Please enter the relative grating move angle as a decimal: ", 'f')
    if value is not None:
        plcc.grating_angle_move_rel(value)

def set_camera_focus_init(plcc):
    plcc.set_camera_focus_init()

def set_camera_focus_reset(plcc):
    plcc.set_camera_focus_reset()

def camera_focus_move_abs(plcc):
    value = get_response("Please enter the focus value as a decimal: ", 'f')
    if value is not None:
        plcc.camera_focus_move_abs(value)

def camera_focus_move_rel(plcc):
    value = get_response("Please enter the relative focus move value as a decimal: ", 'f')
    if value is not None:
        plcc.camera_focus_move_rel(value)

def set_hartman_a(plcc):
    in_beam = get_response("Move the Hartman A lamp in beam (1) or out of beam (0): ", 'b', [0,1])
    if in_beam is not None:
        plcc.set_hartman_a(in_beam)

def set_hartman_b(plcc):
    in_beam = get_response("Move the Hartman B lamp in beam (1) or out of beam (0): ", 'b', [0,1])
    if in_beam is not None:
        plcc.set_hartman_b(in_beam)

def reconnect_serial(plcc):
    plcc.reconnect_serial()

def process_commandline_args():
    """Process the command line arguments."""
    parser = argparse.ArgumentParser(description="Driver program for the camera controller.")
    parser.add_argument('--log', 
                        help="Specify the logging level.", default='info', 
                        choices=['debug', 'info', 'warn', 'error', 'critical'])
    parser.add_argument('--logfile', 
                        help="Specify the log file.")
    parser.add_argument('--tcp_host', dest='tcp_host', help="Specify the tcp host.", default='localhost')
    parser.add_argument('--tcp_port', dest='tcp_port', help="Specify the tcp port.", default="9090")
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

def initialize(plcc):
        """Initialize the plc"""
        print("Calling controller's initialize method")
        plcc.initialize()

    

def run_event_loop(plcc):
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
                get_status(plcc)
            elif act == '1A':
                get_errors(plcc)
            elif act == '1B':
                reset_connection(plcc)
            elif act == '2':
                set_gm_centred(plcc)
            elif act == '3':
                set_gm_inbeam(plcc)
            elif act == '4':
                set_rear_of_slit_mirror(plcc)
            elif act == '5':
                set_arc_mirror(plcc)
            elif act == '6':
                set_arc_lamp_1(plcc)
            elif act == '7':
                set_arc_lamp_2(plcc)
            elif act == '8':
                set_fw_reset(plcc)
            elif act == '9':
                set_fw_init(plcc)
            elif act == '10':
                set_fw_move(plcc)
            elif act == '11':
                set_slit_reset(plcc)
            elif act == '12':
                set_slit_init(plcc)
            elif act == '13':
                set_slit_width(plcc)
            elif act == '14':
                set_slit_shutter(plcc)
            elif act == '15':
                set_slit_illumination_status(plcc)
            elif act == '16':
                set_slit_illumination_value(plcc)
            elif act == '17':
                grating_angle_move_abs(plcc)
            elif act == '18':
                grating_angle_move_rel(plcc)
            elif act == '19':
                set_grating_angle_reset(plcc)
            elif act == '20':
                set_grating_angle_init(plcc)
            elif act == '22':
                camera_focus_move_abs(plcc)
            elif act == '23':
                camera_focus_move_rel(plcc)
            elif act == '24':
                set_camera_focus_reset(plcc)
            elif act == '25':
                set_camera_focus_init(plcc)
            elif act == '27':
                set_hartman_a(plcc)
            elif act == '28':
                set_hartman_b(plcc)
            elif act == '29':
                reconnect_serial(plcc)
            else:
                # Invalid action selected
                print("Invalid option entered. Try again!")
        except TTransportException as e:
            print("Hmm, caught thrift transport exception. Is the PLC server still up?")
        except Exception as e:
            print("Uh oh. Caught an exception: {}".format(e.message))

def main():

    print(welcome_str)
    # Basics
    args = process_commandline_args()
    setup_logging(args)
    log.info("Starting {}.".format(__file__))

    # Create the camera controller class
    log.debug("Instantiating PLC controller at {}:{}.".format(args.tcp_host, args.tcp_port))
    try:
        plcc = PLCController(tcp_host=args.tcp_host, tcp_port=int(args.tcp_port))
    except Exception as e:
        print("Error initializing plc controller. Reason: " + str(e))
        sys.exit(1)
    log.debug("Done instantiating PLC controller at {}:{}.".format(args.tcp_host, args.tcp_port))

    # Run the main event loop. Control stays here until the user quits.
    run_event_loop(plcc)

    log.info("Exiting {}.".format(__file__))
    print("We're done... Goodbye.\n")

if __name__ == '__main__':
    main()
