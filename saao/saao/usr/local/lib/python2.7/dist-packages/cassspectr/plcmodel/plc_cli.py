#!/usr/bin/env python

""" This class provides a command line interface to the PLCDriver code.

The actions possible include status updates, and setting of the various parameters
commands for the specified filter wheels.  

The port the filter wheel controller may be specified from the command
line (default is /dev/ttyUSB0).
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

from cassspectr.plcmodel.plc_driver import PLCDriver, PLCInvalidResultException
import logging as log
import time
import argparse
import sys

welcome_str = """
==================================================
Welcome to the Cassegrain Spectrograph PLC Driver.
==================================================
"""

help_str = """
Options:
--------

  h                              -- Show this help menu
  stat                           -- Display the plc driver controller status
  set <option>                   -- Set one of [gm_centred|gm_inbeam|fw_init|fw_reset|slit_init|slit_reset|grating_angle_init|grating_angle_reset]
  set <option> <boolean>         -- Set one of [arc_mirror|slit_shutter|ros_mirror|hartman_a|hartman_b] with a boolean value
  set <option> <value>           -- Set one of [fw_move|slit_width] with a value
  set <arc_lamp> <1|2> <boolean> -- Set arc_lamp 1 or 2 with a boolean value
  set <slit_illumination> <value> <boolean> -- Set the slit_illumination with a value and a boolean value
  q                              -- Quit or exit the program

In the above, boolean values may be chosen from <on|off> or <inbeam|outofbeam> (as appropriate).
"""

command_opts = ["gm_centred",
                "gm_inbeam",
                "fw_init",
                "fw_reset",
                "fw_move",
                "arc_mirror",
                "arc_lamp",
                "slit_init",
                "slit_reset",
                "slit_width",
                "slit_shutter",
                "slit_illumination",
                "rear_of_slit_mirror",
                "hartman_A",
                "hartman_B",
                "grating_angle_init",
                "grating_angle_reset"]


def process_commandline_args():
    """Process the command line arguments."""
    parser = argparse.ArgumentParser(description="Driver program for the filterwheel controller.")
    parser.add_argument('--serial', default='/dev/ttyUSB0',
                        dest='serial_port',
                        help='Specify the serial port to use.')
    parser.add_argument('--log', 
                        help="Specify the logging level.", default='info', 
                        choices=['debug', 'info', 'warn', 'error', 'critical'])
    parser.add_argument('--diagnostic_mode', help="Enable diagnostic mode.", action="store_true")
    parser.add_argument('--logfile', 
                        help="Specify the log file.")
    args = parser.parse_args()
    return args

def setup_logging(args):
    # Set up logging
    logfile = __file__.split('.py')[0] + ".log"
    if "logfile" in args and args.logfile != None:
        logfile = args.logfile
    print("Writing logs to {}".format(logfile))
    loglevel = args.log
    if args.diagnostic_mode:
        loglevel = "debug"

    log.basicConfig(filename = logfile, 
                    format='%(asctime)s: %(message)s', 
                    level=getattr(log, loglevel.upper()))

def get_status(plcd):
        """Update and display the current status of the filter wheel."""
        plcd.update_status()
        show_status(plcd)

def show_status(plcd):
    """Display the status of the PLC controller."""

    status_string = """
GMCentred             : {GMCentred}
GMInbeam              : {GMInbeam}
GMMoving              : {GMMoving}
GMFailure             : {GMFailure}
FilterInit            : {FilterInit}
FilterCentred         : {FilterCentred}
FilterMoving          : {FilterMoving}
FilterFailure         : {FilterFailure}
ARCMirror             : {ARCMirror}
ARC1                  : {ARC1}
ARC2                  : {ARC2}
SlitShutter           : {SlitShutter}
SlitIllumination      : {SlitIllumination}
RearOfSlitMirror      : {RearOfSlitMirror}
HartmanA              : {HartmanA}
HartmanB              : {HartmanB}
SlitWidthInitPos      : {SlitWidthInitPos}
SlitWidthMoving       : {SlitWidthMoving}
SlitWidthFailure      : {SlitWidthFailure}
GratingAngleInit      : {GratingAngleInit}
GratingAngleLimit1    : {GratingAngleLimit1}
GratingAngleLimit2    : {GratingAngleLimit2}
GratingAngleMoving    : {GratingAngleMoving}
GratingAngleFailure   : {GratingAngleFailure}
CameraFocusInit       : {CameraFocusInit}
CameraFocusLimit1     : {CameraFocusLimit1}
CameraFocusLimit2     : {CameraFocusLimit2}
CameraFocusMoving     : {CameraFocusMoving}
SlitWidthInitReq      : {SlitWidthInitReq}
AngleInitReq          : {AngleInitReq}
FilterwheelPosition   : {FilterwheelPosition}
GratingID             : {GratingID}
GratingInserted       : {GratingInserted}
GratingHatchClosed    : {GratingHatchClosed}
SlitShutterFailure    : {SlitShutterFailure}
ARCMirrorFailure      : {ARCMirrorFailure}
RoSMirrorFailure      : {RoSMirrorFailure}
HartmanFailure        : {HartmanFailure}
GratingAngleLo        : {GratingAngleLo}
GratingAngleHi        : {GratingAngleHi}
FocusPosition         : {FocusPosition}
FocusStepsLo          : {FocusStepsLo}
FocusStepsHi          : {FocusStepsHi}
SlitWidthPosition     : {SlitWidthPosition}
TopCrateInterlock     : {TopCrateInterlock}
FilterInterlock       : {FilterInterlock}
PneumaticsInterlock   : {PneumaticsInterlock}
ARCInterlock          : {ARCInterlock}
SlitWidthInterlock    : {SlitWidthInterlock}
BottomSignalInterlock : {BottomSignalInterlock}
BottomDriveInterlock  : {BottomDriveInterlock}"""
    
    state = plcd.state.get_state()
    print(status_string.format(**state))
    print("")

def run_event_loop(plcd):
    """Prompt for and process user input."""
    while True:
        action = input("Please enter an option: ")
        log.debug("Selected action = {}".format(action))
        if len(action) == 0:
            print("No response entered. Try again!")
            continue
        act = action[0].upper()
        if act == 'Q':
            break
        elif act == 'H':
            print(help_str)
            continue
        if action.upper()[0:4] == 'STAT':
            print("Getting status...")
            try:
                get_status(plcd)
            except PLCInvalidResultException:
                print("Uh oh - there was an error getting the status.")
        elif action.upper()[0:3] == 'SET':
            # TODO: Replace two or more spaces with one 
            args = action.split(' ')
            if len(args) < 2:
                print("You'll need to specify what you want to set. Please try again...")
                continue
            cmd = args[1]
            if cmd in ["fw_move", "slit_width"]:
                if len(args) < 3:
                    print("That choice requires a parameter...")
                    continue
                value = int(args[2])
                command = "set_{}".format(cmd)
                method = getattr(plcd, command)
                method(value)
            elif cmd in ["arc_mirror", "slit_shutter", "rear_of_slit__mirror", "hartman_a", "hartman_b"]:
                if len(args) < 3:
                    print("That choice requires a parameter...")
                    continue
                value = args[2]
                if args[2].upper() == 'INBEAM' or args[2].upper() == 'ON' or args[2].upper() == 'OPEN':
                    value = True
                else:
                    value = False
                command = "set_{}".format(cmd)
                method = getattr(plcd, command)
                method(value)
            elif cmd == "arc_lamp":
                if len(args) < 4:
                    print("That choice requires two parameters...")
                    continue
                value1 = int(args[2])
                if value1 not in [1,2]:
                    print("Please specify arc lamp 1 or 2")
                    continue
                if args[3].upper() == 'ON':
                    value2 = True
                else:
                    value2 = False
                command = "set_{}".format(cmd)
                method = getattr(plcd, command)
                method(value1, value2)
            elif cmd =="slit_illumination":
                if len(args) < 4:
                    print("That choice requires two parameters...")
                    continue
                value1 = int(args[2])
                if args[3].upper() == 'ON':
                    value2 = True
                else:
                    value2 = False
                command = "set_{}".format(cmd)
                method = getattr(plcd, command)
                method(value1, value2)
            elif cmd in command_opts:
                command = "set_{}".format(cmd)
                method = getattr(plcd, command)
                method()
            else:
                print("Unknown command {}".format(cmd))
        else:
            print("Unknown command {}".format(action))

if __name__ == '__main__':
    print(welcome_str)
    # Basics
    args = process_commandline_args()
    setup_logging(args)
    log.info("Starting {}.".format(__file__))

    # Create the PLCDriver class
    log.debug("Instantiating PLCDriver class.")
    try:
        plcd = PLCDriver(port=args.serial_port, test=True)
        # plcd = PLCDriver(port=args.serial_port)
    except Exception as e:
        print("Error initializing PLC driver. Reason: " + str(e))
        sys.exit(1)

    print(help_str)
    # Run the main event loop. Control stays here until the user quits.
    run_event_loop(plcd)

    # Tidy up before exit...
    log.debug("Closing serial port...")
    plcd.ser.close()
    log.info("Exiting {}.".format(__file__))
    print("We're done... Goodbye.\n")

