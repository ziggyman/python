#!/usr/bin/env python
import sys

from thrift import Thrift
#from shared.ttypes import SharedStruct
from thrift.transport import TSocket
from thrift.transport import TTransport
from thrift.protocol import TBinaryProtocol
from thrift.server import TServer

from cassspectr.interface.spectrograph_interface.PLC import Processor
from cassspectr.plcmodel.plc_driver import PLCDriver

import logging as log
import argparse
import sys

def process_commandline_args():
    """Process the command line arguments."""
    parser = argparse.ArgumentParser(description="Driver program for the PLC controller.")
    parser.add_argument('--serial', default='/dev/ttyUSB1',
                        dest='serial_port',
                        help='Specify the serial port to use.')
    parser.add_argument('--log', 
                        help="Specify the logging level.", default='info', 
                        choices=['debug', 'info', 'warn', 'error', 'critical'])
    parser.add_argument('--test_mode', help="Enable test mode.", action="store_true")
    parser.add_argument('--tcp_port', default="9090", dest='tcp_port', help="Specify the TCP port.")
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

    log.basicConfig(filename = logfile, 
                    format='%(asctime)s: %(message)s', 
                    level=getattr(log, loglevel.upper()))


if __name__ == '__main__':
    # Basics
    args = process_commandline_args()
    setup_logging(args)
    log.info("Starting {}.".format(__file__))

    # Create the PLCDriver class
    log.debug("Instantiating PLCDriver class.")
    try:
        if args.test_mode:
            print("Starting in test mode")
            plcd = PLCDriver(port=args.serial_port, test=True)
        else:
            plcd = PLCDriver(port=args.serial_port)
            
    except Exception as e:
        print("Error initializing PLC driver. Reason: " + str(e))
        sys.exit(1)

    handler = plcd
    processor = Processor(handler)
    transport = TSocket.TServerSocket(port=args.tcp_port)
    tfactory = TTransport.TBufferedTransportFactory()
    pfactory = TBinaryProtocol.TBinaryProtocolFactory()

    server = TServer.TSimpleServer(processor, transport, tfactory, pfactory)

    print('Starting the server. Listening on port={}.'.format(args.tcp_port))
    server.serve()
    print("Server exiting... Done.")
