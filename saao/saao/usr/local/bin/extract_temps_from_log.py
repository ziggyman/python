#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
from future.builtins import *
import sys
import time
import datetime

ccd_string = "CCD temp. ="
coldfinger_string = "cold finger temp. ="

def get_file_names():
    today =  time.strftime("%Y%m%d")
    if len(sys.argv) == 2:
        infile = sys.argv[1]
    else:
        infile = input("Please enter the logfile name [cassspectr-ui{}.log]: ".format(today))
        if not infile or infile is '':
            infile = "cassspectr-ui{}.log".format(today)
    outfile = input("Please enter the output filename [temperatures_{}.txt]: ".format(today))
    if not outfile or outfile is '':
        outfile = "temperatures_{}.log".format(today)
    return infile, outfile

def process_log(infile, outfile):
    append = input("Do you want to append the temperatures to the existing output file [Y/n]? ")
    if not append or append == '':
        append = 'a'
    elif append.upper() not in ['Y', 'N']:
        print("Invalid response. Will append data.")
        append = 'a'
    elif append.upper() == 'Y':
        append = 'a'
    elif append.upper() == 'N':
        append = 'w'
    else:
        print("Huh? we shouldn't be here.")
        sys.exit()
    with open(infile, "r") as fin:
        with open(outfile, append) as fout:
            ccd_date = None
            cf_date = None
            for line in fin:
                if ccd_string in line:
                    ccd_dt, ccd_temp = line.strip().split(ccd_string)
                    print("Date = {}  CCD temp = {}".format(ccd_dt, ccd_temp))
                    ccd_dt = ccd_dt.split(',')[0]
                    ccd_date = datetime.datetime.strptime(ccd_dt, "%Y-%m-%d %H:%M:%S")
                elif coldfinger_string in line:
                    cf_dt, cf_temp = line.strip().split(coldfinger_string)
                    print("Date = {}  Coldfinger temp = {}".format(cf_dt, cf_temp))
                    cf_dt = cf_dt.split(',')[0]
                    cf_date = datetime.datetime.strptime(cf_dt, "%Y-%m-%d %H:%M:%S")
                #print("CCD date: {}  CF date == {}".format(ccd_date, cf_date))
                if (ccd_date is not None and cf_date is not None):
                    tdelta = cf_date - ccd_date
                    if tdelta.seconds < 10:
                        print("{} {} {} Delta = {}".format(ccd_date, ccd_temp, cf_temp, tdelta.seconds))
                        fout.write("{} {} {}\n".format(ccd_date, ccd_temp, cf_temp))
                        ccd_date = None

def run():
    infile, outfile = get_file_names()
    print("Input file name {0}  Output file name = {1}".format(infile, outfile))
    process_log(infile, outfile)

if __name__ == '__main__':
    run()
