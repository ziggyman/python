#! /usr/bin/env python

#  
#  Code Name: sedm_createbias
#  
#  Brief Description: 
#     Run the zerocombine in pyRAF with input bias image list to create master bias image 
#  
#  Required Inputs:
#     Bias image list
#     name of output master bias image
#
#  Author: 
#     Chow-Choong Ngeow (National Central University)
#
#  Date Created:
#     15 Feb 2012
#
#  Date Last Modify:
#     15 Feb 2012
#

import sys, os, getopt
from pyraf import iraf
from pyraf.iraf import imred, ccdred

# global parameter for versioning
version = 1.0

def print_help(program_name):
    print program_name, " (Version:",version,")"
    print "Usage:", program_name, "-l <bias_img.list> -o <master_bias.fits> [OPTION]"
    print "   Option:"
    print "       -c <combine type>   Combine type (average|median) [default: median]"
    print "       -r <reject type>    Reject type (none|minmax|ccdclip|crreject|sigclip|avsigclip|pcli) [default: minmax]"
    print "       -v                  Verbose"
    print "       -h                  Print help page"

def main():
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hl:o:c:r:v", ["help", "list=", "output=", "combine=", "reject="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print sys.argv[0],"ERROR:",str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    biaslist = None
    biasfits = None
    combinetype = "median"
    rejecttype = "minmax"
    verbose = False

    # print out the usage
    if len(sys.argv) < 2: 
        print_help(sys.argv[0])
        sys.exit(1) 

    # input the command line arguments
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            print_help(sys.argv[0])
            sys.exit()
        elif o in ("-l", "--list"):
            biaslist = a
        elif o in ("-o", "--output"):
            biasfits = a
        elif o in ("-c", "--combine"):
            combinetype = a
        elif o in ("-r", "--reject"):
            rejecttype = a
        else:
            assert False, "unhandled option"
    
    # Check/Output the biaslist 
    if biaslist == None:
        print sys.argv[0],"ERROR: No input list given"
        sys.exit(1)
    else:
        atlist = "\"@" + biaslist 
        if(verbose == True):
            print sys.argv[0],"STATUS: Using bias list of",biaslist

    # Check/Output the biasfits  
    if biasfits == None:
        print sys.argv[0],"ERROR: No output fits filename given"
        sys.exit(1)
    else:
        if(verbose == True):
            print sys.argv[0],"STATUS: Output master bias image as",biasfits
           
    # Edit the images headers by changing IMAGETYP to "zero"
    if(verbose == True):
        print sys.argv[0],"STATUS: Running hedit to set IMAGETYP=zero"
    if(verbose == True):
        iraf.hedit(atlist, fields="IMAGETYP", value="zero", add="yes", verify="no", update="yes", show="yes")
    else:
        iraf.hedit(atlist, fields="IMAGETYP", value="zero", add="yes", verify="no", update="yes", show="no")

    # Run the zerocombine
    if(verbose == True):
        print sys.argv[0],"STATUS: Running zerocombine on imagelist =",biaslist,"and output as",biasfits,"with combine =",combinetype,"and reject =",rejecttype
    # hardwired some of the parameters at the moment
    iraf.zerocombine(atlist, output=biasfits, combine=combinetype, ccdtype="zero", reject=rejecttype, scale="none", nlow=0, nhigh=1, nkeep=1, mclip="yes", lsigma=3.0, hsigma=3.0, rdnoise="0.", gain ="1.")
    
    # Rename the logfile 
    outlog=biasfits.rstrip('fits')+'log'
    os.rename("logfile",outlog)

    # Finally print the ending message
    print sys.argv[0],"STATUS: zerocombine done"

if __name__ == "__main__":
    main()


#----------------------------------------------------
#  Date      Version    Note
#----------------------------------------------------
# 15/02/2012  1.0       Initial version
    
