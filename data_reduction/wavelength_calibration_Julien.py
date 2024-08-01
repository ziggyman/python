from astropy.coordinates import EarthLocation
from drUtils_Julien import dispCor
import os
import numpy as np

workPath = '/Users/azuri/daten/uni/HKU/interns_projects/julien/PHR1750-1803'
ecFiles = ...#copy each row containing PN spectrum to individual file, do "ls .... > ecFiles"
ecdFiles = ...#copy ecFiles to ecdFiles, replace Ec.fits with Ecd.fits
observatoryLocation = EarthLocation.of_site('SAAO')


def readFileToArr(fname):
    with open(fname, "r") as f:
        lines = f.readlines()

    """remove empty lines"""
    lines = list(filter(len, lines))

    linesOut = [line.strip() for line in lines]
    return linesOut


def getListOfFiles(fname):
    fList = readFileToArr(fname)
    print('fname = ',fname,': fList = ',fList)
    if len(fList) > 0:
        if fList[0].rfind('/') == -1:
            fList = [os.path.join(workPath, fileName) for fileName in fList]
    return fList


inputList = getListOfFiles(os.path.join(workPath,'arcList_otzxfiEc.list'))
wavelengthsOrig = []
for i in range(len(inputList)):
    wLenStr = readFileToArr(inputList[i][:inputList[i].rfind('.')]+'_wLenOrig.dat')
    wLens = [float(wLen) for wLen in wLenStr]
    wavelengthsOrig.append(np.asarray(wLens))

dispCor(ecFiles,
        getListOfFiles(os.path.join(workPath,'arcList_otzxfiEc.list')),
        wavelengthsOrig,
        ecdFiles,
        observatoryLocation,
        'RA',#RA',#
        'DEC',#DEC',#
        'DATE-OBS',
        doHelioCor = True)
