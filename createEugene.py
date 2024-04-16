pythonFileIn = '/Users/azuri/entwicklung/python/data_reduction/drUtils.py'
pythonFileOut = '/Users/azuri/entwicklung/python/data_reduction/drUtilsEugene.py'

copyFunctions = ['getImageData','getWavelengthArr','writeFits1D','normalizeX','sfit','sigmaReject','gauss','gauss_lin','sigmaRej','getHeaderValue','curve_fit','continuum','getHeader','fitLines']

with open(pythonFileIn,'r') as f:
    lines = f.readlines()

functionName = ''
with open(pythonFileOut,'w') as f:
    f.write('import numpy as np\n')
    f.write('import matplotlib.pyplot as plt\n')
    f.write('import astropy.io.fits as pyfits\n')
    f.write('from PyAstronomy import pyasl\n')
    f.write('from scipy.optimize import curve_fit\n')
    f.write('from scipy import exp\n')
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('\n')
    for line in lines:
        if line[0:4] == 'def ':
            functionName = line[4:line.rfind('(')]
            print('functionName = ',functionName)
        if functionName in copyFunctions:
            f.write(line)
