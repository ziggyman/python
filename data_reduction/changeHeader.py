import astropy.io.fits as pyfits
import os
from drUtils import getHeaderValue

path = '/Users/azuri/spectra/saao/saao_mar2014'
subPaths = ['night1','night2','night3','night4','night5','night6','night7',]

def addEXPTYPE():
    for subPath in subPaths:
        p = os.path.join(path,subPath,'allFits.list')
        with open(p,'r') as f:
            lines = f.readlines()
        lines = [os.path.join(path,subPath,line.strip()) for line in lines]
        for inputFileName in lines:
            with pyfits.open(inputFileName, mode='update') as hdulist:
                # Change something in hdul.
                print('writeFits: len(hdulist) = ',len(hdulist),', inputFileName = ',inputFileName)
                print('header = ',hdulist[len(hdulist)-1].header)
                newType = hdulist[len(hdulist)-1].header['IMAGETYP']
                print('newType = <'+newType+'>')
                if newType == 'COMPARISON':
                    newType = 'ARC'
                elif newType == 'zero':
                    newType = 'BIAS'
                elif newType == 'flat':
                    newType = 'FLAT'
                elif newType == 'object':
                    newType = 'SCIENCE'
                hdulist[len(hdulist)-1].header['EXPTYPE'] = newType
                hdulist.flush()  # changes are written back to original.fits

if __name__ == '__main__':
    addEXPTYPE()
