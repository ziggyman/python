import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import os

from drUtils import getWavelengthArr
import csvFree,csvData

path = '/Users/azuri/daten/uni/HKU/Pa30/CAHA_3.5m/'

def plot_Pa30_NE_intspec():
    # all nebula spectra together
    fitsName = os.path.join(path,'Pa30-NE_intspec.fits')

    hdulist = pyfits.open(fitsName)
    print('len(hdulist) = ',len(hdulist))
    print('hdulist = ',hdulist)
    print('hdulist[0].header = ',hdulist[0].header)
    print('hdulist[0].data = ',hdulist[0].data)
    print('hdulist[1].header = ',hdulist[1].header)
    print('hdulist[1].data = ',hdulist[1].data)
    print('hdulist[2].header = ',hdulist[2].header)
    print('hdulist[2].data = ',hdulist[2].data)
    print('hdulist[3].header = ',hdulist[3].header)
    print('hdulist[3].data = ',hdulist[3].data)

    wave = hdulist[1].data
    spec = hdulist[0].data

    print('wave = ',type(wave),': ',len(wave))
    print('spec = ',type(spec),': ',len(spec))

    plt.plot(wave,spec)
    plt.title('average spectrum Pa30-NE')
    plt.show()

# Pa 30 CS
def plot_Pa30_CS():
    fitsName = os.path.join(path,'Pa30_CSPN.fits')

    hdulist = pyfits.open(fitsName)
    print('len(hdulist) = ',len(hdulist))
    print('hdulist = ',hdulist)
    print('hdulist[0].header = ',hdulist[0].header)
    print('hdulist[0].data = ',hdulist[0].data)

    wave = getWavelengthArr(fitsName,0)
    spec = hdulist[0].data

    print('wave = ',type(wave),': ',len(wave))
    print('spec = ',type(spec),': ',len(spec))

    plt.plot(wave,spec)
    plt.title('Pa 30 CS')
    plt.show()

def reduce_Pa30():
    logFileName = os.path.join(path,'231016_PMAS','night[Tel35m-23B-11].log')
    logcvs = csvFree.readCSVFile(logFileName,'||')
    print('logcvs.size() = ',logcvs.size())

if __name__ == '__main__':
    plot_Pa30_CS()
    STOP
    reduce_Pa30()
