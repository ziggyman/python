import numpy as np
from matplotlib  import pyplot as plt

def readGvaramadzeFile():
    lbt_filename = '/Users/azuri/daten/uni/HKU/Pa30/variability/obj_spect_PA256.dat'

    with open(lbt_filename,'r') as f:
        lines = f.readlines()
    wLen = []
    flux = []
    for iLine in range(len(lines)):
        lines[iLine] = lines[iLine].strip().replace('    ',' ').split(' ')
        print('lines[',iLine,'] = ',lines[iLine])
        wLen.append(float(lines[iLine][0]))
        flux.append(float(lines[iLine][1]))
    print('wLen = ',wLen)
    print('flux = ',flux)
    plt.plot(wLen,flux)
    plt.show()
    return [wLen,flux]

def readLBTFiles():
    fNames = ['/Users/azuri/daten/uni/HKU/Pa30/variability/ritter/b2a.dat',
              '/Users/azuri/daten/uni/HKU/Pa30/variability/ritter/b2b.dat',
              '/Users/azuri/daten/uni/HKU/Pa30/variability/ritter/b2c.dat',
              '/Users/azuri/daten/uni/HKU/Pa30/variability/ritter/b2d.dat',
              '/Users/azuri/daten/uni/HKU/Pa30/variability/ritter/b2e.dat',
             ]
    wavelengths = []
    fluxes = []
    for fName in fNames:
        with open(fName,'r') as f:
            lines = f.readlines()
        wLen = []
        flux = []
        for iLine in np.arange(1,len(lines),1):
            lines[iLine] = lines[iLine].strip().replace('  ',' ').split(' ')
            print('lines[',iLine,'] = ',lines[iLine])
            wLen.append(float(lines[iLine][0]))
            flux.append(float(lines[iLine][1]))
        wavelengths.append(wLen)
        fluxes.append(flux)
        plt.plot(wLen,flux)
    
    wLen = []
    flux = []
    for i in range(len(wavelengths[0])):
        wLens = []
        flxs = []
        for j in range(5):
            wLens.append(wavelengths[j][i])
            flxs.append(fluxes[j][i])
        wLen.append(np.mean(wLens))
        print('wLens[',i,'] = ',wLens,': mean = ',wLen[i])
        flux.append(np.mean(flxs))
        print('flxs[',i,'] = ',flxs,': mean = ',flux[i])
    plt.plot(wLen,flux)
    plt.show()

    return [wLen,flux]

if __name__ == '__main__':
    readGvaramadzeFile()
    readLBTFiles()