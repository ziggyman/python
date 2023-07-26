import numpy as np
from drUtils import getImageData
from matplotlib import pyplot as plt

fName = '/Users/azuri/Downloads/O59C01010_X1D.fits'
"""[('SPORDER', '>i2'), ('NELEM', '>i2'), ('WAVELENGTH', '>f8', (1024,)), ('GROSS', '>f4', (1024,)), ('BACKGROUND', '>f4', (1024,)), ('NET', '>f4', (1024,)), ('FLUX', '>f4', (1024,)), ('ERROR', '>f4', (1024,)), ('NET_ERROR', '>f4', (1024,)), ('DQ', '>i2', (1024,)), ('A2CENTER', '>f4'), ('EXTRSIZE', '>f4'), ('MAXSRCH', '>i2'), ('BK1SIZE', '>f4'), ('BK2SIZE', '>f4'), ('BK1OFFST', '>f4'), ('BK2OFFST', '>f4'), ('EXTRLOCY', '>f4', (1024,)), ('OFFSET', '>f4')]"""

imgData = getImageData(fName,1)
wLen = []
flux = []
for i in np.arange(len(imgData)-1,-1,-1):
    wave = imgData[i][2]
    print('i = ',i,': wave = ',len(wave),': ',wave)
    for w in wave:
        wLen.append(w)
    flx = imgData[i][6]
    print('i = ',i,': flx = ',len(flx),': ',flx)
    for f in flx:
        flux.append(f)
#    plt.plot(wave,flx)
#    plt.show()

print('wLen = ',len(wLen),': ',wLen)
print('flux = ',len(flux),': ',flux)
plt.plot(wLen, flux)
plt.show()
