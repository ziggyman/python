import numpy as np
#import pyfits
#import pyraf
import matplotlib.pyplot as plt

from drUtils import getImageData

FName = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_av_x_wl_flt_cal_mSky_obs_not_smoothed_minComb.fits'
ObjectArea = [1123,1866]
imageData = getImageData(FName,0)

print('imageData = ',imageData.shape,': ',imageData)
wLen = np.arange(3347.6939,7695.65,2.12)
print('wLen = ',wLen.shape,': ',wLen)

#erg/s -> W
imageData = imageData * 0.0000001
print('imageData in W / cm^2 / A = ',imageData.shape,': ',imageData)

##integral over dLambda
#imageData = np.array([imageData[:,i] * 2.12 for i in range(imageData.shape[1])])
#print('imageData in W / cm^2 = ',imageData.shape,': ',imageData)

#calculate per arcsec^2 (0.8" slit width, pixel size 0.25")
imageData = np.array([imageData[:,i] / 0.2 for i in range(imageData.shape[1])])
print('imageData in W / cm^2 / A / arcsec^2 = ',imageData.shape,': ',imageData)

#cm^-2 -> m^-2
imageData = np.array([imageData[:,i] * 10000. for i in range(imageData.shape[1])])
print('imageData in W / m^2 / A / arcsec^2 = ',imageData.shape,': ',imageData)

plt.imshow(imageData,vmin=0)#,vmax=1e-17)
plt.show()

imageData = imageData[1554:1630,1144:1880]
print('imageData interesting = ',imageData.shape,': ',imageData)

#subtract background
imageData = np.array([imageData[:,i] - np.median(imageData[:,i]) for i in range(imageData.shape[1])])
print('imageData - background = ',imageData.shape,': ',imageData)

#remove star residual
imageData[np.where(imageData > 2.e-21)] = 0.

fig, ax = plt.subplots()
h = ax.imshow(imageData,vmin=0,vmax=1.5e-21)
ax.set_xlabel('pixel')
ax.set_ylabel('pixel')
cbar = fig.colorbar(h)
cbar.set_label('$B_\lambda\ [\mathrm{W}\ \mathrm{m}^{-2}\ \mathrm{\AA}^{-1}\ \mathrm{arcsec}^{-2}]$')
plt.show()
fig.savefig(FName[:FName.rfind('.')]+'_surfaceBrightness_lambda.png', format='png', bbox_inches='tight', pad_inches=0.1)
