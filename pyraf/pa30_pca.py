import numpy as np
import pyfits
import pyraf
from pyraf import iraf
from sklearn.decomposition import PCA
from myUtils import getWavelength
import matplotlib.pyplot as plt

standardStar = False
cubePlusSky = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/pa30_zdtsEcndr.fits'

cubeMinusSkyOut = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/pa30_zdtsEcndr-skyMean.fits'
cubeMinusSkyOutMedian = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/pa30_zdtsEcndr-skyMedian.fits'
skyOutMean = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/pa30_zdtsEcndr_skyMean.fits'
skyOutMedian = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/pa30_zdtsEcndr_skyMedian.fits'
componentsOutRoot = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/pa30_zdtsEcndr_PCA'

if standardStar:
    cubePlusSky = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/kwb_141015_114545_ori_zdtsEcn1dr.fits'

    cubeMinusSkyOut = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/kwb_141015_114545_ori_zdtsEcn1dr-skyMean.fits'
    cubeMinusSkyOutMedian = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/kwb_141015_114545_ori_zdtsEcn1dr-skyMedian.fits'
    skyOutMean = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/kwb_141015_114545_ori_zdtsEcn1dr_skyMean.fits'
    skyOutMedian = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/kwb_141015_114545_ori_zdtsEcn1dr_skyMedian.fits'
    componentsOutRoot = '/Users/azuri/daten/uni/HKU/Pa30/sparsepak/spectra/kwb_141015_114545_ori_zdtsEcn1dr_PCA'



hdulistCubePlusSky = pyfits.open(cubePlusSky)
wavelength = getWavelength(hdulistCubePlusSky[0].header, 1)
print("wavelength = ",wavelength)
dataCubePlusSky = hdulistCubePlusSky[0].data

# create the PCA instance
nComponents = 10
pca = PCA(nComponents)
# fit on data
pca.fit(dataCubePlusSky)
# access values and vectors
print('pca.components_ = ',pca.components_.shape,': ',pca.components_)
print('pca_explained_variance_ = ',pca.explained_variance_.shape,': ',pca.explained_variance_)
# transform data
B = pca.transform(dataCubePlusSky)
print('B = ',B.shape,': ',B)

iFiber = 31

componentsSum = pca.components_[0,:] * B[iFiber,0]

plt.plot(wavelength, dataCubePlusSky[iFiber,:], 'g-', label = str(iFiber))
plt.plot(wavelength, pca.mean_, 'c-', label = 'mean')
for i in range(nComponents):
    plt.plot(wavelength, pca.components_[i,:], label = 'PCA component '+str(i))
#    plt.plot(wavelength, (pca.components_[i,:] * B[iFiber,i]) + pca.mean_, label = 'PCA component '+str(i))
    if i > 0:
        componentsSum += pca.components_[i,:] * B[iFiber,i]
plt.plot(wavelength, componentsSum + pca.mean_, 'b-', label = 'sum')

inverseTransform = pca.inverse_transform(B)
plt.plot(wavelength, inverseTransform[iFiber,:], 'r-', label = 'fiber1 inverse transform')

plt.legend()
plt.show()

hdulistCubeMinusSkyOut = hdulistCubePlusSky
hdulistCubeMinusSkyOut.data = pca.mean_
hdulistCubeMinusSkyOut.writeto(componentsOutRoot+'mean.fits', clobber=True)
hdulistCubeMinusSkyOut.writeto(skyOutMean, clobber=True)

plt.plot(wavelength, pca.mean_, label = 'PCA mean')
#plt.xlim(5690., 6760.)
plt.legend()
plt.show()

for i in range(nComponents):
    hdulistCubeMinusSkyOut.data = pca.components_[i,:]
    hdulistCubeMinusSkyOut.writeto(componentsOutRoot+str(i)+'.fits', clobber=True)
    plt.plot(wavelength, pca.components_[i,:], label = 'PCA component '+str(i))
    plt.legend()
    plt.show()

# subtrace mean as sky
hdulistCubeMinusSkyOut.data = dataCubePlusSky
for i in np.arange(0,dataCubePlusSky.shape[0],1):
    hdulistCubeMinusSkyOut.data[i,:] = dataCubePlusSky[i,:] - pca.mean_
hdulistCubeMinusSkyOut.writeto(cubeMinusSkyOut, clobber=True)

# subtrace mean as sky
hdulistCubeMinusSkyOut.data = dataCubePlusSky
skyArr = np.zeros(dataCubePlusSky.shape,dtype=np.float32)
for i in np.arange(0,dataCubePlusSky.shape[1],1):
    col = dataCubePlusSky[:,i]
    sky = np.median(dataCubePlusSky[:,i])
    skyArr[:,i] = sky
    print('column ',i,': col = ',col.shape,': ',col,', sky = ',sky.shape,': ',sky)
    hdulistCubeMinusSkyOut.data[:,i] = col - sky
hdulistCubeMinusSkyOut.writeto(cubeMinusSkyOutMedian, clobber=True)

hdulistCubeMinusSkyOut.data = skyArr
hdulistCubeMinusSkyOut.writeto(skyOutMedian, clobber=True)

#zap.process(cubePlusSky, outcubefits=cubeOut)