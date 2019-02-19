import numpy as np
import pyfits
import pyraf
from pyraf import iraf
from sklearn.decomposition import PCA
execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate

standardStar = True
cubePlusSky = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum.fits'

cubeMinusSkyOut = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum-skyMean.fits'
skyOut = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum_sky.fits'

componentsOutRoot = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Pa30_WIYN2014-10-15_botzfxsEcBld_sum_PCA'
if standardStar:
    cubePlusSky = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Hiltner102_botzfxsEcBld_combined.fits'

    cubeMinusSkyOut = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Hiltner102_botzfxsEcBld_combined-skyMean.fits'
    skyOut = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Hiltner102_botzfxsEcBld_combined_sky.fits'

    componentsOutRoot = '/Volumes/obiwan/azuri/spectra/sparsepak/stella/Hiltner102_botzfxsEcBld_combined_PCA'



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

iFiber = 51

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
hdulistCubeMinusSkyOut.writeto(skyOut, clobber=True)

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

hdulistCubeMinusSkyOut.data = dataCubePlusSky
for i in range(dataCubePlusSky.shape[0]):
    hdulistCubeMinusSkyOut.data[i,:] = dataCubePlusSky[i,:] - pca.mean_
hdulistCubeMinusSkyOut.writeto(cubeMinusSkyOut, clobber=True)
#zap.process(cubePlusSky, outcubefits=cubeOut)