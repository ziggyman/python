import numpy as np
import pyfits
import pyraf
from pyraf import iraf
from sklearn.decomposition import PCA
#import zap

cubeMinusSky = '/Volumes/obiwan/azuri/spectra/sparsepak/2014B-0176/kwb_141015_080034_ori_otxzd.ms.fits'
sky = '/Volumes/obiwan/azuri/spectra/sparsepak/2014B-0176/skykwb_141015_080034_ori_otxzd.fits'

cubePlusSky = '/Volumes/obiwan/azuri/spectra/sparsepak/2014B-0176/kwb_141015_080034_ori_otxzd_ms+sky.fits'
cubeOut = '/Volumes/obiwan/azuri/spectra/sparsepak/2014B-0176/kwb_141015_080034_ori_otxzd_ms+sky_zap.fits'

hdulistCubeMinusSky = pyfits.open(cubeMinusSky)
hdulistSky = pyfits.open(sky)

dataCubeMinusSky = hdulistCubeMinusSky[0].data
dataSky = hdulistSky[0].data

print("dataCubeMinusSky.shape = ",dataCubeMinusSky.shape)
print("dataSky.shape = ",dataSky.shape)

dataCubePlusSky = np.ndarray(shape=dataCubeMinusSky.shape, dtype=type(dataSky[0]))
for fiber in range(dataCubeMinusSky.shape[0]):
    dataCubePlusSky[fiber,:] = dataCubeMinusSky[fiber,:] + dataSky[:]

hdulistCubeMinusSky[0].data = dataCubePlusSky
hdulistCubeMinusSky.writeto(cubePlusSky, clobber=True)

# create the PCA instance
pca = PCA(2)
# fit on data
pca.fit(dataCubePlusSky)
# access values and vectors
print('pca.components_ = ',pca.components_)
print('pca_explained_variance_ = ',pca.explained_variance_)
# transform data
B = pca.transform(dataCubePlusSky)
print('B = ',B)


#zap.process(cubePlusSky, outcubefits=cubeOut)