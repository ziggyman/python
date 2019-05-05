from drUtils import combine, separateFileList, silentRemove
import numpy as np
from shutil import copyfile
from myUtils import readFileToArr

def test_separateFileList():
    inList = '/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/allFits.list'
    suffixes = ['','z','zf']

    copyfile(inList, inList+'bak')
    silentRemove(inList[:inList.rfind('/')+1]+'*.list')
    copyfile(inList+'.bak', inList)
    separateFileList(inList, suffixes, True)

def test_combine():
    #combine(ccdImages,
    #        combinerMethod='median',
    #        clippingMethod='None',
    #        clippingParameters=None,
    #        scaling=False,
    #        fitsOutName=None)
    inFiles = readFileToArr('/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/bias.list')
    combinedImage = combine(inFiles, fitsOutName='/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/combinedBias.fits')
    print('combinedImage = ',type(combinedImage))
    print('dir(combinedImage) = ',dir(combinedImage))
    print('combinedImage.shape = ',combinedImage.shape)
    print('combinedImage.ndim = ',combinedImage.ndim)
    print('combinedImage.size = ',combinedImage.size)
    print('mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='minmax',
                    clippingParameters={'min_clip':0.,'max_clip':700.},
                    scaling=True,
                    fitsOutName='/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average minmax: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='extrema',
                    clippingParameters={'nlow':2,'nhigh':2},
                    scaling=True,
                    fitsOutName='/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average extrema: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':0,'low_thresh':-3.,'high_thresh':3.,'func':np.ma.median},
                    scaling=True,
                    fitsOutName='/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average sigma 0: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':1,'low_thresh':-3.,'high_thresh':3.,'func':np.ma.median},
                    scaling=True,
                    fitsOutName='/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average sigma 1: mean(combinedImage) = ',np.mean(combinedImage))

    combinedImage = combine(inFiles,
                    combinerMethod='average',
                    clippingMethod='sigma',
                    clippingParameters={'niter':3,'low_thresh':-3.,'high_thresh':3.,'func':np.ma.median},
                    scaling=True,
                    fitsOutName='/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/combinedBias_av.fits')
    print('average sigma 3: mean(combinedImage) = ',np.mean(combinedImage))

#test_separateFileList()
test_combine()
