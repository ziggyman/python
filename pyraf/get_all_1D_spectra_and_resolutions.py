import os
from astropy.io import fits
from shutil import copyfile

mainDir = '/Volumes/Untitled/mashtun/SOURCES/FULL_DATASETS'
outDir = '/Volumes/discovery/azuri/data/hash_spectra'

subDirs = [x[0] for x in os.walk(mainDir)]
print('subDirs = ',subDirs)

spectra = []
for subDir in subDirs:
    for root, dirs, files in os.walk(subDir):
        for file in files:
            if file.endswith(".fits"):
                print(os.path.join(root, file))
                spectra.append(os.path.join(root, file))
print('spectra = ',spectra)

spectra1D = []
i = 0
for spectrum in spectra:
    i += 1
    try:
        hdul = fits.open(spectrum)
        header = hdul[0].header
        #data, header = fits.getdata(spectrum, header=True)
#        print('data = ',data)
#        print('type(data) = ',type(data))
#        print('dir(data) = ',dir(data))
#        print('data.shape = ',data.shape)
#        print('len(data.shape) = ',len(data.shape))
        if header['NAXIS'] == 1:
            if 'CDELT1' in header:
                print('spectrum = ',spectrum,': CDELT1 = ',header['CDELT1'])
                if header['CDELT1'] < 1.6:
                    print(len(spectra1D),': ',spectrum,' is good, ',len(spectra)-i,' to go')
                    copyfile(spectrum,os.path.join(outDir,spectrum[spectrum.rfind('/')+1:]))
                    spectra1D.append(spectrum)

    except:
        pass
#print('spectra1D = ',len(spectra1D),': ',spectra1D)
