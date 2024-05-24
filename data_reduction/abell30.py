import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
import os

from scipy.optimize import curve_fit

def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

path = '/Users/azuri/daten/uni/HKU/Kamila/OIII_ratios/Red and blue cube/'
fitsName = os.path.join(path,'Blue_final_cube_angstroms.fits')

cs_spectrum_1977_name = '/Users/azuri/daten/uni/HKU/Kamila/Abell30_spec1977.txt'
with open(cs_spectrum_1977_name,'r') as fn:
    lines = fn.readlines()
wave_1977 = []
spec_1977 = []
for line in lines:
    wave_1977.append(float(line[:line.find(' ')]))
    spec_1977.append(float(line[line.find(' ')+1:line.rfind('\n')]))
spec_1977 = spec_1977 / np.mean(np.array(spec_1977))

cs_spectrum_1964_name = '/Users/azuri/daten/uni/HKU/Kamila/Abell30_spec_im.txt'
with open(cs_spectrum_1964_name,'r') as fn:
    lines = fn.readlines()
wave_1964 = []
spec_1964 = []
for line in lines:
    wave_1964.append(float(line[:line.find(' ')]))
    spec_1964.append(float(line[line.find(' ')+1:line.rfind('\n')]))
#print('wave_1964 = ',wave_1964)
#idx_1964 = np.where(np.array(wave_1964) > 4342.)[0]
#print('idx_1964 = ',idx_1964)
#idx_1977 = np.where(np.array(wave_1977) < 5071.)[0]
#print('idx_1977 = ',idx_1977)
spec_1964 = np.array(spec_1964)  / np.mean(np.array(spec_1964))# * np.mean(np.array(spec_1964)[idx_1964])
plt.plot(wave_1964,spec_1964)
plt.plot(wave_1977,spec_1977)
plt.show()

cs_spectrum_name = '/Users/azuri/daten/uni/HKU/Kamila/Abell30_CS_centre_sum_spec_new.fits'

pa30_cs_spectrum_name = '/Users/azuri/daten/uni/HKU/Pa30/Pa30_GT080716_cal_sum_cleaned_scaled.fits'

def getImageData(fname,hduNum=1):
    hdulist = pyfits.open(fname)
    scidata = hdulist[hduNum].data
    hdulist.close()
    return np.array(scidata)

def getWavelengthArr(fname,hduNum=1,dim='3'):
    hdulist = pyfits.open(fname)
    header = hdulist[hduNum].header
    hdulist.close()
    if 'CDELT'+dim in header.keys():
        cdelt = header['CDELT'+dim]
        wLen = ((np.arange(header['NAXIS'+dim]) + 1.0) - header['CRPIX'+dim]) * cdelt + header['CRVAL'+dim]
    elif 'CD1_'+dim in header.keys():
        cdelt = header['CD1_'+dim]
        wLen = ((np.arange(header['NAXIS'+dim]) + 1.0) - header['CRPIX'+dim]) * cdelt + header['CRVAL'+dim]
    else:
        print('WARNING: neither CDELT'+dim+' nor CD1_'+dim+' found in header of file <'+fname+'>')
        wLen = np.arange(header['NAXIS'+dim]) + 1.0
    return wLen

def getHeaderValue(fname, keyword, hduNum=0):
    hdulist = pyfits.open(fname)
    header = hdulist[hduNum].header
    hdulist.close()
    try:
        return header[keyword]
    except:
        return None

# CS Spectrum
cs_wave = getWavelengthArr(cs_spectrum_name,hduNum=0,dim='1')
cs_spec = getImageData(cs_spectrum_name,0)
cs_spec = cs_spec / np.mean(cs_spec)
plt.plot(cs_wave,cs_spec)

pa30_wave = getWavelengthArr(pa30_cs_spectrum_name,0,dim='1')
pa30_spec = getImageData(pa30_cs_spectrum_name,0)
#idx = np.where((pa30_wave > 3763.) & (pa30_wave < 5235.))[0]
pa30_spec = pa30_spec / np.mean(pa30_spec)# * np.mean(cs_spec)
plt.plot(pa30_wave,pa30_spec)
plt.plot(wave_1964,spec_1964)
plt.plot(wave_1977,spec_1977)
plt.show()
STOP

hdulist = pyfits.open(fitsName)
print('len(hdulist) = ',len(hdulist))
#print('hdulist = ',hdulist)
print('hdulist[0].header = ',hdulist[0].header)
#print('hdulist[0].data = ',hdulist[0].data)
print('hdulist[1].header = ',hdulist[1].header)
#print('hdulist[1].data = ',hdulist[1].data)

wLen = getWavelengthArr(fitsName)
print('wLen = ',wLen.shape,': ',wLen)

regionOfInterest = np.where((wLen > 4948.) & (wLen < 5030.))[0]
regionOfInterest4959 = np.where((wLen > 4950.) & (wLen < 4962.))[0]
regionOfInterest5007 = np.where((wLen > 4998.) & (wLen < 5010.))[0]

idxForContinuum = np.where(np.logical_or(np.logical_or(((wLen > 4940.) & (wLen < 4950.)),((wLen > 4962.) & (wLen < 4998.))),((wLen > 5020.) & (wLen < 5030.))))[0]
wLenForContinuum = wLen[idxForContinuum]

data = hdulist[1].data
print('data = ',data.shape,': ',data)

oIII4959 = np.zeros([data.shape[1],data.shape[2]])
oIII5007 = np.zeros([data.shape[1],data.shape[2]])
integratedSpec = np.zeros([data.shape[1],data.shape[2]])
intensities = np.zeros([data.shape[1],data.shape[2]])

for i in range(data.shape[1]):
    for j in range(data.shape[2]):
        spec = data[:,i,j]
        integratedSpec[i,j] = np.sum(spec)
#        print('i = ',i,', j = ',j,': spec = ',spec)
#        print('idxForContinuum = ',idxForContinuum)
        specForContinuum = spec[idxForContinuum]
        if not np.isnan(specForContinuum[0]):
#            print('i = ',i,', j = ',j,': wLenForContinuum = ',wLenForContinuum)
#            print('i = ',i,', j = ',j,': specForContinuum = ',specForContinuum)
            popt, pcov = curve_fit(f, wLenForContinuum, specForContinuum) # your data x, y to fit
            spec[regionOfInterest] = spec[regionOfInterest] - f(wLen[regionOfInterest],popt[0],popt[1])
            plt.plot(wLen[regionOfInterest],spec[regionOfInterest])
            oIII4959[i,j] = trapz(spec[regionOfInterest4959],dx=getHeaderValue(fitsName,'CDELT3',1))
            oIII5007[i,j] = trapz(spec[regionOfInterest5007],dx=getHeaderValue(fitsName,'CDELT3',1))
            intensities[i,j] = np.max(spec[regionOfInterest5007])
#        print(wLenForContinuum = ',wLenForContinuum)
#        STOP

plt.show()

print('oIII4959 = ',oIII4959)
print('oIII5007 = ',oIII5007)

ratio = oIII5007/oIII4959
print('ratio = ',ratio)

fig, (ax1, ax2) = plt.subplots(figsize=(13, 6), ncols=2)
cb = ax1.imshow(ratio, vmin=2., vmax=4.)
fig.colorbar(cb,ax=ax1)
integ = ax2.imshow(integratedSpec, vmin=1.5e3, vmax=2.5e3)
fig.colorbar(integ,ax=ax2)
plt.show()

hist,limits = np.histogram(ratio,bins=20,range=[2.,4.])
print('hist = ',len(hist),': ',hist)
print('limits = ',len(limits),': ',limits)
plt.hist(limits[:-1], limits, weights=hist)
plt.show()

plt.scatter(np.ndarray.flatten(intensities),np.ndarray.flatten(ratio))
plt.xlabel('intensity [OIII] 5007')
plt.ylabel('ratio [OIII] 5007/4959')
plt.show()
