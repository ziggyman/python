import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from PyAstronomy import pyasl
import astropy
from astropy import units as u
from astropy.nddata import CCDData
import numpy as np

colours = ['k',
           'grey',
           'lightcoral',
           'maroon',
           'r',
           'orangered',
           'peru',
           'gold',
           'olive',
           'yellow',
           'lawngreen',
           'g',
           'springgreen',
           'c',
           'deepskyblue',
           'navy',
           'blueviolet',
           'deeppink']

with open('/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/sensFuncs/fluxstds.list','r') as f:
    fluxstd_names = f.readlines()

for i in range(int(len(fluxstd_names)/2)):
    sensfunc_b = astropy.io.ascii.read(fluxstd_names[2*i].strip())
    sensfunc_r = astropy.io.ascii.read(fluxstd_names[2*i+1].strip())
    ss_b = np.argsort(sensfunc_b['wave'])
    ss_r = np.argsort(sensfunc_r['wave'])
    # interpolate the sensfunc onto the observed wavelength axis
    plt.plot(sensfunc_b['wave'][ss_b], sensfunc_b['S'][ss_b]/np.max(sensfunc_r['S'][ss_r]), 'r' if 'bs01153' in fluxstd_names[i] else colours[i], label=fluxstd_names[i][fluxstd_names[2*i].rfind('/')+1:])
    plt.plot(sensfunc_r['wave'][ss_r], sensfunc_r['S'][ss_r] / np.max(sensfunc_r['S'][ss_r]), 'r' if 'bs01153' in fluxstd_names[i] else colours[i])#label=fluxstd_names[i][fluxstd_names[2*i].rfind('/')+1:])
plt.legend()
plt.show()


for i in range(int(len(fluxstd_names)/2)):
    sensfunc_b = astropy.io.ascii.read(fluxstd_names[2*i].strip())
    sensfunc_r = astropy.io.ascii.read(fluxstd_names[2*i+1].strip())
    ss_b = np.argsort(sensfunc_b['wave'])
    ss_r = np.argsort(sensfunc_r['wave'])
    idx_b = np.where(np.abs(sensfunc_b['wave'][ss_b] - 4861.) < 20)[0]
    print('idx_b = ',idx_b)
    idx_r = np.where(np.abs(sensfunc_r['wave'][ss_r] - 6563.) < 20)[0]
    print('idx_r = ',idx_r)
    point = np.mean(sensfunc_r['S'][ss_r][idx_r]) / np.mean(sensfunc_b['S'][ss_b][idx_b])
    print('i = ',type(i),': ',i,', point = ',type(point),': ',point)
    plt.scatter(float(i),float(point),c=('r' if 'bs01153' in fluxstd_names[i] else colours[i]),label=fluxstd_names[i][fluxstd_names[2*i].rfind('/')+1:])
    # interpolate the sensfunc onto the observed wavelength axis
    #plt.plot(sensfunc_b['wave'][ss_b], sensfunc_b['S'][ss_b], colours[i], label=fluxstd_names[i][fluxstd_names[2*i].rfind('/')+1:])
    #plt.plot(sensfunc_r['wave'][ss_r], sensfunc_r['S'][ss_r] / np.mean(sensfunc_b['S'][ss_b]), 'r' if 'bs01153' in fluxstd_names[i] else colours[i])#label=fluxstd_names[i][fluxstd_names[2*i].rfind('/')+1:])
plt.legend()
plt.show()
