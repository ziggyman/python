import os

import pyfits
import gzip
import numpy as np

"""
For the description of the datamodel see
https://gaia.esac.esa.int/documentation/GDR1/datamodel/Ch1/gaia_source.html
"""

dataDir = '/Volumes/external/azuri/data/gaia/cdn.gea.esac.esa.int/Gaia/gaia_source/fits/'
fileNameRoot = 'GaiaSource_%03d-%03d-%03d.fits'#% (release, tract, patch)

solution_id = 1
source_id  = 2
random_index = 3
ref_epoch = 4
ra = 5
ra_error = 6
dec = 7
dec_error = 8
parallax = 9
parallax_error = 10
pmra = 11
pmra_error = 12
pmdec = 13
pmdec_error = 14
ra_dec_corr = 15
ra_parallax_corr = 16
ra_pmra_corr = 17
ra_pmdec_corr = 18
dec_parallax_corr = 19
dec_pmra_corr = 20
dec_pmdec_corr = 21
parallax_pmra_corr = 22
parallax_pmdec_corr = 23
pmra_pmdec_corr = 24
astrometric_n_obs_al = 25
astrometric_n_obs_ac = 26
astrometric_n_good_obs_al = 27
astrometric_n_good_obs_ac = 28
astrometric_n_bad_obs_al = 29
astrometric_n_bad_obs_ac = 30
astrometric_delta_q = 31
astrometric_excess_noise = 32
astrometric_excess_noise_sig = 33
astrometric_primary_flag = 34
astrometric_relegation_factor = 35
astrometric_weight_al = 36
astrometric_weight_ac = 37
astrometric_priors_used = 38
matched_observations = 39
duplicated_source = 40
scan_direction_strength_k1 = 41
scan_direction_strength_k2 = 42
scan_direction_strength_k3 = 43
scan_direction_strength_k4 = 44
scan_direction_mean_k1 = 45
scan_direction_mean_k2 = 46
scan_direction_mean_k3 = 47
scan_direction_mean_k4 = 48
phot_g_n_obs = 49
phot_g_mean_flux = 50
phot_g_mean_flux_error = 51
phot_g_mean_mag = 52
phot_variable_flag = 53
l = 54
b = 55
ecl_lon = 56
ecl_lat = 57

filename = os.path.join(dataDir, fileNameRoot % (0, 0, 0))

with pyfits.open(filename) as hdulist:
    print 'hdulist.info() = ',hdulist.info()
    print 'hdulist[0].header[solution_id] = ',hdulist[0].header[solution_id]
    print 'hdulist[1].header = ',hdulist[1].header
    print 'hdulist[1].header[TTYPE1] = ',hdulist[1].header['TTYPE1']
    print 'hdulist[solution_id=',solution_id,'].data = ',hdulist[solution_id].data
    print 'hdulist[0].data.shape = ',hdulist[0].data.shape
    print 'type(hdulist[0].data[0]) = ',type(hdulist[0].data[0])
    data = hdulist[1].data
    print 'data[0] = ',data[0]
    print 'data[0][0] = ',data[0][0]
    ls = np.array([dat[l-1] for dat in data])
    print 'ls = ',ls
    print 'type(ls) = ',type(ls)
    print 'type(ls[0]) = ',type(ls[0])
    print 'ls.shape = ',ls.shape
    print 'np.min(hdulist[1].data[:][l]) = ',np.min(ls)
    print 'np.max(hdulist[1].data[:][l]) = ',np.max(ls)
    bs = np.array([dat[b-1] for dat in data])
    print 'np.min(hdulist[1].data[:][b]) = ',np.min(bs)
    print 'np.max(hdulist[1].data[:][b]) = ',np.max(bs)

longitudes = np.arange(0,361,10)
print 'longitudes = ',longitudes
latitudes = np.arange(-90,91,10)
print 'latitudes = ',latitudes
STOP
#    file_content = f.read()
#    file_lines = file_content.splitlines()
#    header = file_lines[0].split(',')
#print 'header = ',header

fileLines = []
dataList = []
previous = False
prePrevious = False
isThere = False
for iRelease in range(3):
    for iTract in np.arange(0, 1000):
        for iPatch in np.arange(0, 1000):
            prePrevious = previous
            previous = isThere
            filename = os.path.join(dataDir, fileNameRoot % (iRelease, iTract, iPatch))
            if os.path.exists(filename):
                isThere = True
                if prePrevious and not previous:
                    print 'previous file is missing'
                with gzip.open(filename, 'rb') as f:
                    file_content = f.read()
                    file_lines = file_content.splitlines()
                    header = file_lines[0].split(',')
                    del file_lines[0]
    #                fileLines.append(file_lines)
                    for line in file_lines:
                        line_split = line.split(',')
#                        print 'line_split = ',line_split
                        data = {}
                        for iCol in range(len(header)):
                            data[header[iCol]] = line_split[iCol]
#                            print 'data[',header[iCol],'] = ',data[header[iCol]]
    #                    print 'data = ',data
#                        print "data['phot_g_mean_mag'] = ",data['phot_g_mean_mag']
                        dataList.append(np.float32(data['phot_g_mean_mag']))
                print 'len(dataList) = ',len(dataList)
                print 'filename ',filename,' read'
#print 'dataList = ',dataList
#print 'type(dataList) = ',type(dataList)
#print 'type(dataList[0]) = ',type(dataList[0])
print 'np.min(dataList) = ',np.min(dataList)
print 'np.max(dataList) = ',np.max(dataList)
#dataList = []
#for line in fileLines:
#    line_split = line.split(',')
#    data = {}
#    for iCol in range(len(header)):
#        data[header[iCol]] = line[iCol]
#    print 'data = ',data
#    dataList.append(data)
#print 'len(dataList) = ',len(dataList)
