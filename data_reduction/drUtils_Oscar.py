import astropy.io.fits as pyfits
from PyAstronomy import pyasl
import astropy
from astropy import units as u
from astropy.nddata import CCDData
from specutils import Spectrum1D
from astropy.time import Time
from specreduce import fluxcal,calibration_data
from fluxcal import apply_sensfunc, obs_extinction, airmass_cor
import numpy as np
import os
import matplotlib.pyplot as plt
from drUtils import getHeader,getClosestInTime

def getHeaderValue(fname, keyword, hduNum=0):
    hdulist = pyfits.open(fname)
    header = hdulist[hduNum].header
    hdulist.close()
    try:
        return header[keyword]
    except:
        return None

def getWavelengthArr(fname, hduNum=0, dim=1):
    hdulist = pyfits.open(fname)
    header = hdulist[hduNum].header
    hdulist.close()
    if 'CDELT%d' % (dim) in header.keys():
        cdelt = header['CDELT%d' % (dim)]
        wLen = ((np.arange(header['NAXIS%d' % (dim)]) + 1.0) - header['CRPIX%d' % (dim)]) * cdelt + header['CRVAL%d' % (dim)]
    elif 'CD1_%d' % (dim) in header.keys():
        cdelt = header['CD1_%d' % (dim)]
        wLen = ((np.arange(header['NAXIS%d' % (dim)]) + 1.0) - header['CRPIX%d' % (dim)]) * cdelt + header['CRVAL%d' % (dim)]
    else:
        print('WARNING: neither CDELT1 nor CD1_1 found in header of file <'+fname+'>')
        wLen = np.arange(header['NAXIS%d' % (dim)]) + 1.0
    return wLen

def writeFits1D(flux, outFileName, wavelength=None, header=None, CRVAL1=None, CRPIX1=None, CDELT1=None):
    head = None
    if not header is None:
        if isinstance(header,str):
#            refFileName = header
            hdulist = pyfits.open(header)
            head = hdulist[0].header
            hdulist.close()
        else:
            head = header
        if 'NAXIS2' in head.keys():
            del head['NAXIS2']
        print('writeFits1D: dir(head) = ',dir(head))
        head['NAXIS1'] = flux.shape[0]
        for key in head:
            print('writeFits1D: ',key+': <',head[key],'>: ',type(head[key]))
            if isinstance(head[key],str):
                if '\n' in head[key]:
                    head[key].replace('\n','')
            elif (key == 'COMMENT'):# and (head[key] == ''):
                del head[key]
            elif '\n' in str(head[key]):
                del head[key]
#            elif isinstance(head[key],'astropy.io.fits.header._HeaderCommentaryCards'):
#                del head[key]


    waveParams = None
    if not CRVAL1 is None:
        waveParams = {'CRVAL1': CRVAL1,
                      'CRPIX1': CRPIX1,
                      'CDELT1': CDELT1}
    #print('writeFits1D: flux = ',flux)
    #print('writeFits1D: wavelength = ',wavelength)
    #print('writeFits1D: waveParams = ',waveParams)
    #print('writeFits1D: head = ',head)
    print('writeFits1D: writing file <'+outFileName+'>')
    print('wavelength = ',len(wavelength),': ',wavelength)
    pyasl.write1dFitsSpec(outFileName, flux, wvl=wavelength, waveParams=waveParams, fluxErr=None, header=head, clobber=True, refFileName=None, refFileExt=0)

#@params:
# objectSpectrIn: list of input spectra file names
# objectSpectraOut: list of output spectra file names
# stdsDLam: either list of deltaLambda in same order as sensFuncs or file name which will be read
# sensFuncs: list of sensitivity functions (files can be read with sens = astropy.io.ascii.read('<standard star name>_sens.ecsv'))
# airmassExtCor: file to use for extinction correction depending on air mass
def applySensFuncs(lam,spec,header, objectSpectraOut, sensFuncFName, airmassExtCor='apoextinct.dat'):
    sensFunc = astropy.io.ascii.read(sensFuncFName)
    print('sensFunc = ',sensFunc)
    img = CCDData(spec, unit=u.adu)
    # put in units of ADU/s
    img.data = img.data / float(header['EXPTIME'])
    img.unit = u.adu / u.s
#        print('applySensFuncs: img.header = ',img.header)
#        print('applySensFuncs: dir(img.header) = ',dir(img.header))
#        print('applySensFuncs: img.header.keys = ',img.header.keys)

    wLen = lam * u.angstrom

    obj_spectrum = Spectrum1D(spectral_axis=wLen, flux=img.data * img.unit)#,
#                                 uncertainty=StdDevUncertainty(ex_tbl['fluxerr']))
    Xfile = obs_extinction(airmassExtCor)

    try:
        AIRVAL = float(header['AIRMASS'])#img.header['AIRMASS'])
    except:
        AIRVAL = float(header['SECZ'])
#        print('applySensFuncs: AIRMASS = ',AIRVAL)
#        print('applySensFuncs: obj_spectrum = ',obj_spectrum)
    obj_spectrum = airmass_cor(obj_spectrum, AIRVAL, Xfile)
#        print('applySensFuncs: obj_spectrum = ',obj_spectrum)
#        print('applySensFuncs: sensFuncs[0] = ',sensFuncs[0])
    objectSpectrumFluxCalibrated = apply_sensfunc(obj_spectrum, sensFunc)
#        print('applySensFuncs: objectSpectrumFluxCalibrated.data = ',objectSpectrumFluxCalibrated.data)
#        print('applySensFuncs: dir(objectSpectrumFluxCalibrated.data) = ',dir(objectSpectrumFluxCalibrated.data))
#        print('applySensFuncs: img.header.keys = ',img.header.keys)
    crval = obj_spectrum.wavelength[0]
    cdelt = obj_spectrum.wavelength[1]-obj_spectrum.wavelength[0]
    crpix = 1
#        print('applySensFuncs: crval = ',crval)
#        print('applySensFuncs: dir(crval) = ',dir(crval))
#        print('applySensFuncs: crval.value = ',crval.value)
#        print('applySensFuncs: dir(crval.value) = ',dir(crval.value))
#        print('applySensFuncs: cdelt = ',cdelt)
#        print('applySensFuncs: dir(cdelt) = ',dir(cdelt))
#        print('applySensFuncs: crpix = ',crpix)
#        print('applySensFuncs: objectSpectrumFluxCalibrated.data = ',objectSpectrumFluxCalibrated.data)
    #print('objectSpectraIn = ',len(objectSpectraIn),': ',objectSpectraIn,', iSpec = ',iSpec,', objectsSpectraIn[',iSpec,'] = ',objectSpectraIn[iSpec])
    #print('objectSpectraOut = ',len(objectSpectraOut),': ',objectSpectraOut,', iSpec = ',iSpec)
    #print('applySensFuncs: writing file ',objectSpectraOut)
    #writeFits1D(objectSpectrumFluxCalibrated.data,
    #            objectSpectraOut,
    #            wavelength=obj_spectrum.wavelength,
    #            header=header,
    #            CRVAL1=None,#crval.value,
    #            CRPIX1=None,#crpix,
    #            CDELT1=None#cdelt.value,
    #            )
    if (wLen[0] != obj_spectrum.wavelength[0]) or (wLen[len(wLen)-1] != obj_spectrum.wavelength[len(wLen)-1]) or (len(wLen) != len(obj_spectrum.wavelength)):
        print('ERROR')
        STOP
    return objectSpectrumFluxCalibrated.data

import pickle

# write list to binary file
def store_list(list_name,file_name):
    # store list in binary file so 'wb' mode
    with open(file_name, "wb") as fp:
        pickle.dump(list_name, fp)

# read list from pickle
def read_list(file_name):
    with open(file_name, "rb") as fp:
        return pickle.load(fp)

data_table_b_raw = read_list("/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/data_b raw")
data_table_r_raw = read_list("/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/data_r raw")

objects_b = np.array(read_list("/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/object_b"))
objects_r = np.array(read_list("/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/object_r"))

sensFuncsListName = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/sensFuncs/fluxstds.list'
with open(sensFuncsListName,'r') as f:
    sensFuncsNames = f.readlines()
sensFuncsNames = [s.strip() for s in sensFuncsNames]

all_otzxfifListName = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/sensFuncs/all_otzxfif.list'
with open(all_otzxfifListName,'r') as f:
    all_otzxfif = f.readlines()
all_otzxfif = [a.strip() for a in all_otzxfif]

keyWord = 'DATE-OBS'
found = False
sensTimes = []
for i in range(len(sensFuncsNames)):
    found = False
    for j in range(len(all_otzxfif)):
        if sensFuncsNames[i][sensFuncsNames[i].find('dbs'):sensFuncsNames[i].find('dbs')+8]+'b' in all_otzxfif[j]:
            print('found ',sensFuncsNames[i][sensFuncsNames[i].find('dbs'):sensFuncsNames[i].find('dbs')+8],' in all_otzxfif[',j,'] = ',all_otzxfif[j])
            found = True
            hdulist = pyfits.open(all_otzxfif[j])
            headerSc = hdulist[0].header
        #    for key in headerSc.keys():
        #        print('getClosestArc: headerSc[',key,'] = ',headerSc[key])
            try:
                if keyWord == 'DATE-OBS':
                    specTimeTemp = Time(headerSc[keyWord], format='isot', scale='utc')
                    specTime = specTimeTemp.mjd
                else:
                    specTime = float(headerSc[keyWord])
            except:
                print('getClosestArcs: ERROR: keyWord <'+keyWord+'> not found in '+all_otzxfif[j])
                if keyWord == 'MJD-OBS':
                    try:
                        specTimeTemp = Time(headerSc['DATE-OBS'], format='isot', scale='utc')
                        specTime = specTimeTemp.mjd
                    except:
                        print('getClosestArcs: ERROR: keyWord <DATE-OBS> not found in '+all_otzxfif[j])
                        STOP
                if keyWord == 'HJD-OBS':
                    try:
                        specTimeTemp = Time(headerSc['DATE-OBS'], format='isot', scale='utc')
                        specTime = specTimeTemp.hjd
                    except:
                        print('getClosestArcs: ERROR: keyWord <DATE-OBS> not found in '+all_otzxfif[j])
                        STOP
            hdulist.close()
            sensTimes.append(specTime)
            break
    if not found:
        print('ERROR: could not find sensFuncsNames[',i,'] = ',sensFuncsNames[i])
        STOP

print('sensFuncsNames = ',len(sensFuncsNames),': ',sensFuncsNames)
print('sensTimes = ',len(sensTimes),': ',sensTimes)
#STOP
display = False

if True:
    with open('/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/filenames.txt','r') as f:
        fileNames = f.readlines()
    fNames_blue = []
    fNames_red = []
    for fileName in fileNames:
        if len(fNames_blue) == len(fNames_red):
            fNames_blue.append(fileName.strip())
        else:
            fNames_red.append(fileName.strip())
    print('fNames_blue = ',len(fNames_blue),': ',fNames_blue)
    print('fNames_red = ',len(fNames_red),': ',fNames_red)

    spectraDirBlue = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/B_data/1DspectraBlue/'
    spectraDirRed = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/R_data/1DspectraRed/'

    print('len(data_table_b_raw) = ',len(data_table_b_raw))
    print('data_table_b_raw[0][0] = ',data_table_b_raw[0][0])

    sorted_idx_blue = np.argsort(objects_b)
    print('sorted_idx_blue = ',sorted_idx_blue)
    sorted_idx_red = np.argsort(objects_r)
    print('sorted_idx_red = ',sorted_idx_red)

    fNames_blue_sorted = []
    fNames_red_sorted = []
    for i in range(len(sorted_idx_red)):
        for j in range(len(fNames_blue)):
            if objects_b[sorted_idx_blue[i]] in fNames_blue[j]:
                fNames_blue_sorted.append(fNames_blue[j])
            if objects_r[sorted_idx_red[i]] in fNames_red[j]:
                fNames_red_sorted.append(fNames_red[j])
    fNames_blue = np.array(fNames_blue_sorted)
    print('fNames_blue sorted = ',fNames_blue)

    fNames_red = np.array(fNames_red_sorted)
    print('fNames_red sorted = ',fNames_red)

    data_b = []
    data_r = []
    len_blue = 0
    len_red = 0
    for i in range(len(data_table_b_raw)):
        print('fNames_blue[',i,'] = ',fNames_blue[i],': len(data_table_b_raw[',sorted_idx_blue[i],']) = ',len(data_table_b_raw[sorted_idx_blue[i]]))
        print('fNames_red[',i,'] = ',fNames_red[i],': len(data_table_r_raw[',sorted_idx_red[i],']) = ',len(data_table_r_raw[sorted_idx_red[i]]))
        data_b_raw = [data_table_b_raw[sorted_idx_blue[i]][0]]
        data_r_raw = [data_table_r_raw[sorted_idx_red[i]][0]]
        len_blue += len(data_table_b_raw[sorted_idx_blue[i]])
        len_red += len(data_table_r_raw[sorted_idx_red[i]])
        for j in np.arange(1,len(data_table_b_raw[sorted_idx_blue[i]]),1):
            data_b_raw.append(data_table_b_raw[sorted_idx_blue[i]][j].data)
            data_r_raw.append(data_table_r_raw[sorted_idx_red[i]][j].data)
        data_b.append(data_b_raw)
        data_r.append(data_r_raw)
    print('len_blue = ',len_blue)
    print('len_red = ',len_red)

    if display:
        for i in range(len(data_table_b_raw)):
            fig,(ax1,ax2) = plt.subplots(1,2,figsize=(20,10))
            ax1.set_title(fNames_blue[i])
            ax2.set_title(fNames_red[i])
            vmin = np.min([np.min(data_b[i][1:]),np.min(data_r[i][1:])])
            vmax = np.max([np.max(data_b[i][1:]),np.max(data_r[i][1:])])
            ax1.imshow(data_b[i][1:],vmin=vmin,vmax=vmax)
            ax2.imshow(data_r[i][1:],vmin=vmin,vmax=vmax)
            plt.show()

            plt.plot(data_b[i][0],data_b[i][int(len(data_b[i])/2)])
            plt.plot(data_r[i][0],data_r[i][int(len(data_r[i])/2)])
            plt.title(fNames_blue[i])
            plt.show()


    data_table_b_cal = []
    data_table_r_cal = []

    for i in range(len(fNames_blue)):
        wLen = data_b[i][0]
        header = getHeader(os.path.join(spectraDirBlue,fNames_blue[i][:fNames_blue[i].rfind('-sky')]+'Ecd.fits'))
        try:
            if keyWord == 'DATE-OBS':
                specTimeTemp = Time(header[keyWord], format='isot', scale='utc')
                specTime = specTimeTemp.mjd
            else:
                specTime = float(header[keyWord])
        except:
            print('getClosestArcs: ERROR: keyWord <'+keyWord+'> not found in '+all_otzxfif[j])
            if keyWord == 'MJD-OBS':
                try:
                    specTimeTemp = Time(header['DATE-OBS'], format='isot', scale='utc')
                    specTime = specTimeTemp.mjd
                except:
                    print('getClosestArcs: ERROR: keyWord <DATE-OBS> not found in '+all_otzxfif[j])
                    STOP
            if keyWord == 'HJD-OBS':
                try:
                    specTimeTemp = Time(header['DATE-OBS'], format='isot', scale='utc')
                    specTime = specTimeTemp.hjd
                except:
                    print('getClosestArcs: ERROR: keyWord <DATE-OBS> not found in '+all_otzxfif[j])
                    STOP
        hdulist.close()
        closestSensFunc = getClosestInTime(specTime,sensTimes)
        print('closestSensFunc = ',closestSensFunc)
        print('fNames_blue[i] = ',fNames_blue[i],': closestSensFunc = ',sensFuncsNames[closestSensFunc[0]],': ',sensTimes[closestSensFunc[0]])
        sensFuncName = sensFuncsNames[closestSensFunc[0]]
        if sensFuncName[sensFuncName.find('_dbs')+9] == 'r':
            sensFuncNamels = list(sensFuncName)
            sensFuncNamels[sensFuncName.find('_dbs')+9] = 'b'
            sensFuncName = ''.join(sensFuncNamels)
#        STOP
        for j in np.arange(1,len(data_b[i]),1):
            spec = data_b[i][j]
            print('wLen = ',type(wLen),': ',wLen)
            print('spec = ',type(spec),': ',dir(spec))
            print('i = ',i,': len(data_b[',i,']) = ',len(data_b[i]))
            print('fNames_blue[',i,'] = ',fNames_blue[i])
#            if display and (j == int(len(data_table_b_raw[i])/2) and i == 0):
            if display and (i == 0):
                plt.plot(wLen,spec)
                plt.show()
            data_table_b_raw[sorted_idx_blue[i]][j].data = applySensFuncs(wLen,spec,header,os.path.join(spectraDirBlue,fNames_blue[i][:fNames_blue[i].rfind('.')])+'dF%d.fits' % (j-1),sensFuncName)#'/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/blue/SCIENCE_LTT7379_dbs00811b_otzxfif_sens.ecsv')#SCIENCE_EG274_dbs01153b_otzxfif_sens.ecsv')
#            if display and (j == int(len(data_table_b_raw[i])/2) and i == 0):
            if display and (i == 0):
                plt.plot(wLen,data_table_b_raw[sorted_idx_blue[i]][j].data)
                plt.show()
        data_table_b_cal.append(data_table_b_raw[sorted_idx_blue[i]])

        wLen = data_r[i][0]
        sensFuncName = sensFuncsNames[closestSensFunc[0]]
        if sensFuncName[sensFuncName.find('_dbs')+9] == 'b':
            sensFuncNamels = list(sensFuncName)
            sensFuncNamels[sensFuncName.find('_dbs')+9] = 'r'
            sensFuncName = ''.join(sensFuncNamels)
        header = getHeader(os.path.join(spectraDirRed,fNames_red[i][:fNames_red[i].rfind('-sky')]+'Ecd.fits'))
        for j in np.arange(1,len(data_r[i]),1):
            #plt.plot(wLen,data_table_r_raw[i][j])
            spec = data_r[i][j]
            data_table_r_raw[sorted_idx_red[i]][j].data = applySensFuncs(wLen,spec,header,os.path.join(spectraDirRed,fNames_red[i][:fNames_red[i].rfind('.')])+'dF%d.fits' % (j-1),'/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/red/SCIENCE_LTT7379_dbs00811r_otzxfif_sens.ecsv')#SCIENCE_EG274_dbs01153r_otzxfif_sens.ecsv')
            #plt.plot(wLen,data_table_r_raw[i][j])
            #plt.show()
        data_table_r_cal.append(data_table_r_raw[sorted_idx_red[i]])

    store_list(data_table_b_cal,"/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/data_b_flux_calibrated.pkl")
    store_list(data_table_r_cal,"/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/data_r_flux_calibrated.pkl")

with open('/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/filenames_sorted.txt','w') as f:
    for i in range(len(fNames_blue)):
        f.write(fNames_blue[i]+'\n')

data_table_b_cal = read_list("/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/data_b_flux_calibrated.pkl")
print('len(data_table_b_raw[0]) = ',len(data_table_b_raw[0]))
print('len(data_table_b_cal[0]) = ',len(data_table_b_cal[0]))
data_table_r_cal = read_list("/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/data_r_flux_calibrated.pkl")

data_b = []
data_r = []
for i in range(len(data_table_b_cal)):
    print('fNames_blue[',i,'] = ',fNames_blue[i],': len(data_table_b_cal[',i,']) = ',len(data_table_b_cal[i]))
    print('fNames_red[',i,'] = ',fNames_red[i],': len(data_table_r_cal[',i,']) = ',len(data_table_r_cal[i]))
    data_b_cal = [data_table_b_cal[i][0]]
    data_r_cal = [data_table_r_cal[sorted_idx_red[i]][0]]
    len_blue += len(data_table_b_cal[i])
    len_red += len(data_table_r_cal[i])
    for j in np.arange(1,len(data_table_b_cal[i]),1):
        data_b_cal.append(data_table_b_cal[i][j].data)
        data_r_cal.append(data_table_r_cal[i][j].data)
    data_b.append(data_b_cal)
    data_r.append(data_r_cal)

for i in range(len(data_b)):
    fig,axs = plt.subplots(3,2,figsize=(15,10))
    axs[0,0].set_title(fNames_blue[i])
    axs[0,1].set_title(fNames_red[i])
    print('len(data_b[i][1]) = ',len(data_b[i][1]))
    min_b = np.min(np.array(np.min([x[1750:2000] for x in data_b[i][1:]])))
    min_r = np.min(np.array(np.min([x[500:750] for x in data_r[i][1:]])))
    print('min_b = ',min_b)
    vmin = np.min([min_b,min_r])

    max_b = np.max(np.array(np.max([x[1750:2000] for x in data_b[i][1:]])))
    max_r = np.max(np.array(np.max([x[500:750] for x in data_r[i][1:]])))
    vmax = np.max([max_b,max_r])
    spec_b = [fl[1750:2000] for fl in data_b[i][1:]]
    spec_r = [fl[500:750] for fl in data_r[i][1:]]
    axs[0,0].imshow(spec_b,vmin=vmin,vmax=vmax)
    axs[0,1].imshow(spec_r,vmin=vmin,vmax=vmax)

    specSum_blue = np.sum(data_b[i][1:],0)
    specSum_red = np.sum(data_r[i][1:],0)

    axs[1,0].plot(specSum_blue)
    axs[1,1].plot(specSum_red)
    axs[1,0].set_ylim(np.min([np.min(specSum_blue),np.min(specSum_red)]),np.max([np.max(specSum_blue),np.max(specSum_red)]))
    axs[1,1].set_ylim(np.min([np.min(specSum_blue),np.min(specSum_red)]),np.max([np.max(specSum_blue),np.max(specSum_red)]))

    collapsed_blue = np.sum(data_b[i][1:],1)
    collapsed_red = np.sum(data_r[i][1:],1)
    axs[2,0].plot(collapsed_blue)
    axs[2,1].plot(collapsed_red)
    axs[2,0].set_ylim(np.min([np.min(collapsed_blue),np.min(collapsed_red)]),np.max([np.max(collapsed_blue),np.max(collapsed_red)]))
    axs[2,1].set_ylim(np.min([np.min(collapsed_blue),np.min(collapsed_red)]),np.max([np.max(collapsed_blue),np.max(collapsed_red)]))
    plt.savefig('/Users/azuri/daten/uni/HKU/interns_projects/oscar/2D/fluxcal/'+fNames_blue[i][:fNames_blue[i].rfind('.')]+'.png')
    plt.show()

#    plt.plot(data_table_b_cal[i][0],data_table_b_cal[i][int(len(data_table_b_cal[i])/2)])
#    plt.plot(data_table_r_cal[i][0],data_table_r_cal[i][int(len(data_table_r_cal[i])/2)])
#    plt.title(fNames_blue[i])
#    plt.show()
