import astropy.io.fits as pyfits
import numpy as np
import os

import csvFree,csvData

saao_path = '/Volumes/work/azuri/spectra/saao'
saao_all_fits_files = os.path.join(saao_path,'saao_all_fits_files.list')

def createLists():
    with open(saao_all_fits_files,'r') as f:
        all_fits_files = f.readlines()
    all_fits_files = [os.path.join(saao_path,x.rstrip('\n')) for x in all_fits_files]

    objects = []
    raDecs = []

    checkObs = 'BMP0715-0826'
    for f in all_fits_files:
        print('f = <'+f+'>')
        hdulist = pyfits.open(f)
        header = hdulist[0].header
    #    print('header = ',dir(header))
    #    print('header = ',header.keys)
        obs = header['OBJECT']
        obs = obs.replace('_',' ')
        obs = obs.replace('-',' ')
        print('header[OBJECT] = <'+obs+'>')
        if ((obs.lower() not in ['arc','dummy','skyflat','domeflat','test','arcSkyFlat','','hartmansequence','bias','"'])
            and (obs not in objects)):
            objects.append(obs)
    #    STOP
        if ('TARG-RA' in header) and ('TARG-DEC' in header):
            raDec = header['TARG-RA']+' '+header['TARG-DEC']
            print('raDec = ',raDec)
            if raDec not in raDecs:
                raDecs.append(raDec)
        if obs.replace(' ','') == checkObs:
            STOP
    print('objects = ',objects)
    with open('/Users/azuri/daten/uni/HKU/observing/already_observed_May062020_RA_DEC_temp.list','w') as f:
        for i in raDecs:
            f.write(i+'\n')

    # compare result to old list
    oldList = '/Users/azuri/daten/uni/HKU/observing/already_observed.list'
    with open(oldList,'r') as f:
        old_fits_files = f.readlines()
    old_fits_files = [x.rstrip('\n') for x in old_fits_files]
    #old_fits_files = [os.path.join(saao_path,x.rstrip('\n')) for x in old_fits_files]


    objects = [o.replace(' ','') for o in objects]
    for obs in old_fits_files:
        if obs.replace(' ','') in objects:
            print('old name <'+obs+'> found in new names')
        else:
            print('old name <'+obs+'> NOT found in new names')
            objects.append(obs)
        if obs.replace(' ','') == checkObs:
            STOP

    with open('/Users/azuri/daten/uni/HKU/observing/already_observed_May062020_temp.list','w') as f:
        for obs in objects:
            f.write(obs.replace(' ','')+'\n')


if __name__ == '__main__':
    createLists()