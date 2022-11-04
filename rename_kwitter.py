from pickle import STOP
from drUtils import getHeader
import os

fitsListFileName = '/Users/azuri/daten/uni/HKU/HASH/KarenKwitter/new_fits/fitsfiles/fitsfiles.list'
with open(fitsListFileName,'r') as f:
    fitsList = f.readlines()
dates = []
observatories = []
newFileNames = []
for fitsFile in fitsList:
    fitsFile = fitsFile.strip()
    header = getHeader(fitsFile,0)
    print('header = ',header)
    dateObs = header['DATE-OBS']
    print('dateObs = ',dateObs)
    date_obs = dateObs
    if 'T' in dateObs:
        date_obs = dateObs[:dateObs.find('T')]
    print('date_obs = ',date_obs)
    try:
        year,month,day = date_obs.split('-')
        date_obs = day+month+year[2:]
    except:
        day,month,year = date_obs.split('/')
        date_obs = day+month+year
    print(date_obs)
    if len(date_obs) < 6:
        raise('wrong date')
    dates.append(date_obs)
    try:
        observatory = header['OBSERVAT']
    except:
        print('dir(header) = ',dir(header))
        print('ERROR: did not find OBSERVAT in header. header.keys = ',header.cards)
    if observatory == 'KPNO':
        observatory = 'KP'
    elif observatory == 'APO':
        observatory = 'AP'
    elif observatory == 'CTIO':
        observatory = 'CT'
    print('observatory = ',observatory)
    if observatory not in observatories:
        observatories.append(observatory)
    newFileName = fitsFile[:fitsFile.rfind('.')]+'_'+observatory+date_obs+'.fits'
    print('newFileName = <'+newFileName+'>')
    if newFileName in newFileNames:
        raise('name already exists')
    newFileNames.append(newFileName)
    os.rename(fitsFile,newFileName)
print('dates = ',dates)
print('observatories = ',observatories)
print('newFileNames = ',newFileNames)
