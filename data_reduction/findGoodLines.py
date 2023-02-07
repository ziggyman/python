import numpy as np
import os
from drUtils import reidentify,getLineProfiles,getBestLineProfile,findGoodLines,extractSum

path = '/Users/azuri/stella/referenceFiles/dbs'
arc = 'refProfApDef_DBS_BLUE_otzf.fits'
arcInterp = 'refProfApDef_DBS_BLUE_otzfif.fits'

lineProfiles = getLineProfiles(os.path.join(path,arc))
bestLineProfile = getBestLineProfile(lineProfiles)
oneDSpecInterp = extractSum(os.path.join(path,arcInterp),'row')

xSpec = np.arange(0,oneDSpecInterp.shape[0],1.)
print('xSpec = ',xSpec)
findGoodLines(xSpec,
              oneDSpecInterp,
              bestLineProfile,
              outFileNameAllLines='/Users/azuri/spectra/MSSSO_2m3_DBS_aug07/BLUE/night3/refSpec_lines.dat',
              outFileNameGoodLines='/Users/azuri/spectra/MSSSO_2m3_DBS_aug07/BLUE/night3/refSpec_lines_good.dat')
