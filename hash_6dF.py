import csv
import numpy as np
import os

path = '/Users/azuri/spectra/6dF/556/hash/6dF_Aug2003'
fitsFile = '/Users/azuri/daten/uni/HKU/HASH/hash_fitsfiles_020221.csv'
sqlFileOut = '/Users/azuri/daten/uni/HKU/HASH/hash_do.sql'

idRange = np.arange(10382,10449,1)
ids = [3079,3103,11311,3119,3077,3126,3093,3102,3100,3141,3088,3112,3085,3113,3140,3117,3111,3139,3094,3158,3096,11317,3120,3155,3107,3193,3162,3142,3153,3191,3184,3169,11323,3122,]

fits = csv.DictReader(open(fitsFile))

with open(sqlFileOut,'w') as f:
    for fitsRow in fits:
        idFits = fitsRow['idFitsFiles']
        if int(idFits) in idRange:
            if int(fitsRow['idPNMain']) in ids:
                f.write("DELETE FROM `PNSpectra_Sources`.`FitsFiles` WHERE `idFitsFiles`='"+idFits+"';\n")
                fName = fitsRow['fileName']
                print('removing file ',fName)
                os.remove(os.path.join(path,fName))
