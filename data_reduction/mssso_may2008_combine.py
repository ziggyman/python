import numpy as np
from drUtils import scombine
date = '2008-05-08'
year = date[2:4]
month = date[5:7]
day = date[-2:]
print('day = <'+day+'>, month = <'+month+'>, year = <'+year+'>')
pathRoot = '/Users/azuri/spectra/MSO/MSSSO_2m3_DBS_may08/RAW/CentralStars/'

with open(pathRoot+'B_data/'+date+'/areas.csv','r') as f:
    b_data = f.readlines()

with open(pathRoot+'R_data/'+date+'/areas.csv','r') as f:
    r_data = f.readlines()

if len(b_data) != len(r_data):
    print('len(b_data) = ',len(b_data),' != len(r_data) = ',len(r_data))
    STOP

for i in np.arange(1,len(b_data),1):
    if 'SCIENCE' in b_data[i]:
        b_objName = b_data[i][b_data[i].find('SCIENCE_')+8:]
#        print('b_objName = <'+b_objName+'>')
        b_objName = b_objName[:b_objName.find('_')]
#        print('b_objName = <'+b_objName+'>')

        r_objName = r_data[i][r_data[i].find('SCIENCE_')+8:]
#        print('r_objName = <'+r_objName+'>')
        r_objName = r_objName[:r_objName.find('_')]
#        print('r_objName = <'+r_objName+'>')

        if b_objName != r_objName:
            print('b_objName = ',b_objName,' != r_objName = ',r_objName)
            STOP

        with open(pathRoot+'combine.list','w') as f:
            f.write(b_data[i][:b_data[i].find('.fits')]+'-skyEcdF_clean.fits'+'\n')
            f.write(r_data[i][:r_data[i].find('.fits')]+'-skyEcdF_clean.fits'+'\n')

        scombine(pathRoot+'combine.list',
                 pathRoot+r_objName+'_MS'+day+month+year+'.fits',
                 display=True)
