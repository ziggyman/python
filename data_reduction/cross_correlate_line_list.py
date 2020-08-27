import numpy as np
from drUtils import xCor, getImageData, xCorFindMinimum

inputLineList = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_gr7_angle16_3_lines_identified_good.dat'
outputLineList = '/Users/azuri/stella/referenceFiles/spupnic/saao_refspec_gr7_angle16_3_may2020_lines_identified_good.dat'

oldArcSpec = '/Users/azuri/spectra/saao/saao_sep2019/20190907/refArc_spupnic_gr7_16_30_otzxfifEc.fits'
newArcSpec = '/Users/azuri/spectra/saao/saao_may2020/20200519/refARC_spupnic_grat7_16_30_otzxf_may2020iEc.fits'

oldSpec = getImageData(oldArcSpec,0)
oldSpec /= np.mean(oldSpec)
newSpec = getImageData(newArcSpec,0)
newSpec /= np.mean(newSpec)

oldX = np.arange(0,len(oldSpec),1)
newX = np.arange(0,len(newSpec),1)

static = [oldX,oldSpec]
moving = [newX[40:newX.shape[0]-40],newSpec[40:newX.shape[0]-40]]

print('static[0][1] = ',static[0][1])
print('moving[0][1] = ',moving[0][1])

result = xCor(static, moving)
print('result = ',result)

print('result[0] = ',result[0])
print('result[1] = ',result[1])

minPos = xCorFindMinimum(result[0], result[1])
print('minPos = ',minPos)

with open(inputLineList,'r') as f:
    lines = f.readlines()
lines = [line.strip('\n').split(' ') for line in lines]
pix = [float(line[0]) for line in lines]
wLen = [float(line[1]) for line  in lines]
print("pix = ",pix)
print('wLen = ',wLen)

newPix = pix - minPos

with open(outputLineList,'w') as f:
    for i in range(len(pix)):
        f.write('%.5f %.5f\n' % (newPix[i], wLen[i]))
        print('%.5f %.5f\n' % (newPix[i], wLen[i]))
