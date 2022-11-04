import astropy.io.fits as pyfits
import lineid_plot
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import numpy as np
from drUtils import getWavelengthArr,getImageData,readMyPNLineList
import csvFree,csvData

csFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/hash_PNMain_tbCSCoords_240222.csv'

fitsFileName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/IPHASGTC_Ju1_GT110616_id4408.fits'#IPHASGTC_KK26_GT090316_id4393.fits'
lineListName = '/Users/azuri/entwicklung/python/data_reduction/pnLineList_Ju1.dat'#KK26.dat'

wLen = getWavelengthArr(fitsFileName,0)
spec = getImageData(fitsFileName,0)

lineWave,lineLabel = readMyPNLineList(lineListName)

fig, ax = plt.subplots()
ax.plot(wLen,spec)

ax.set_xlabel('Wavelength [$\mathrm{\AA}$]')
ax.set_ylabel('$\mathrm{F}_\lambda$ [$\mathrm{erg}$ $\mathrm{s^{-1}} \mathrm{cm^{-2}} \mathrm{\AA^{-1}}]$')
ax.set_xlim([4000,7200])
ax.set_ylim([-0.1e-15,3.1e-15])
ak = lineid_plot.initial_annotate_kwargs()
ak['arrowprops']['arrowstyle'] = "->"

pk = lineid_plot.initial_plot_kwargs()
pk['color'] = "red"

lineid_plot.plot_line_ids(
    wLen,spec,lineWave,lineLabel,ax=ax,arrow_tip=2.4e-15,box_loc=2.7e-15, annotate_kwargs=ak, plot_kwargs=pk)
#lineid_plot.plot_line_ids(wLen,spec,lineWave,lineLabel,ax=ax,arrow_tip=2.4e-15,box_loc=2.7e-15)
plt.savefig(fitsFileName[:fitsFileName.rfind('.')]+'.pdf', bbox_inches='tight')
plt.show()


css = csvFree.readCSVFile(csFileName)
print('css.header = ',css.header)
hashIDsCS = css.getData('idPNMain')

tableFName = '/Users/azuri/daten/uni/HKU/IPHAS-GTC/table.dat'
with open(tableFName,'r') as f:
    lines = f.readlines()

names = []
pngs = []
ids = []
cs = []
for line in lines:
    line = line.split('&')
    names.append(line[0].strip())
    pngs.append(line[1].strip())
    ids.append(line[2].strip())
    if ids[len(ids)-1] in hashIDsCS:
        cs.append('y')
    else:
        cs.append('n')
    print('names[',len(names)-1,'] = ',names[len(names)-1],': PNG '+pngs[len(pngs)-1],', CS '+cs[len(cs)-1])
