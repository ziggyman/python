import numpy as np
import matplotlib.pyplot as plt
from lineid_plot import plot_line_ids
import csvFree,csvData
from fits_fit_2gauss import gauss

lineListFileName = '/Users/azuri/stella/referenceFiles/dbs/dbs_near.dat'
data = csvFree.readCSVFile(lineListFileName)

nistFileName = '/Users/azuri/stella/referenceFiles/dbs/NIST_lines.csv'
nist = csvFree.readCSVFile(nistFileName)
print('nist.header = ',nist.header)
nistWLen = nist.getData('obs_wl_air(A)')
nistIntens = nist.getData('intens')
print('nistWLen = ',nistWLen)
print('nistIntens = ',nistIntens)
i = 0
while i < len(nistWLen):
    try:
        nistWLen[i] = float(nistWLen[i].replace('=',''))
        nistIntens[i] = float(nistIntens[i].replace('=','').replace('w','').replace('*','').replace('h','').replace('l','').replace('s','').replace(';','').replace('i','').replace('?','').replace('b','').replace('p','').replace('c','').replace('a','').replace('f',''))
        i += 1
    except:
        print('nistWLen[',i,'] = ',nistWLen[i],', nistIntens[',i,'] = ',nistIntens[i])
        del nistWLen[i]
        del nistIntens[i]
#        i -= 1

x = np.arange(4000.,6000.,1.)
intens = np.zeros(x.shape)

for i in range(len(nistWLen)):
    intens += gauss(x,nistIntens[i],nistWLen[i],1.)

wLen = [float(x) for x in data.getData('lambda')]

for x0 in wLen:
    intens += gauss(x,1.,x0,1.)
#    plt.plot(x,intens)
#    plt.show()
#    STOP

#plt.plot(x,intens)
plot_line_ids(x,intens,wLen,data.getData('lambda'))
plt.show()
