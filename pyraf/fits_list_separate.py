import astropy.io.fits as pyfits
import os

inList = '/Volumes/obiwan/azuri/spectra/saao/saao_may2019/20190501/allFits.list'
suffixes = ['','_z','_zf']

def silentRemove(filename):
    if os.path.exists(filename): os.remove(filename)

lines = []
path = inList[:inList.rfind('/')]
with open(inList,'r') as f:
    lines = [os.path.join(path,line.strip('\n')) for line in f]

listOutName = ''

for suffix in suffixes:
    objectNames = []
    nonZerosOutName = os.path.join(path,'nonZeros'+suffix+'.list')
    nonFlatsOutName = os.path.join(path,'nonFlats'+suffix+'.list')
    for filename in [nonZerosOutName, nonFlatsOutName]:
        silentRemove(filename)
    for line in lines:
        hdulist = pyfits.open(line)
        objectName = hdulist[0].header['OBJECT']
        print('line = ',line,': object = ',objectName)
        if objectName not in objectNames:
            objectNames.append(objectName)
            listOutName = os.path.join(path,objectName+suffix+'.list')
            print('writing to listOutName = <'+listOutName+'>')
            silentRemove(listOutName)
        with open(listOutName,'a') as f:
            f.write(line[:line.rfind('.')]+suffix+'.fits'+'\n')
        if objectName != 'Bias':
            with open(nonZerosOutName,'a') as f:
                f.write(line[:line.rfind('.')]+suffix+'.fits'+'\n')
        if objectName not in ['DomeFlat','Skyflat','Bias']:
            with open(nonFlatsOutName,'a') as f:
                f.write(line[:line.rfind('.')]+suffix+'.fits'+'\n')
