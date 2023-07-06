import matplotlib.pyplot as plt
import os
from spectraUtils import readFileToArr,getImageData,getWavelengthArr
import numpy as np

path = '/Users/azuri/spectra/MSSSO_2m3_DBS_aug07/'
hashSpectra = readFileToArr(os.path.join(path,'hash_red.list'))
newSpectra = readFileToArr(os.path.join(path,'objects_red.list'))

hashLength = len(hashSpectra)
newLength = len(newSpectra)

if hashLength != newLength:
    for i in range(np.min([hashLength,newLength])):
        print('hashSpectra[',i,'] = ',hashSpectra[i],', newSpectra[',i,'] = ',newSpectra[i])
    print('hashLength = ',hashLength,' != newLength = ',newLength)
    STOP

for i in range(hashLength):
    objectName = hashSpectra[i][hashSpectra[i].rfind('/')+1:]
    objectName = objectName[:objectName.find('_')]
    fName = os.path.join(path,hashSpectra[i])
    #plt.figure(figsize=(14,6))
    fig,axs = plt.subplots(2, figsize=(15, 9))
    image = getImageData(os.path.join(path,newSpectra[i][:newSpectra[i].rfind('Ecd')]+'-sky.fits'),0)
    vmax = np.max([1.5*np.mean(image),1.])
    print('vmax = ',vmax)
    axs[0].imshow(image,vmin=0.,vmax=vmax)
    axs[1].plot(getWavelengthArr(fName,0),getImageData(fName,0))
    newName = newSpectra[i][newSpectra[i].rfind('/')+1:]
    newName = newName[newName.find('_')+1:]
    newName = newName[:newName.find('_')]
    if objectName != newName:
        print('i = ',i,': objectName = <'+objectName+'>')
        print('i = ',i,': newName = ',newName)
        print('names do not agree')
        STOP
    fName = os.path.join(path,newSpectra[i])
    print('fName = ',fName)
    axs[1].plot(getWavelengthArr(fName,0),getImageData(fName,0))
    axs[1].set_title(newName)
    plt.show()
