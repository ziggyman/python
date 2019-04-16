#execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...
import imageio
import numpy as np
import matplotlib.pyplot as plt

imFile1964 = '/Users/azuri/daten/uni/HKU/Kamila/Abell30_spec_im.png'
imFile1977 = '/Users/azuri/daten/uni/HKU/Kamila/Abell30_spec1977.png'

def rgb2grey(rgb):

    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray

def trace(image, searchRadius=5):
    axis = 0
    print('image.shape = ', image.shape)

    minPos = []#np.argmin(image,axis=1)
    startFound = False
    minIndex = 0
    for row in range(image.shape[axis]):
        addValue = True
        if not startFound:
            minIndex = np.argmin(image[row,:])
            print('row = ',row,': minIndex = ',minIndex,': image[',row,',',minIndex,'] = ',image[row,minIndex])
            if image[row,minIndex] < 10:
                startFound = True
                minPos.append(minIndex)
        else:
            searchStart = minPos[len(minPos)-1]-searchRadius
            if searchStart < 0:
                searchStart = 0

            searchEnd = minPos[len(minPos)-1]+searchRadius
            if searchEnd >= image.shape[axis]:
                searchEnd = image.shape[axis]-1
            print('row ',row,': searchStart = ',searchStart,', searchEnd = ',searchEnd,': ',image[row,searchStart:searchEnd])
            minIndex = np.argmin(image[row,searchStart:searchEnd])
            print('minIndex = ',minIndex)
            searchRadiusTemp = searchRadius
            iRun = 0
            while (minIndex == 0) and (image[row,searchStart + minIndex] == image[row,searchStart + minIndex + int(searchRadius/2)]):
                searchRadiusTemp = 2 * searchRadiusTemp
                searchStart = minPos[len(minPos)-1]-searchRadiusTemp
                if searchStart < 0:
                    searchStart = 0

                searchEnd = minPos[len(minPos)-1]+searchRadiusTemp
                if searchEnd >= image.shape[axis]:
                    searchEnd = image.shape[axis]-1
                print('row ',row,': searchStart = ',searchStart,', searchEnd = ',searchEnd,': ',image[row,searchStart:searchEnd])
                minIndex = np.argmin(image[row,searchStart:searchEnd])
                print('minIndex = ',minIndex)
                iRun += 1
                if iRun > 10:
                    addValue = False
                    break
            if addValue:
                minPos.append(searchStart + minIndex)
        if len(minPos) > 1:
            print('minPos[',len(minPos)-1,'] = ',minPos[len(minPos)-1])

    print('minPos = ',len(minPos),': ',minPos)
    minPos = np.array(minPos)
    minPos = image.shape[1] - np.flip(minPos,0)
    print('minPos = ',minPos.shape,': ',minPos)
    plt.plot(minPos)
    plt.show()
#    for row in range(image.shape[0]):
#        minPos.append()
#    print(image[:,100])
    return minPos

im = rgb2grey(imageio.imread(imFile1964))
print(im.shape)
print(im[100,:])

im = 255.0 - im
print(im[100,:])

emission = np.flip(np.sum(im[:,0:14], axis=1),0)
calib = np.flip(np.sum(im[:,16:91], axis=1),0)
spec = np.flip(np.sum(im[:,99:134],axis=1),0)
print('spec = ',spec.shape,': ',spec)

pix=[[39.20,5007.],[71.9, 4959.],[269.6, 4686.],[293.5,4650.],[439.25,4442.],[838.05, 3868.0],[867.852,3834.24],[880.441,3811.35],[1145.5,3434.]]
x = np.array([a[0] for a in pix])
print('x = ',type(x),': ',x)
y = np.array([a[1] for a in pix])
coeffs = np.polyfit(x,y,5,full=True)
print('polynomial fit: ',coeffs)
wavelength = np.flip(np.polyval(coeffs[0],range(len(spec))),0)

#3434, doublet, 3868, 4442, 4650, 4686, 4959, 5007

#plt.imshow(im)
plt.plot(wavelength, calib,'g-',label='calib')
plt.plot(wavelength, emission,'b-',label='emission')
plt.plot(wavelength, spec,'r-',label='Abell 30')
plt.legend()
plt.show()

f = open(imFile1964[0:imFile1964.rfind('.')]+'.txt','w')
for row in range(spec.shape[0]):
    f.write('%.6f %.6f\n' % (wavelength[row],spec[row]))
f.close()


# 1977 publication
im = rgb2grey(imageio.imread(imFile1977))

plt.imshow(im)
plt.show()

imTrace = trace(im,searchRadius=20)
wavelengthRange = [4342.,6660.]
x = np.arange(wavelengthRange[0],wavelengthRange[1],(wavelengthRange[1] - wavelengthRange[0]) / imTrace.shape[0])
imTraceMin = np.min(imTrace)
imTraceMax = np.max(imTrace)
print('imTraceMin = ',imTraceMin,', imTraceMax = ',imTraceMax)
y = 0.35 + ((imTrace - imTraceMin) / ((imTraceMax - imTraceMin) / 0.65))
print('x.shape = ',x.shape,', imTrace.shape = ',imTrace.shape)
plt.plot(x,y)
plt.show()

f = open(imFile1977[0:imFile1977.rfind('.')]+'.txt','w')
for row in range(y.shape[0]):
    f.write('%.6f %.6f\n' % (x[row],y[row]))
f.close()
