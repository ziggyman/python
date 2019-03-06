execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...
import imageio

im = rgb2grey(imageio.imread('/Users/azuri/daten/uni/HKU/Kamila/Abell30_spec_im.png'))
print(im.shape)
print(im[100,:])

im = 255.0 - im
print(im[100,:])

emission = np.sum(im[:,0:14], axis=1)
calib = np.sum(im[:,16:91], axis=1)
spec = np.sum(im[:,99:134],axis=1)
print('spec = ',spec.shape,': ',spec)

pix=[[39.0,],[867.852,3834.24],[838.05, 3868.0],[880.441,3811.35]]

#plt.imshow(im)
plt.plot(calib,'g-',label='calib')
plt.plot(emission,'b-',label='emission')
plt.plot(spec,'r-',label='Abell 30')
plt.legend()
plt.show()

im = rgb2grey(imageio.imread('/Users/azuri/daten/uni/HKU/Kamila/Abell30_spec1977.png'))

plt.imshow(im)
plt.show()

imTrace = trace(im,searchRadius=20)
x = np.arange(4342.,6660.,(6660. - 4342.) / imTrace.shape[0])
imTraceMin = np.min(imTrace)
imTraceMax = np.max(imTrace)
print('imTraceMin = ',imTraceMin,', imTraceMax = ',imTraceMax)
y = 0.35 + ((imTrace - imTraceMin) / ((imTraceMax - imTraceMin) / 0.65))
print('x.shape = ',x.shape,', imTrace.shape = ',imTrace.shape)
plt.plot(x,y)
plt.show()

f = open('/Users/azuri/daten/uni/HKU/Kamila/Abell30_spec.txt','w')
for row in range(y.shape[0]):
    f.write('%.6f %.6f\n' % (x[row],y[row]))
f.close()
