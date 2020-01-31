import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

imA = '/Users/azuri/daten/uni/HKU/line_ratios/Wyse-1942-plateXII.png'
imB = '/Users/azuri/daten/uni/HKU/line_ratios/Wyse-1942-plateXIII.png'

im = Image.open(imA) # Can be many different formats.
pix = im.load()
print(im.size)  # Get the width and hight of the image for iterating over
print(pix[0,0])  # Get the RGBA Value of the a pixel of an image
plt.imshow(im)
plt.show()

print('pix[',im.size[0]-1,',',im.size[1]-1,'] = ',pix[im.size[0]-1,im.size[1]-1])


ratio5007Hbeta = []
ratio50074959 = []
for row in np.arange(20,75,1):
    r = []
    for col in np.arange(0,im.size[0],1):
#        print('col = ',col)
        r.append(pix[int(col),int(row)][0])
    r = 255-np.array(r)
    plt.plot(r)
    if (np.amax(r[106:122]) < 254) and (np.amax(r[50:66]) < 254):
        ratio5007Hbeta.append(np.amax(r[106:122]) / np.amax(r[50:66]))
    else:
        ratio5007Hbeta.append(0)
    if (np.amax(r[106:122]) < 254) and (np.amax(r[88:102]) < 254):
        ratio50074959.append(np.amax(r[106:122]) / np.amax(r[88:102]))
    else:
        ratio50074959.append(0)
    print('row = ',row,': 5007 = ',r[106:122])
    print('row = ',row,': Hbeta = ',r[50:66])
    print('row = ',row,': 4959 = ',r[88:102])
plt.show()
plt.plot(ratio5007Hbeta,label='plateXII 5007/Hbeta')
plt.legend()
plt.show()
plt.plot(ratio50074959,label='plateXII 5007/4959')
plt.legend()
plt.show()

img=mpimg.imread(imB)
print('img = ',img)
imgplot = plt.imshow(img)
plt.show()

ratio5007Hbeta = []
ratio50074959 = []
for row in np.arange(15,79,1):
    r = []
    for col in np.arange(227,430,1):
        r.append(img[int(row),int(col)][0])
    r = 255 - np.array(r)
    plt.plot(r)
    if (np.amax(r[175:200]) < 254) and (np.amax(r[4:28]) < 254):
        ratio5007Hbeta.append(np.amax(r[175:200]) / np.amax(r[4:28]))
    else:
        ratio5007Hbeta.append(0)
    if (np.amax(r[175:200]) < 254) and (np.amax(r[122:140]) < 254):
        ratio50074959.append(np.amax(r[175:200]) / np.amax(r[122:140]))
    else:
        ratio50074959.append(0)
    print('row = ',row,': 5007 = ',r[175:200])
    print('row = ',row,': Hbeta = ',r[4:28])
    print('row = ',row,': 4959 = ',r[122:140])
plt.show()

plt.plot(ratio5007Hbeta,label='plateXIII 5007/Hbeta')
plt.legend()
plt.show()
plt.plot(ratio50074959,label='plateXIII 5007/4959')
plt.legend()
plt.show()

