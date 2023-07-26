
from spectraUtils import getImageData
import matplotlib.pyplot as plt
import numpy as np
from tkinter import simpledialog
import time
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("Qt5Agg")
#    matplotlib.use("TkAgg")

fName = 'image.fits'
image = getImageData(fName)
data_position=[]
global first_time
first_time=0.0

def on_press(event):
    global first_time
    time_space=time.time()-first_time
    print(time_space)
    first_time=time.time()
    if time_space<0.5:
#        try:
        if True:
            dlg = simpledialog.SimpleDialog(None, 'Please Select', ['o', 'l', 'u'])
            typessss=dlg.go()
            if typessss==0:
                typesss='o'
            elif typessss==1:
                typesss = 'l'
            elif typessss == 2:
                typesss = 'u'
            print(typessss)
            with open('data.txt', 'a') as f:
                f.write(typesss + " " + str(int(event.xdata)) + " " + str(int(event.ydata)) + "\n")
                f.close()
#        except:
#            return


vmax = np.max([1.5 * np.mean(image), 1.])
print('vmax = ', vmax)
fig = plt.figure()
plt.imshow(image, vmin=0., vmax=vmax)
fig.canvas.mpl_connect('button_press_event', on_press)
plt.show()
