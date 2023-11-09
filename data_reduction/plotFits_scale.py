from spectraUtils import getImageData
import matplotlib.pyplot as plt
import numpy as np
from tkinter import simpledialog
import time
from matplotlib.widgets import Slider

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
        try:
            dlg = simpledialog.SimpleDialog(None, 'Please Select', ['o', 'l', 'u'])
            typessss=dlg.go()
            if typessss==0:
                typesss='o'
            elif typessss==1:
                typesss = 'l'
            elif typessss == 2:
                typesss = 'u'
            print(typessss)
            with open('data.txt', 'a') as f:  # set document object
                f.write(typesss + " " + str(int(event.xdata)) + " " + str(int(event.ydata)) + "\n")
                f.close()
        except:
            return




fig,ax = plt.subplots()
max_val=0
min_val=0

def update_max(val):
    print(val)
    max_val=val
    ax.clear()
    ax.imshow(image, vmin=min_val, vmax=max_val)
    fig.canvas.draw_idle()

def update_min(val):
    print(val)
    min_val=val
    ax.clear()
    ax.imshow(image, vmin=min_val, vmax=max_val)
    fig.canvas.draw_idle()



vmax = np.max([1.5 * np.mean(image), 1.])
vmin=np.min([1.5 * np.mean(image), 1.])
plt.imshow(image, vmin=vmin, vmax=vmax)
fig.canvas.mpl_connect('button_press_event', on_press)

#Edit the parameters
d_val=(vmax-vmin)*1

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)
sfreq = Slider(axfreq, 'vmax', vmax-d_val, vmax,valfmt='% .2f', valinit=0, valstep=0.01)
sfreq.on_changed(update_max)
sfreq.reset()
axfreq1 = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor=axcolor)
sfreq1 = Slider(axfreq1, 'vmin', vmax-d_val, vmax,valfmt='% .2f', valinit=0, valstep=0.01)
sfreq1.on_changed(update_min)
sfreq1.reset()
max_val=vmax
min_val=vmin
sfreq.set_val(vmax)
sfreq1.set_val(min_val)

plt.show()
