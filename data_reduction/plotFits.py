from spectraUtils import getImageData
import matplotlib.pyplot as plt
import numpy as np
from tkinter import simpledialog
import time
from matplotlib.widgets import Slider,Button

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
            dlg = simpledialog.SimpleDialog(None, 'please select', ['o', 'l', 'u'])
            typessss=dlg.go()
            if typessss==0:
                typesss='o'
            elif typessss==1:
                typesss = 'l'
            elif typessss == 2:
                typesss = 'u'
            print(typessss)
            with open('data.txt', 'a') as f:
                print('event.xdata',event.xdata)
                print('event.ydata', event.ydata)
                f.write(typesss + " " + str(int(event.xdata)) + " " + str(int(event.ydata)) + "\n")
                f.close()
        except:
            return



fig,ax = plt.subplots()
# plt.subplots_adjust(bottom=0.2)
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

def turn_coler_withY(a):
    with open('data.txt', 'r') as f:
        f_line_s=f.readlines()
        datas_y=[]
        for i in f_line_s:
            row=i[:len(i)-1].split(' ')
            row1=[]
            row1.append(int(row[2]))
            row1.append(int(row[1]))
            if row[0]=='l':
                row1.append(0)
            elif row[0]=='o':
                row1.append(1)
            elif row[0]=='u':
                row1.append(2)
            datas_y.append(row1)
        datas_y_sort=np.array(sorted(datas_y,key=lambda y:y[0]))
        print('datas_y_sort[:,2]',datas_y_sort[:,2])

        ip_array=datas_y_sort[:, 2][:-1]-datas_y_sort[:,2][1:]
        print('ip_array',ip_array)
        ip_list=np.where(ip_array==0)[0]
        for i in ip_list:
            if datas_y_sort[i][2]==0:
                image[datas_y_sort[i][0]:datas_y_sort[i+1][0]]=np.ones(image[datas_y_sort[i][0]:datas_y_sort[i+1][0]].shape)*max_val*1.2
            elif datas_y_sort[i][2]==1:
                image[datas_y_sort[i][0]:datas_y_sort[i+1][0]]=np.ones(image[datas_y_sort[i][0]:datas_y_sort[i+1][0]].shape)*(max_val+min_val)/2
            if datas_y_sort[i][2]==2:
                image[datas_y_sort[i][0]:datas_y_sort[i+1][0]]=np.ones(image[datas_y_sort[i][0]:datas_y_sort[i+1][0]].shape)*(max_val+min_val)/4
        f.close()
        ax.clear()
        ax.imshow(image, vmin=min_val, vmax=max_val)
        fig.canvas.draw_idle()

vmax = np.max([1.5 * np.mean(image), 1.])
vmin=np.min([1.5 * np.mean(image), 1.])
plt.imshow(image, origin='lower', vmin=0., vmax=vmax)
# print(image[0])
fig.canvas.mpl_connect('button_press_event', on_press)


d_val=(vmax-vmin)*1

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.15, 0.05, 0.55, 0.03], facecolor=axcolor)
sfreq = Slider(axfreq, 'vmax', vmin, vmax,valfmt='% .2f', valinit=0, valstep=0.01)
sfreq.on_changed(update_max)
sfreq.reset()
axfreq1 = plt.axes([0.15, 0.1, 0.55, 0.03], facecolor=axcolor)
sfreq1 = Slider(axfreq1, 'vmin', vmin, vmax,valfmt='% .2f', valinit=0, valstep=0.01)
sfreq1.on_changed(update_min)
sfreq1.reset()
# font = FontProperties(fname=r"C:\WINDOWS\Fonts\STKAITI.TTF", size=14)

ax_normal = plt.axes([0.8, 0.05, 0.1, 0.08])
btn_normal = Button(ax_normal, 'Select')
btn_normal.on_clicked(func=turn_coler_withY)

max_val=vmax
min_val=vmin
sfreq.set_val(vmax)
sfreq1.set_val(min_val)

plt.show()
