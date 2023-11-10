import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets
import matplotlib.patches as patches

import numpy as np

from spectraUtils import getImageData


cid = None
rid = None
coords = []
skyBelow = []
skyAbove = []
obs = []
rectSkyBelow = None
rectSkyAbove = None
rectObs = None

def onKeyPress(event):
    global cid
    global rid
    global coords

    print('You pressed the %s button' % (event.key))

    def onClick(event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))
        coords.append([event.xdata,event.ydata])
    def onRelease(event):
        print('%s release: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))
        coords.append([event.xdata,event.ydata])
    if cid is None:
        cid = fig.canvas.mpl_connect('button_press_event', onClick)
    if rid is None:
        rid = fig.canvas.mpl_connect('button_release_event', onRelease)
    print('coords = ',coords)

def onKeyRelease(event):
    global cid
    global rid
    global coords
    global skyBelow
    global skyAbove
    global obs
    global rectSkyAbove
    global rectSkyBelow
    global rectObs
    print('You released the %s button' % (event.key))
    fig.canvas.mpl_disconnect(cid)
    fig.canvas.mpl_disconnect(rid)
    cid = None
    rid = None
    xmin, xmax = ax.get_xlim()
    print('xmin = ',xmin,', xmax = ',xmax)
    if event.key == 'b':
        skyBelow = coords[len(coords)-2:]
        if rectSkyBelow is not None:
            rectSkyBelow.remove()
        rectSkyBelow = patches.Rectangle( ( xmin,skyBelow[0][1] ), xmax-xmin, skyBelow[1][1]-skyBelow[0][1], alpha = 0.5, ec = "gray", fc = "CornflowerBlue", visible = True)
        ax.add_patch(rectSkyBelow)
        fig.canvas.draw_idle()
    elif event.key == 'u':
        skyAbove = coords[len(coords)-2:]
        if rectSkyAbove is not None:
            rectSkyAbove.remove()
        rectSkyAbove = patches.Rectangle( ( xmin,skyAbove[0][1] ), xmax-xmin, skyAbove[1][1]-skyAbove[0][1], alpha = 0.5, ec = "gray", fc = "green", visible = True)
        ax.add_patch(rectSkyAbove)
        fig.canvas.draw_idle()
    elif event.key == 'o':
        obs = coords[len(coords)-2:]
        if rectObs is not None:
            rectObs.remove()
        rectObs = patches.Rectangle( ( xmin,obs[0][1] ), xmax-xmin, obs[1][1]-obs[0][1], alpha = 0.5, ec = "gray", fc = "orange", visible = True)
        ax.add_patch(rectObs)
        fig.canvas.draw_idle()
    coords = []
    print('skyBelow = ',skyBelow)
    print('skyAbove = ',skyAbove)
    print('obs = ',obs)

fName = 'image.fits'
image = getImageData(fName)

fig, ax = plt.subplots()
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

kpid = fig.canvas.mpl_connect('key_press_event', onKeyPress)
krid = fig.canvas.mpl_connect('key_release_event', onKeyRelease)

d_val=(vmax-vmin)*1

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)
sfreq = mwidgets.Slider(axfreq, 'vmax', vmax-d_val, vmax,valfmt='% .2f', valinit=0, valstep=0.01)
sfreq.on_changed(update_max)
sfreq.reset()
axfreq1 = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor=axcolor)
sfreq1 = mwidgets.Slider(axfreq1, 'vmin', vmax-d_val, vmax,valfmt='% .2f', valinit=0, valstep=0.01)
sfreq1.on_changed(update_min)
sfreq1.reset()
max_val=vmax
min_val=vmin
sfreq.set_val(vmax)
sfreq1.set_val(min_val)

plt.show()
