import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets
import matplotlib.patches as patches

import numpy as np
import os

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
areaType = 'sky below'

fNameOut = 'areas.csv'
if os.path.exists(fNameOut):
    os.remove(fNameOut)

fName = 'image.fits'
image = getImageData(fName)

fig = plt.figure(figsize=(15,6))
axMainRect = [0.14,0.3,0.85,0.7]
axMain = plt.axes(axMainRect)
axAreaType = plt.axes([0.01,0.55,0.08,0.13])
axVMin = plt.axes([0.1,0.2,0.88,0.09])
axVMax = plt.axes([0.1,0.1,0.88,0.09])
axDone = plt.axes([0.01,0.9,0.1,0.1])
radio_background = 'lightgoldenrodyellow'

axAreaType.set_facecolor(radio_background)

max_val=0
min_val=0
def onClick(event):
    global coords
    global axMainRect
    if event.inaxes is axMain:
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))
        coords.append(event.ydata)
def onRelease(event):
    global coords
    global axMainRect
    global rectSkyAbove
    global rectSkyBelow
    global rectObs
    global areaType
    global obs
    global skyBelow
    global skyAbove
    if event.inaxes is axMain:
        print('%s release: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))
        print('fig.canvas.toolbar.mode = <',fig.canvas.toolbar.mode,'>')
        if fig.canvas.toolbar.mode == '':
            coords.append(event.ydata)
            xmin, xmax = axMain.get_xlim()
            print('xmin = ',xmin,', xmax = ',xmax)
            if areaType == 'sky below':
                skyBelow = coords[len(coords)-2:]
                skyBelow.sort()
                skyBelow = [int(skyBelow[0]),int(skyBelow[1])+1]
                print('skyBelow set to ',skyBelow)
                if rectSkyBelow is not None:
                    rectSkyBelow.remove()
                rectSkyBelow = patches.Rectangle( ( xmin,skyBelow[0] ), xmax-xmin, skyBelow[1]-skyBelow[0], alpha = 0.5, ec = "gray", fc = "CornflowerBlue", visible = True)
                axMain.add_patch(rectSkyBelow)
                fig.canvas.draw_idle()
            elif areaType == 'sky above':
                skyAbove = coords[len(coords)-2:]
                skyAbove.sort()
                skyAbove = [int(skyAbove[0]),int(skyAbove[1])+1]
                print('skyAbove set to ',skyAbove)
                if rectSkyAbove is not None:
                    rectSkyAbove.remove()
                rectSkyAbove = patches.Rectangle( ( xmin,skyAbove[0] ), xmax-xmin, skyAbove[1]-skyAbove[0], alpha = 0.5, ec = "gray", fc = "green", visible = True)
                axMain.add_patch(rectSkyAbove)
                fig.canvas.draw_idle()
            elif areaType == 'object':
                obs = coords[len(coords)-2:]
                obs.sort()
                obs = [int(obs[0]),int(obs[1])+1]
                print('obs set to ',obs)
                if rectObs is not None:
                    rectObs.remove()
                rectObs = patches.Rectangle( ( xmin,obs[0] ), xmax-xmin, obs[1]-obs[0], alpha = 0.5, ec = "gray", fc = "orange", visible = True)
                axMain.add_patch(rectObs)
                fig.canvas.draw_idle()
            else:
                print('WARNING: areaType(=',areaType,') not recognised')
            coords = []
        #    print('skyBelow = ',skyBelow)
        #    print('skyAbove = ',skyAbove)
        #    print('obs = ',obs)

cid = fig.canvas.mpl_connect('button_press_event', onClick)
rid = fig.canvas.mpl_connect('button_release_event', onRelease)
print('coords = ',coords)


vmax = np.max([1.5 * np.mean(image), 1.])
vmin=np.min([1.5 * np.mean(image), 1.])
im = axMain.imshow(image, origin='lower', vmin=vmin, vmax=vmax)

def update_max(val):
    print(val)
    max_val=val
#    axMain.clear()
    im.set_clim(vmax=max_val)
    fig.canvas.draw_idle()

def update_min(val):
    print(val)
    min_val=val
    im.set_clim(vmin=min_val)
    fig.canvas.draw_idle()


radio = mwidgets.RadioButtons(axAreaType, ('sky below', 'sky above', 'object'), active=0)
def setAreaType(label):
    global areaType
    areaType = label
radio.on_clicked(setAreaType)

d_val=(vmax-vmin)*1

#axcolor = 'lightgoldenrodyellow'
#axfreq = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)
vMaxS = mwidgets.Slider(axVMax, 'vmax', vmax-d_val, vmax,valfmt='% .2f', valinit=0, valstep=0.01)
vMaxS.on_changed(update_max)
vMaxS.reset()
#axfreq1 = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor=axcolor)
vMinS = mwidgets.Slider(axVMin, 'vmin', vmax-d_val, vmax,valfmt='% .2f', valinit=0, valstep=0.01)
vMinS.on_changed(update_min)
vMinS.reset()
max_val=vmax
min_val=vmin
vMaxS.set_val(vmax)
vMinS.set_val(min_val)


def done(event):
    global skyAbove
    global skyBelow
    global obs
    global fNameOut
    print('obs = ',obs,', skyBelow = ',skyBelow,', skyAbove = ',skyAbove)
    if not os.path.exists(fNameOut):
        with open(fNameOut,'w') as f:
            f.write('fName,object,skyBelow,skyAbove,method,notes\n')
    with open(fNameOut,'a') as f:
        f.write('%s,[%d:%d],[%d:%d],[%d:%d],sum,\n' % (fName,obs[0],obs[1],skyBelow[0],skyBelow[1],skyAbove[0],skyAbove[1]))

bDone = mwidgets.Button(axDone, 'Done')
bDone.on_clicked(done)

plt.show()
