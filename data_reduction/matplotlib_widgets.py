import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets
import numpy as np

cid = None
rid = None
coords = []
skyBelow = []
skyAbove = []
obs = []

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
    print('You released the %s button' % (event.key))
    fig.canvas.mpl_disconnect(cid)
    fig.canvas.mpl_disconnect(rid)
    cid = None
    rid = None
    if event.key == 'l':
        skyBelow = coords[len(coords)-2:]
    elif event.key == 'u':
        skyAbove = coords[len(coords)-2:]
    elif event.key == 'o':
        obs = coords[len(coords)-2:]
    coords = []
    print('skyBelow = ',skyBelow)
    print('skyAbove = ',skyAbove)
    print('obs = ',obs)

fig, ax = plt.subplots()
ax.plot(np.random.rand(10))

kpid = fig.canvas.mpl_connect('key_press_event', onKeyPress)
krid = fig.canvas.mpl_connect('key_release_event', onKeyRelease)

plt.show()
