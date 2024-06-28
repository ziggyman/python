import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets
import matplotlib.patches as patches

import numpy as np
import os

from spectraUtils import getImageData

coords = []
skyBelow = []
skyAbove = []
obs = []
areaType = 'sky below'
axMainRect = [0.14,0.3,0.85,0.7]
rectSkyAbove = None
rectSkyBelow = None
rectObs = None
note = ''
fNameOut = ''
areasWritten = False
min_val = 0.
max_val = 0.

def defineRegions(folderName):
    global coords
    global skyAbove
    global skyBelow
    global obs
    global axMainRect
    global rectSkyAbove
    global rectSkyBelow
    global rectObs
    global note
    global fNameOut
    global areasWritten
    global min_val
    global max_val
    global areaType

    toExtractFName = os.path.join(folderName,'to_extract.csv')
    f = open(toExtractFName,'r')
    toExtract = f.readlines()
    toExtract = toExtract[1:]
    fNameOut = os.path.join(folderName,'areas.csv')

#    if os.path.exists(fNameOut):
#        os.remove(fNameOut)

    for line in toExtract:
        areasWritten = False
        fName = line[:line.find(',')]
        note = line[line.find(',')+1:]
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

        image = getImageData(fName)

        fig = plt.figure(figsize=(15,6))
        axMainRect = [0.14,0.3,0.85,0.7]
        axMain = plt.axes(axMainRect)
        axAreaType = plt.axes([0.01,0.55,0.08,0.13])
        axVMin = plt.axes([0.1,0.2,0.88,0.09])
        axVMax = plt.axes([0.1,0.1,0.88,0.09])
        axDone = plt.axes([0.01,0.9,0.1,0.1])
        axNote = plt.axes([0.1,0.01,0.88,0.08])
        #axSubmitNote = plt.axes([0.91,0.01,0.08,0.08])
        radio_background = 'lightgoldenrodyellow'

        axAreaType.set_facecolor(radio_background)
        #axSubmitNote.set_facecolor(radio_background)

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
                print('%s release: button=%d, x=%d, y=%d, xdata=%f, ydata=%f, areaType=%s' %
                    ('double' if event.dblclick else 'single', event.button,
                    event.x, event.y, event.xdata, event.ydata, areaType))
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
                        print('rectObs = ',rectObs,': obs = ',obs)
                        print('rectSkyBelow = ',rectSkyBelow,': skyBelow set to ',skyBelow)
                        print('rectSkyAbove = ',rectSkyAbove,': skyAbove = ',skyAbove)
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
                        print('rectObs = ',rectObs,': obs = ',obs)
                        print('rectSkyBelow = ',rectSkyBelow,': skyBelow = ',skyBelow)
                        print('rectSkyAbove = ',rectSkyAbove,': skyAbove set to ',skyAbove)
                        fig.canvas.draw_idle()
                    elif areaType == 'object':
                        obs = coords[len(coords)-2:]
                        obs.sort()
                        obs = [int(obs[0]),int(obs[1])+1]
                        if rectObs is not None:
                            rectObs.remove()
                        rectObs = patches.Rectangle( ( xmin,obs[0] ), xmax-xmin, obs[1]-obs[0], alpha = 0.5, ec = "gray", fc = "orange", visible = True)
                        print('rectObs = ',rectObs,': obs set to ',obs)
                        print('rectSkyBelow = ',rectSkyBelow,': skyBelow = ',skyBelow)
                        print('rectSkyAbove = ',rectSkyAbove,': skyAbove = ',skyAbove)
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

        imMin = np.min(image)
        print('imMin = ',imMin)
        imMax = np.max(image)
        print('imMax = ',imMax)
        imMean = np.mean(image)
        print('imMean = ',imMean)
        vmin = imMean - (imMean - imMin)/2.
        vmax = imMean + (imMax - imMean)/2.
        im = axMain.imshow(image, origin='lower', vmin=vmin, vmax=vmax)

        def update_max(val):
            max_val=val
        #    axMain.clear()
            im.set_clim(vmax=max_val)
            fig.canvas.draw_idle()

        def update_min(val):
            min_val=val
            im.set_clim(vmin=min_val)
            fig.canvas.draw_idle()


        radio = mwidgets.RadioButtons(axAreaType, ('sky below', 'sky above', 'object'), active=0)
        def setAreaType(label):
            global areaType
            areaType = label
        radio.on_clicked(setAreaType)

        def submit(text):
            global note
            note = text
            print('note set to <'+note+'>')
        noteBox = mwidgets.TextBox(axNote,'Notes',initial=note)
        noteBox.on_submit(submit)

        #d_val=(vmax-vmin)*1

        #axcolor = 'lightgoldenrodyellow'
        #axfreq = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)
        vMaxS = mwidgets.Slider(axVMax, 'vmax', vmin, vmax,valfmt='% .2f', valinit=vmax, valstep=(vmax-vmin)/1000.)
        vMaxS.on_changed(update_max)
        #vMaxS.reset()
        #axfreq1 = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor=axcolor)
        vMinS = mwidgets.Slider(axVMin, 'vmin', vmin, vmax,valfmt='% .2f', valinit=0, valstep=(vmax-vmin)/1000.)
        vMinS.on_changed(update_min)
        #vMinS.reset()
        max_val=vmax
        min_val=vmin
        vMaxS.set_val(vmax)
        vMinS.set_val(min_val)

        def done(event):
            global skyAbove
            global skyBelow
            global obs
            global fNameOut
            global areasWritten
            print('obs = ',obs,', skyBelow = ',skyBelow,', skyAbove = ',skyAbove)
            if not os.path.exists(fNameOut):
                with open(fNameOut,'w') as f:
                    f.write('fName,object,skyBelow,skyAbove,method,notes\n')
            if (len(obs) < 2) or (len(skyAbove) < 2) or (len(skyBelow) < 2):
                print('WARNING: NOT ALL # AREAS HAVE BEEN DEFINED!')
            else:
                with open(fNameOut,'a') as f:
                    f.write('%s,[%d:%d],[%d:%d],[%d:%d],sum,%s%s' % (fName,obs[0],obs[1],skyBelow[0],skyBelow[1],skyAbove[0],skyAbove[1],note,'\n' if '\n' not in note else ''))
                    print('%s,[%d:%d],[%d:%d],[%d:%d],sum,%s%s written to areas.csv' % (fName,obs[0],obs[1],skyBelow[0],skyBelow[1],skyAbove[0],skyAbove[1],note,'\n' if '\n' not in note else ''))
#            plt.close()
            areasWritten = True

        bDone = mwidgets.Button(axDone, 'Done')
        bDone.on_clicked(done)

        plt.show()
        if not areasWritten:
            if not os.path.exists(fNameOut):
                with open(fNameOut,'w') as f:
                    f.write('fName,object,skyBelow,skyAbove,method,notes\n')
            if (len(obs) < 2) or (len(skyAbove) < 2) or (len(skyBelow) < 2):
                print('WARNING: NOT ALL # AREAS HAVE BEEN DEFINED!')
            else:
                if not areasWritten:
                    with open(fNameOut,'a') as f:
                        f.write('%s,[%d:%d],[%d:%d],[%d:%d],sum,%s%s' % (fName,obs[0],obs[1],skyBelow[0],skyBelow[1],skyAbove[0],skyAbove[1],note,'\n' if '\n' not in note else ''))
                        print('%s,[%d:%d],[%d:%d],[%d:%d],sum,%s%s written to areas.csv' % (fName,obs[0],obs[1],skyBelow[0],skyBelow[1],skyAbove[0],skyAbove[1],note,'\n' if '\n' not in note else ''))
