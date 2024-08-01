def cleanImage(inputImage, outputImage):
    from matplotlib.widgets import AxesWidget, RadioButtons, Slider
    import matplotlib.colors as colors
    from scipy import interpolate# import RegularGridInterpolator

    global cleanType
    global image
    global cleanRange

    fig = plt.figure(figsize=(15,9))
#    axMainRect = [0.04,0.2,0.95,0.5]
    ax2DRect = plt.axes([0.04,0.2,0.95,0.8])
#    axMain = plt.axes(axMainRect)
    axTrimClean = plt.axes([0.01,0.01,0.08,0.1])
    axVMin = plt.axes([0.15,0.06,0.84,0.04])
    axVMax = plt.axes([0.15,0.01,0.84,0.04])
    max_val=0
    min_val=0


#    spec = getImageData(inputSpec1D,0)
#    wLen = getWavelengthArr(inputSpec1D,0)
    image = np.array(getImageData(inputImage,0))
    vmax = np.max(image)#np.max([1.5 * np.mean(image), 1.])
    vmin = np.min(image)#np.min([1.5 * np.mean(image), 1.])
    im = ax2DRect.imshow(image, origin='lower', norm=colors.SymLogNorm(linthresh=0.5, linscale=1,
                                              vmin=vmin, vmax=vmax, base=10))#,vmin=vmin, vmax=vmax)#,cmap='gist_rainbow')
    ax2DRect.set_title(inputImage[inputImage.rfind('/')+1:inputImage.rfind('.')])
    cleanRange = []
    cleanType = 'trim'

    radio = RadioButtons(axTrimClean, ('trim', 'clean'), active=0)
    def setCleanType(label):
        global cleanType
        cleanType = label

    radio.on_clicked(setCleanType)


    def update_max(val):
        max_val=val
    #    axMain.clear()
        im.set_clim(vmax=max_val)
        fig.canvas.draw_idle()

    def update_min(val):
        min_val=val
        im.set_clim(vmin=min_val)
        fig.canvas.draw_idle()


    d_val=(vmax-vmin)/1000

    #axcolor = 'lightgoldenrodyellow'
    #axfreq = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)
    vMaxS = Slider(axVMax, 'vmax', vmin, vmax,valfmt='% .2f', valinit=0, valstep=d_val)
    vMaxS.on_changed(update_max)
    vMaxS.reset()
    #axfreq1 = plt.axes([0.15, 0.1, 0.65, 0.03], facecolor=axcolor)
    vMinS = Slider(axVMin, 'vmin', vmin, vmax,valfmt='% .2f', valinit=0, valstep=d_val)
    vMinS.on_changed(update_min)
    vMinS.reset()
    max_val=vmax
    min_val=vmin
    vMaxS.set_val(vmax)
    vMinS.set_val(min_val)

    def onClick(event):
        global image
        global cleanRange

        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            ('double' if event.dblclick else 'single', event.button,
            event.x, event.y, event.xdata, event.ydata))
        if event.inaxes is ax2DRect:
            print('fig.canvas.toolbar.mode = ',fig.canvas.toolbar.mode)
            if fig.canvas.toolbar.mode == '':
                print('cleanRange = ',cleanRange)
#                if cleanType == 'trim':
#                    if event.xdata < wLen[int(len(wLen)/2.)]:
#                        idx = np.where(wLen > event.xdata)
#                    else:
#                        idx = np.where(wLen < event.xdata)
#                    wLen = wLen[idx]
#                    spec = spec[idx]
#                    axMain.plot(wLen,spec)
#                    fig.canvas.draw_idle()
#                elif cleanType == 'clean':
                if len(cleanRange) == 1:
                    xlim = ax2DRect.get_xlim()
                    ylim = ax2DRect.get_ylim()
                    cleanRange.append([int(event.xdata),int(event.ydata)])
                    xRange = np.arange(0,image.shape[1],1)
                    yRange = np.arange(0,image.shape[0],1)
                    idx = np.where((xRange >= cleanRange[0][0]) & (xRange <= cleanRange[1][0]))[0]
#                    print('clean: idx = ',idx)
                    idy = np.where((yRange >= cleanRange[0][1]) & (yRange <= cleanRange[1][1]))[0]
#                    print('clean: idy = ',idy)
                    for x in idx:
                        for y in idy:
                            image[y,x] = np.nan
#                            print('image[',y,',',x,'] = ',image[y,x])
#                    image[idx,idy] = np.nan
                    ax2DRect.cla()
#                    ax2DRect.imshow(image, origin='lower', norm=colors.SymLogNorm(linthresh=0.5, linscale=1,
#                                              vmin=vmin, vmax=vmax, base=10))#,vmin=vmin, vmax=vmax)#,cmap='gist_rainbow')
                    x = np.arange(0,image.shape[1])
                    y = np.arange(0,image.shape[0])
                    image = np.ma.masked_invalid(image)
                    xx, yy = np.meshgrid(x,y)
#                    print('xx = ',xx.shape,': ',xx)
#                    print('xx = ',yy.shape,': ',yy)
                    x1 = xx[~image.mask]
#                    print('x1 = ',x1)
                    y1 = yy[~image.mask]
#                    print('y1 = ',y1)
                    newImage = image[~image.mask]
#                    print('newImage = ',newImage.shape,': ',newImage)
#                    ax2DRect.imshow(newImage, origin='lower', norm=colors.SymLogNorm(linthresh=0.5, linscale=1,
#                                              vmin=vmin, vmax=vmax, base=10))#,vmin=vmin, vmax=vmax)#,cmap='gist_rainbow')
                    image = interpolate.griddata((x1,y1), newImage.ravel(),(xx,yy),method='linear')
#                    grid = np.zeros((len(idx)+2,len(idy)+2))
#                    grid[:,:] = np.nan
#                    print('len(idx) = ',len(idx),', len(idy) = ',len(idy),': grid.shape = ',grid.shape)
#                    grid[:,0] = image[idy[0]-1,idx[0]-1:idx[len(idx)-1]+2]
#                    print('grid = ',grid)
#                    grid[:,len(idy)+1] = image[idy[len(idy)-1]+1,idx[0]-1:idx[len(idx)-1]+2]
#                    print('grid = ',grid)
#                    grid[0,:] = image[idy[0]-1:idy[len(idy)-1]+2,idx[0]-1]
#                    print('grid = ',grid)
#                    grid[len(idx)+1,:] = image[idy[0]-1:idy[len(idy)-1]+2,idx[len(idx)-1]+1]
#                    print('grid = ',grid)
#                    image[idy[0]-1:idy[len(idy)-1]+2,idx[0]-1:idx[len(idx)-1]+2] = np.transpose(grid)
                    ax2DRect.imshow(image, origin='lower', norm=colors.SymLogNorm(linthresh=0.5, linscale=1,
                                              vmin=vmin, vmax=vmax, base=10))#,vmin=vmin, vmax=vmax)#,cmap='gist_rainbow')
                    ax2DRect.set_xlim(xlim)
                    ax2DRect.set_ylim(ylim)
                    #plt.show()
                    fig.canvas.draw_idle()
#                    x = np.arange(0,grid.shape[1])
#                    y =
#                    for x in np.arange(1,grid.shape[0]-1,1):
#                        for y in np.arange(1,grid.shape[1]-1,1):
#                            smallGrid = np.array([[grid[0,y],grid[len(idx)+1,y]],[grid[x,0],grid[x,len(idy)+1]]])
#                            print('x = ',x,', y = ',y,': smallGrid = ',smallGrid)
#                            interp = RegularGridInterpolator((np.array([0,1]),np.array([0,1])), smallGrid)
#                            grid[x,y] = interp(np.array([x/grid.shape[0],y/grid.shape[1]]))
#                            print('xx=',x/grid.shape[0],', yy=',y/grid.shape[1],': grid[',x,',',y,'] = ',grid[x,y])
#                    print('image.shape = ',image.shape)
#                    image[idy[0]:idy[len(idy)-1]+1,idx[0]:idx[len(idx)-1]+1] = np.transpose(grid[1:len(idx)+1,1:len(idy)+1])
#                    print('image[',idx[0],':',idx[len(idx)-1],',',idy[0],':',idy[len(idy)-1],'] = ',image[idx[0]:idx[len(idx)-1],idy[0]:idy[len(idy)-1]])
#                    ax2DRect.imshow(grid, origin='lower', norm=colors.SymLogNorm(linthresh=0.5, linscale=1,
#                                              vmin=vmin, vmax=vmax, base=10))#,vmin=vmin, vmax=vmax)#,cmap='gist_rainbow')
#                    ax2DRect.imshow(image, origin='lower', norm=colors.SymLogNorm(linthresh=0.5, linscale=1,
#                                              vmin=vmin, vmax=vmax, base=10))#,vmin=vmin, vmax=vmax)#,cmap='gist_rainbow')
#                    ax2DRect.set_xlim(xlim)
#                    ax2DRect.set_ylim(ylim)
                    #plt.show()
#                    fig.canvas.draw_idle()
                else:
                    cleanRange = [[int(event.xdata),int(event.ydata)]]
#                else:
#                    print('ERROR: cleanType <'+cleanType+'> not supported')

    cid = fig.canvas.mpl_connect('button_press_event', onClick)

    plt.show()
    writeFits(image,inputImage,outputImage)

