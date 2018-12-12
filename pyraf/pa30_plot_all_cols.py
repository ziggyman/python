import os
import pyfits
import pyraf
from pyraf import iraf

#import gi
#gi.require_version('Gtk', '3.0')
#from gi.repository import Gtk

#from matplotlib.figure import Figure
#from matplotlib.backends.backend_gtk3agg import FigureCanvas
#from matplotlib.backends.backend_gtk3 import (
#    NavigationToolbar2GTK3 as NavigationToolbar)

execfile("/Users/azuri/entwicklung/python/myUtils.py")# import getDate, findClosestDate,...

#win = Gtk.Window()
#win.connect("destroy", lambda x: Gtk.main_quit())
#win.set_default_size(400,300)
#win.set_title("Embedding in GTK")

#vbox = Gtk.VBox()
#win.add(vbox)

#fig = Figure(figsize=(5,4), dpi=100)
#ax = fig.add_subplot(111)

wcsFitsName = '/Users/azuri/daten/uni/HKU/Pa30/gtc_object_wcs.fits'
center = [1469.676, 994.]
plotCols = np.arange(1323,1709)#1401)

tab=[{'FName': '/Users/azuri/daten/uni/HKU/Pa30/Pa30_combined.fits', 'ylim' : [0,0], 'OName': 'Pa30', 'SkyLeft':[0,0], 'SkyRight':[0,0], 'ObjectAreas':[[1323,1341],[1387,1401],[1434,1490],[1479,1561],[1650,1709]]},
    ]
obs = tab[0]
directory = obs['FName'][0:obs['FName'].rfind('/')]
print('directory = <'+directory+'>')
print(dir(iraf))
#    iraf.user.reduce_gtcmos(indirs=directory)
image_file = obs['FName']
try:
    hdulist = pyfits.open(image_file)
except:
    print('file <'+image_file+'> not found')
    image_file = image_file[0:image_file.rfind('_x')]+image_file[image_file.rfind('_x')+2:]
    print('trying <'+image_file+'>')
    try:
        hdulist = pyfits.open(image_file)
    except:
        print('file <'+image_file+'> not found')

#    print 'type(hdulist) = ',type(hdulist)
header = hdulist[0].header
if header['OBJECT'] != '':
    obs['OName'] = header['OBJECT']

wavelength = getWavelength(header)
image_data = fits.getdata(image_file)
maxSpec = 0.
yShift0 = getArcsecDistance(wcsFitsName, center[0], center[1], center[0]+1, center[1])
for area in tab[0]['ObjectAreas']:
    print('area = ',area)
    for column in np.arange(area[0], area[1]):
        spectrum = image_data[:,column]
        if maxSpec == 0.0:
            maxSpec = np.max(spectrum) / 10.
        print('max(spectrum) = ',np.max(spectrum))
        yShift = getArcsecDistance(wcsFitsName, center[0], center[1], column, center[1])# * maxSpec / yShift0
        if column < center[0]:
            yShift = 0. - yShift
        print('column = ',column,': yShift = ',yShift)
        specOutName = image_file[0:image_file.rfind('/')+1]+'Pa30_col'+str(column)+'_yShift'+str(yShift)+'.fits'
        print('specOutName = <'+specOutName+'>')
        hdulist[0].data = spectrum
#        hdulist[0].header['NAXIS1'] = header['NAXIS2']
        hdulist[0].header['CRPIX1'] = header['CRPIX2']
        hdulist[0].header['CRVAL1'] = wavelength[0]
        hdulist[0].header['CDELT1'] = header['CDELT2']
        hdulist[0].header['CD1_1'] = header['CD2_2']
        hdulist[0].header['WAT1_001'] = header['WAT2_001']
    #        hdulist[0].header[''] = header['']
        obsdate = getDateTime(header['DATE-OBS'])
        datestr = obsdate.strftime('%d%m%y')
        print 'datestr = <'+datestr+'>'

        hdulist.writeto(specOutName, clobber=True)

        yShift = yShift * maxSpec / yShift0
        if column in plotCols:
            plt.plot(wavelength, spectrum + yShift)
plt.xlim(6600., 6800.)
plt.ylim(-7.0e-17,9.3e-17)
plt.show()
#canvas = FigureCanvas(fig)  # a Gtk.DrawingArea
#vbox.pack_start(canvas, True, True, 0)
#toolbar = NavigationToolbar(canvas, win)
#vbox.pack_start(toolbar, False, False, 0)

#win.show_all()
#Gtk.main()
