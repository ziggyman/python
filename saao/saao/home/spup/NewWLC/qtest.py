import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import warnings
from astropy.io import ascii,fits

from glob import glob
import json

def procChris(filename='a0061330.fits'):
    datadir='/home/ccd/data/20151204/'
    im2d,hdr=fits.getdata(datadir+filename,header=True)
    jfile = open('/home/ccd/.spup/qlsettings.json','r')
    foo = json.load(jfile)
    extractWin_x0 = foo['extractWin_x0']
    extractWin_x1 = foo['extractWin_x1']

    im1d = np.sum(im2d[extractWin_x0:extractWin_x1,:],0)
    good = np.arange(98,len(im1d)-3)
    return im1d[good],hdr['GR-ANGLE']

def wlcChris(cf,cang,pfit):
    x=np.arange(len(cf))
    yy = np.repeat(cang,len(x))
    lamfit = pfit(x,yy)
    return lamfit

files = glob('gr4*dc.txt')
plt.figure()
lam=[]
xpix=[]
order=[]
for f in files:
    angle = float(f.split('_')[1])
    x,y = np.loadtxt(f,unpack=True)
    xx=np.arange(len(x))
    plt.plot(xx,x,label=f)
    lam = np.append(lam,x)
    xpix = np.append(xpix,xx)
    order=np.append(order,np.repeat(angle,len(x)))
plt.legend()
plt.show()



z = lam
x = xpix
y = order


xord = 3; yord = 3

# Fit the data using astropy.modeling
p_init = models.Legendre2D(x_degree=xord,y_degree=yord)
fit_p = fitting.LevMarLSQFitter()

with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter
    warnings.simplefilter('ignore')
    #*** need to define correct equation here:
    #p = fit_p(p_init, x, y, z)
    p = fit_p(p_init, x, y, z)#*y)

yy = np.arange(order.min(),order.max(),0.2)
nord = len(yy)
xx = np.arange(0,2000,100)
nx = len(xx)
yyy = np.repeat(yy,nx)
xxx = np.tile(xx,nord)

#plt.figure()
#plt.plot(xxx,p(xxx,yyy)/yyy)
#plt.show()

# Use this to predict a different grating angle:

##ang=-3.5

# test new:
#ang=3.5 ; file1d = 'a0053713_1D.txt'
ang=-2.0 ; file1d = 'a0053711_1D.txt'




xx = np.arange(0,2048)
yy = np.repeat(ang,len(xx))
f = p(xx,yy)#/yy

plt.plot(xx,f,lw=2)
plt.show()

#stop()
# check residuals:for f in files:
for f in files:

#    if f=='gr4_2.00_dc.txt':
#        print 'skipping %s'%f
#        continue

    angle = float(f.split('_')[1])
    x,y = np.loadtxt(f,unpack=True)

    yy = np.repeat(angle,len(x))
    xx=np.arange(len(x))
    f = p(xx,yy)#/yy
    resid = f-x
    print angle,np.max(np.abs(resid))
    plt.plot(xx,f,'--')
plt.show()


a,b,c,d = np.loadtxt(file1d,unpack=True)
a = a*0. # make sure not to use these dummy wavelengths
xx=np.arange(len(x))
yy = np.repeat(ang,len(x))
f = p(xx,yy)#/yy

# Nice! seems to correctly predict intermediate solutions! (just doesn't work so well on extrapolations)
# Try some more***

ok=np.where(np.isnan(b)==False)[0]
#plt.plot(x,f,':')

plt.figure()
plt.plot(f,b[ok],'-',alpha=0.2,lw=3)
plt.show()

gdat = ascii.read('GCALcuar_short.dat') # This seems to be a list of weights(?) rather than our intensities
##wok = np.where( (gdat['col1']>np.min(a[ok])) & (gdat['col1']<np.max(a[ok])) )[0]
wok = np.arange(len(gdat)) #all
for i in range(len(wok)):
#    plt.axvline(gdat['col1'][wok[i]],gdat['col2'][wok[i]],color='gray')
     plt.plot(np.array([gdat['col1'][wok[i]],gdat['col1'][wok[i]]]),650.+0.1*np.array([0.,gdat['col2'][wok[i]]]),color='gray')

plt.show()

#4880 5220
#plt.plot(f+515.,b[ok],'-')
plt.plot(f+0.,b[ok],'-')
plt.show()


x,y=np.loadtxt('gr4_-2.00_dc.txt',unpack=True)

##plt.figure()
plt.plot(x,y)
wok = np.arange(len(gdat)) #all
for i in range(len(wok)):
#    plt.axvline(gdat['col1'][wok[i]],gdat['col2'][wok[i]],color='gray')
#     plt.plot(np.array([gdat['col1'][wok[i]],gdat['col1'][wok[i]]]),650.+0.1*np.array([0.,gdat['col2'][wok[i]]]),color='gray')
     plt.plot(np.array([gdat['col1'][wok[i]],gdat['col1'][wok[i]]]),650.+0.1*np.array([0.,gdat['col2'][wok[i]]]),color='red',lw=2)

plt.show()

##cf,cang = procChris()
#junk,cf = np.loadtxt('a0053715_1D.txt',unpack=True,usecols=(0,1)) ; cang=6.00
junk,cf = np.loadtxt('a0053713_1D.txt',unpack=True,usecols=(0,1)) ; cang=2.00
cang=float(cang)
clam = wlcChris(cf,cang,p)
plt.figure()
#plt.plot(clam,cf-10000.)
plt.plot(clam,cf)

wok = np.where( (gdat['col1']>np.min(clam[ok])) & (gdat['col1']<np.max(clam[ok])) )[0]
for i in range(len(wok)):
     plt.plot(np.array([gdat['col1'][wok[i]],gdat['col1'][wok[i]]]),650.+0.1*np.array([0.,gdat['col2'][wok[i]]]),color='red',lw=2)


plt.show()
