
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import warnings
from astropy.io import ascii,fits

from glob import glob
import json,cPickle


files = glob('gr4*dc.txt')
plt.figure()
for i,f in enumerate(np.sort(files)):
    angle = float(f.split('_')[1])
    x,y = np.loadtxt(f,unpack=True)
    xx=np.arange(len(x))

    if i==0:
        aangle=np.copy(angle)
        alam=np.copy(x)
        aflux=np.copy(y)
    else:
        aangle=np.append(aangle,angle)
        alam=np.vstack([alam,x])
        aflux=np.vstack([aflux,y])
    plt.plot(xx,x,label=f)
    plt.plot(xx,x-np.median(x))
    coeff=np.polyfit(xx,x,1)
#    print angle, coeff
    yfit=np.polyval(coeff,np.arange(len(xx)))
    plt.plot(xx,x-yfit)
    coeff3 = np.polyfit(xx,x-yfit,3)
    print angle, coeff3

plt.legend()
plt.show()

#plot(alam.T,aflux.T)

# how well does it work to approximate fit by interpolated central wavelength and dispersion?
ss=np.argsort(aangle)
aangle=aangle[ss]
alam=alam[ss]
aflux=aflux[ss]
plt.figure()
plt.plot(aangle,np.median(alam,1),'-o')
plt.show()
plt.figure()
plt.plot(aangle,(np.max(alam,1)-np.min(alam,1))/2048.,'-o')
plt.plot(aangle,(np.max(alam[:,900:1100],1)-np.min(alam[:,900:1100],1))/200.,'-o') # ave. central dispersion
plt.show()

# linear approx.
# let's try to see how close we get to WLC by using linear approx.
if(0):
    dispns = (np.max(alam,1)-np.min(alam,1))/2048.
    cenlams = np.median(alam,1)

    angle=-6.
    ok=np.where(aangle == angle)[0][0]
    xxx = np.arange(2047)
    lamlinfit =  dispns[ok]*xxx 
    lamlinfit -= np.median(lamlinfit)
    lamlinfit +=cenlams[ok]

    plt.figure()
    plt.plot(alam[ok,:],aflux[ok,:])
    plt.plot(lamlinfit,aflux[ok,:])
    plt.show()


#linear w. average cubic correction:
if(0):
    coeff3 = np.array([ -6.90984694e-09 , -4.43208745e-06 ,  3.51051695e-02 , -1.49240823e+01]) # from ang=+2.0

    dispns = (np.max(alam,1)-np.min(alam,1))/2048.
    cenlams = np.median(alam,1)

    angle=-6.
    ok=np.where(aangle == angle)[0][0]
    xxx = np.arange(2047)
    lamlinfit =  dispns[ok]*xxx 
    lamlinfit -= np.median(lamlinfit)
    lamlinfit +=cenlams[ok]
    
    poly3 = np.polyval(coeff3,xxx)
    lamave3fit = lamlinfit + poly3 #-np.median(poly3)

    plt.figure()
    plt.plot(alam[ok,:],aflux[ok,:],label='ref')
    plt.plot(lamlinfit,aflux[ok,:],label='linear')
    plt.plot(lamave3fit,aflux[ok,:],label='linear + ave. cubic')
    plt.legend()
    plt.show()

# Now try full 2D fit. Is this any better than above?



# Fit the data using astropy.modeling
#xord = 3; yord = 3
#p_init = models.Legendre2D(x_degree=xord,y_degree=yord)
degree=3
p_init = models.Polynomial2D(degree=degree)
fit_p = fitting.LevMarLSQFitter()

# want to fit wavelength as a function of x pixel and angle

with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter
    warnings.simplefilter('ignore')
    #*** need to define correct equation here:
    xxx = np.arange(2047)

    angle=4.00

    # skip the angle we're trying to predict in the fitting procedure!:
    ##usefit = np.where(aangle <> angle)[0]
    usefit = np.arange(len(aangle))
    eaangle = aangle[usefit]
    ealam = alam[usefit,:] 
###    p2D = fit_p(p_init, np.tile(xxx,len(aangle)), np.repeat(aangle,len(xxx)), alam.flatten())
    p2D = fit_p(p_init, np.tile(xxx,len(eaangle)), np.repeat(eaangle,len(xxx)), ealam.flatten())

    
    lamp2Dfit = p2D(xxx,np.repeat(angle,len(xxx)))
    ok=np.where(aangle == angle)[0][0]
    plt.figure()
    plt.plot(alam[ok,:],aflux[ok,:],label='ref')
#    plt.plot(lamlinfit,aflux[ok,:],label='linear')
    plt.plot(lamp2Dfit,aflux[ok,:],label='2D poly fit')
    plt.legend()
    plt.show()

    # okay! 2D poly looks much better!!!

    # == let's store these coeffs and move on to trying to fit data from another night ==

    # hmm, should prob also store which are good pixels from orig spec:
    aOrigxPix = np.arange(98,2145,1) # yep. c.f.:
    # tag = -np.isnan(foo[:,1])
    # ok=np.where(tag)[0]
    # aOrigxPix = ok


    # store measurements/WLC solns:
    np.savez('gr4_CuAr_WLC.npz', aangle=aangle, alam=alam, aflux=aflux, aOrigxPix=aOrigxPix)
    
    # pickle coeffs (and fitting function):
    pf = open('gr4_CuAr.pkl','wb')
    cPickle.dump(p2D,pf)
    pf.close()

