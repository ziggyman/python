# Attempt at a new, from-scratch WLC fitting package for SpUpNIC
# Use combination of predicted cal (empirical model) and blind fitting

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits,ascii
from astropy.modeling import models, fitting
import warnings,cPickle
from scipy.interpolate import interp1d

#from robust_linefit import robust_linefit
#from numpy import linalg
#import scikits.statsmodels.api as sm
#import matplotlib
#from scipy.stats import norm
#import matplotlib.mlab as mlab
#from scipy.optimize import curve_fit
#from scipy.interpolate import interp1d

from glob import glob
import json

#####################################################################
# just gr4 for now

grating = 'gr4'
lamp = 'CuAr'


# load libraries for relevant setup
f = np.load('%s_%s_WLC.npz'%(grating,lamp))
# just in case we want to compare with original (by-hand) fits


pf = open('%s_%s.pkl'%(grating,lamp),'rb')
p2D = cPickle.load(pf)
pf.close()







def xcor3(inspec,refspec,maxshift,dx=1.,quiet='false'):#       tshift=xcor(tflux,tay0,tay2,maxshift)
    # Apply quadratic interpolation to the CCF around the peak estimated from the discrete CC.
    # This is equivalent to upsampling the input data.
    # http://dsp.stackexchange.com/questions/2321/is-up-sampling-prior-to-cross-correlation-useless

    # make running vector from -maxshift,+maxhsift
    ivec=np.arange(-maxshift,maxshift)

    refvec=refspec

    xcorval=np.zeros(len(ivec))
    for ii in range(len(ivec)):
        shifrefvec=np.roll(refvec,int(ivec[ii]))

#        print np.shape(shifrefvec)
#        print np.shape(inspec)
        xcorval[ii]=np.sum(shifrefvec * inspec)
#        if(quiet=='false'):
#            plt.plot(ivec[ii],xcorval[ii]/np.sum(inspec)/np.sum(refvec)*10000.,'ok')
#            print ii,ivec[ii],xcorval[ii]#/np.sum(inspec)/np.sum(refvec)

    # interpolate here:
    # ivec, xcorval
    fn = interp1d(ivec,xcorval, kind='quadratic')

    iivec = np.arange(-maxshift,maxshift,dx)
    ixcorval = fn(iivec)

    pk=np.argmax(ixcorval)
#    stop()
    if(quiet=='false'):
        if ((pk==0) | (pk==(len(ixcorval)-1))):
            print
            print 'WARNING: xcor2 peak at edge of window!'
            print 'CHECK RESULTS'
            print
    return iivec[pk],ixcorval[pk]

# load data to calibrate:

def procChris(filename='a0061330.fits'):
    datadir='/home/ccd/data/20151204/'
    im2d,hdr=fits.getdata(datadir+filename,header=True)
    jfile = open('/home/ccd/.spup/qlsettings.json','r')
    foo = json.load(jfile)
    extractWin_x0 = foo['extractWin_x0']
    extractWin_x1 = foo['extractWin_x1']

    im1d = np.sum(im2d[extractWin_x0:extractWin_x1,:],0)
    im1d=im1d[::-1] # reverse!
    good = np.arange(98,len(im1d)-3)
    return im1d[good],float(hdr['GR-ANGLE'])

def wlcChris(cf,cang,pfit):
    x=np.arange(len(cf))
    yy = np.repeat(cang,len(x))
    lamfit = pfit(x,yy)
    return lamfit

def buildRefSpec(lamp2Dfit,f):
    # XXX this interpolating overlapping spectra leads to weird effects! prob shouldn't use

    # Maybe this is fine (use some smoothing?). Gives same results as finding nearest unique spec

    # Construct reference spectrum by taking full wavelength library arcs, combining, and interpolating to same pixels as predicted WLC for real data:
#    for i in range(f.f.alam.shape[0]):
#        tirefflux = np.interp(lamp2Dfit,f.f.alam[i,:],f.f.aflux[i,:])
    flam = f.f.alam.flatten()
#    fflam = np.interp1d(

    fflux = f.f.aflux.flatten()
    ss=np.argsort(flam)
    flam=flam[ss]
    fflux=fflux[ss]

    tirefflux = np.interp(lamp2Dfit,flam,fflux)
    return tirefflux

def findRefSpec(lamp2Dfit,f):
    # find which spec has best overlap with current data
    besti=np.NaN
    bestlen=0
    for i in range(f.f.alam.shape[0]):
#        tirefflux = np.interp(lamp2Dfit,f.f.alam[i,:],fflux)
        olap = np.where( (f.f.alam[i,:]>lamp2Dfit.min()) & (f.f.alam[i,:]<=lamp2Dfit.max()) )[0]
        print i,len(olap)
        if (len(olap)>bestlen):
            besti = np.copy(i)
            bestlen = np.copy(len(olap))
    print besti
    tirefflux = np.interp(lamp2Dfit,f.f.alam[besti,:],f.f.aflux[besti,:])
    return tirefflux

    


#=====================================================================================
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Currently returns two lists of tuples, but maybe arrays would be better
    """
    maxtab = []
    mintab = []

    if x is None:
        x = np.arange(len(v))

    v = np.asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN

    lookformax = True
    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return maxtab, mintab
#=====================================================================================
def peakdetection(y, thres, min_dist):
    '''Peak detection routine.

    Finds the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks.

    Parameters
    ----------
    y : ndarray
        1D amplitude data to search for peaks.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.

    Returns
    -------
    ndarray
        Array containing the indexes of the peaks that were detected
    '''
    thres *= np.max(y) - np.min(y)

    # find the peaks by using the first order difference
    dy = np.diff(y)
    peaks = np.where((np.hstack([dy, 0.]) < 0.)
                     & (np.hstack([0., dy]) > 0.)
                     & (y > thres))[0]

    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak-min_dist), peak+min_dist+1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]

    return peaks


#=====================================================================================
def polyfitr(x, y, order, clip, xlim=None, ylim=None, mask=None, debug=False):
    """ Fit a polynomial to data, rejecting outliers.

    Fits a polynomial f(x) to data, x,y.  Finds standard deviation of
    y - f(x) and removes points that differ from f(x) by more than
    clip*stddev, then refits.  This repeats until no points are
    removed.

    Inputs
    ------
    x,y:
        Data points to be fitted.  They must have the same length.
    order: int (2)
        Order of polynomial to be fitted.
    clip: float (6)
        After each iteration data further than this many standard
        deviations away from the fit will be discarded.
    xlim: tuple of maximum and minimum x values, optional
        Data outside these x limits will not be used in the fit.
    ylim: tuple of maximum and minimum y values, optional
        As for xlim, but for y data.
    mask: sequence of pairs, optional
        A list of minimum and maximum x values (e.g. [(3, 4), (8, 9)])
        giving regions to be excluded from the fit.
    debug: boolean, default False
        If True, plots the fit at each iteration in matplotlib.

    Returns
    -------
    coeff, x, y:
        x, y are the data points contributing to the final fit. coeff
        gives the coefficients of the final polynomial fit (use
        np.polyval(coeff,x)).

    Examples
    --------
    >>> x = np.linspace(0,4)
    >>> np.random.seed(13)
    >>> y = x**2 + np.random.randn(50)
    >>> coeff, x1, y1 = polyfitr(x, y)
    >>> np.allclose(coeff, [1.05228393, -0.31855442, 0.4957111])
    True
    >>> coeff, x1, y1 = polyfitr(x, y, order=1, xlim=(0.5,3.5), ylim=(1,10))
    >>> np.allclose(coeff, [3.23959627, -1.81635911])
    True
    >>> coeff, x1, y1 = polyfitr(x, y, mask=[(1, 2), (3, 3.5)])
    >>> np.allclose(coeff, [1.08044631, -0.37032771, 0.42847982])
    True
    """

    x = np.asanyarray(x)
    y = np.asanyarray(y)
    isort = x.argsort()
    x, y = x[isort], y[isort]

    keep = np.ones(len(x), bool)
    if xlim is not None:
        keep &= (xlim[0] < x) & (x < xlim[1])
    if ylim is not None:
        keep &= (ylim[0] < y) & (y < ylim[1])
    if mask is not None:
        badpts = np.zeros(len(x), bool)
        for x0,x1 in mask:
            badpts |=  (x0 < x) & (x < x1)
        keep &= ~badpts

    x,y = x[keep], y[keep]
    if debug:
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.plot(x,y,'.')
        ax.set_autoscale_on(0)
        pl.show()

    coeff = np.polyfit(x, y, order)
    if debug:
        pts, = ax.plot(x, y, '.')
        poly, = ax.plot(x, np.polyval(coeff, x), lw=2)
        pl.show()
        raw_input('Enter to continue')
    norm = np.abs(y - np.polyval(coeff, x))
    stdev = np.std(norm)
    condition =  norm < clip * stdev
    y = y[condition]
    x = x[condition]
    while norm.max() > clip * stdev:
        if len(y) < order + 1:
            raise Exception('Too few points left to fit!')
        coeff = np.polyfit(x, y, order)
        if debug:
            pts.set_data(x, y)
            poly.set_data(x, np.polyval(coeff, x))
            pl.show()
            raw_input('Enter to continue')
        norm = np.abs(y - np.polyval(coeff, x))
        stdev = norm.std()
        condition =  norm < clip * stdev
        y = y[condition]
        x = x[condition]

    return coeff,x,y
#----------------------------------------------------------------





im1d, grang = procChris() # a0061330.fits    
im1d, grang = procChris(filename='a0061303.fits')
print len(im1d)
# should really check length against library file***

# get grating angle (and also check for config: grating, lamp***):
xxx = np.arange(f.f.alam.shape[1])

lamp2Dfit = p2D(xxx,np.repeat(grang,len(xxx)))

plt.figure()
plt.plot(lamp2Dfit,im1d/100.+200.,alpha=0.3,label='data')
#plt.plot(lamp2Dfit+35,im1d/100.+300.)


    
#*** should we apply some padding and go beyond end of data with template??
refspec = buildRefSpec(lamp2Dfit,f)
##refspec = findRefSpec(lamp2Dfit,f)
# This isn't actually used in finding solution:
"""
offset, rval = xcor3(im1d,refspec,300.,dx=1.0)
print offset
plt.plot(lamp2Dfit-offset,im1d/100.+300.,label='best offset')
# This seems ~okay for finding very crude offset
"""
plt.plot(lamp2Dfit,refspec-100,'r-',label='ref')
plt.legend()
plt.show()

#for i in range(f.f.alam.shape[0]):
#    plt.plot(f.f.alam[i,:],f.f.aflux[i,:],'-',alpha=0.3)
#plt.show()


# -- Tunable pars:
# peak detn:
fracThresh = 0.3#5
minDis = 50.
# peak binning:
maxOffset = 500. #pix
binSize = 20. # pix
# --



# ---- Find peaks:
#a, b = peakdet(im1d,100)
xIndPeaks = peakdetection(im1d,fracThresh,minDis)
nPeaks = len(xIndPeaks)
# sort in order of decreasing amplitude:
ss = np.argsort(im1d[xIndPeaks])[::-1]
xIndPeaks = xIndPeaks[ss]


xRefIndPeaks = peakdetection(refspec,fracThresh,minDis)
nRefPeaks = len(xRefIndPeaks)
# sort in order of decreasing amplitude:
sss = np.argsort(refspec[xRefIndPeaks])[::-1]
xRefIndPeaks = xRefIndPeaks[sss]


# build array of reference hist distances:
bins = np.arange(-1.*maxOffset,1.0*maxOffset,binSize)
nbins = len(bins)-1
refHists = np.zeros( (nRefPeaks,nbins) )
for j in range(nRefPeaks):
    tRefHist = np.histogram( xRefIndPeaks - xRefIndPeaks[j], bins=bins )
    refHists[j,:] = tRefHist[0]


# build histograms of distances for each line:
# work in pixels
x1=np.zeros(nPeaks)
x2=np.zeros(nPeaks)
for i in range(nPeaks):
    tHist = np.histogram( xIndPeaks - xIndPeaks[i], bins=bins)
#    print tHist[0]
    
    plt.axvline(lamp2Dfit[xIndPeaks[i]],color='b')
    # find nearest match in reference spec:
    diffs = np.ones(nRefPeaks)*99.0
    for j in range(nRefPeaks):
        tdiff = np.sum(np.abs(tHist[0]-refHists[j,:]))
        diffs[j]=tdiff
    bestMatch = np.argmin(diffs)
    x1[i] = xIndPeaks[i]
    x2[i] = xRefIndPeaks[bestMatch]
    plt.axvline(lamp2Dfit[xRefIndPeaks[bestMatch]],color='r')
                      
                            

#    if i==3: stop()

plt.figure()
plt.plot(x1,x2-x1,'ok')
plt.show()
medOff = np.median(x2-x1)
good = np.where(np.abs(x2-x1 - medOff) < 10.)[0]
plt.plot(x1[good],x2[good]-x1[good],'or')
plt.show()

order=3
clip=3.0
#coeffs,x,y = polyfitr(x1[good],x2[good]-x1[good],order,clip)
coeffs,x,y = polyfitr(x1,x2-x1,order,clip)
#
xxx=np.arange(len(lamp2Dfit))
plt.plot(xxx,np.polyval(coeffs,xxx))
plt.show()


# apply this correction to WLC:
newInd = xxx.astype('int')+np.polyval(coeffs,xxx).astype('int')
bad = np.where(newInd <0)[0]
newInd[bad]=0
bad = np.where(newInd> len(newInd)-1)[0]
newInd[bad]=len(newInd)-1
lamcal = lamp2Dfit[newInd]
# this isn't strictly monotonic. Should probably refit function to do this properly***

# Looks good, though!

# This may only work for small shifts against reference solutions? Need to check***
# *** should probably also fit full solution so can extrapolate outside overlap region

plt.figure()
plt.plot(lamp2Dfit,refspec-100,color='r')
plt.plot(lamcal,im1d/100.+400.)
plt.show()

# Should we have another iteration of centroiding against a line list?***


