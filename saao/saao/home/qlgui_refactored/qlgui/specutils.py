
import numpy as np
import pylab as plt
#import pyfits,sys,glob,re
from astropy.io import fits as pyfits

#from robust_linefit import robust_linefit
from numpy import linalg
import scikits.statsmodels.api as sm
import matplotlib
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

#idebug=0




def robust_linefit(x,y):
        
        import scikits.statsmodels.api as sm
        from scipy import polyval
        import numpy as np

        X=np.c_[x,np.ones(len(x))]
        y2=y
        # resrlm2 = sm.RLM(y2, X,).fit()
        resrlm2 = sm.RLM(y2, X,M=sm.robust.norms.TukeyBiweight()).fit()
        coeff=resrlm2.params

        return coeff

#xxx=np.arange(0.1,0.6,0.1)
#yfr=polyval(coeff,xxx)
#pl.plot(xxx,yfr,'-g')


#from numpy.polynomial.polynomial import polyvander
#from numpy.polynomial.polynomial import polyvander2d

"""
;    Xi and Yi are expressed as polynomials of Xo, Yo:
;        Xi = Kx[i,j] * Xo^j * Yo^i   Summed for i,j = 0 to degree.
;    And
;        Yi = Ky[i,j] * Xo^j * Yo^i.
;
;    This coordinate transformation may be then used to
;    map from Xo, Yo coordinates into Xi, Yi coordinates.
"""

# NOTE: USE:
# xf,yf=applywarp(xo,yo,kx,ky)

#** This hasn't been exhaustively tested, but i think it's ok **
def applywarp(xo,yo,kx,ky):
    np1=np.shape(kx)[0]
    n=np1#-1
    m=np.size(xo)
    xi=np.zeros(m)
    yi=np.zeros(m)
    for i in range(n):
        for j in range(n):
            xi=xi+kx[j,i]*(xo**j)*(yo**i)
            #print xi
            yi=yi+ky[j,i]*xo**j*yo**i
    # this python indexing is really confusing!
    # Seems to work now, though!
    # xi and yi are reversed
    
    return xi,yi
#    return yi,xi

#def polyvander2d(x,y,deg):
#    # from future version of numpy
#    ideg = [int(d) for d in deg]
#    is_valid = [id == d and id >= 0 for id, d in zip(ideg, deg)]
#    if is_valid != [1, 1]:
#        raise ValueError("degrees must be non-negative integers")
#    degx, degy = ideg
#    x, y = np.array((x, y), copy=0) + 0.0
#
#    vx = polyvander(x, degx)
#    vy = polyvander(y, degy)
#    v = vx[..., None]*vy[..., None, :]
#    # einsum bug
#    #v = np.einsum("...i,...j->...ij", vx, vy)
#    return v.reshape(v.shape[:-2] + (-1,))



def polywarp(xi,yi,xo,yo, degree):
#if (1):
#    xi=np.arange(10)+5
#    yi=np.arange(10)
#    xo=np.arange(10)#+5.
#    yo=np.arange(10)
#    
#    xi = np.array([24, 35, 102, 92])
#    yi = np.array([81, 24, 25, 92])
#    xo = np.array([61, 62, 143, 133])
#    yo = np.array([89, 34, 38, 105])
#     
#    degree=1
    
    
    m = np.size(xi) # no. of points
# ** error checking. check no.s of elts.
    if (np.size(xo) <> m):
        print 'no. of elements must be same in xo,yo,xi,yi'
        stop()
    
# 
    n = degree        #;use halls notation
    n2=(n+1)**2
    #** if n2 gt m then message, '# of points must be ge (degree+1)^2.'
    if (n2>m):
        print # of points must be >= (degree+1)^2.'
        stop()
        

    x = np.transpose(np.array([np.transpose(xi[:]),np.transpose(yi[:])]))
    u = np.transpose(np.array([np.transpose(xo[:]),np.transpose(yo[:])]))

    ut=np.zeros((m,n2)) #;transpose of U
    u2i = np.zeros(n+1)    #;[1,u2i,u2i^2,...]
    
    for i in np.arange(0,m):
        u2i[0]=1 # ;init u2i
        zz = u[i,1]
        for j in np.arange(1,n+1): u2i[j]=u2i[j-1]*zz
        ut[i,0:n+1]=u2i # ;evaluate 0 th power separately
        for j in np.arange(1,n+1): ut[i,j*(n+1):j*(n+1)+n+1]=u2i*u[i,0]**j # ;fill ut=u0i^j * U2i
#    #
    # above is correct, by c.f. IDL

#    ut=polyvander2d(x,u,(n,n))
#    mut=np.matrix(ut)
    
    uu = np.dot(np.transpose(ut),ut) # ;big u
 #   uu=mut*mut.transpose()
    
    kk = np.linalg.inv(uu) #**** need to add error checking
    okk=kk
    
    kk = np.dot(ut,kk)
    #fkk = np.matrix(kk)*np.matrix(ut)
    #kk=0
    #kk=fkk
    #fkk=0
    
    #kx = np.zeros(n+1,n+1) + (np.multiply(kk,np.transpose(x[:,0])))
    #ky = np.zer
    nkk=np.matrix(np.reshape(kk,-1))
    xT=np.matrix(x).transpose()
    #kx = np.zeros(n+1,n+1) + (xT[0,:]*kk)
    #ky = np.zeros(n+1,n+1) + (xT[1,:]*kk)
    
    
    kx=np.reshape( np.matrix(kk.T)*np.matrix(xT[0,:].T) ,(n+1,n+1))
    ky=np.reshape( np.matrix(kk.T)*np.matrix(xT[1,:].T), (n+1,n+1))
    
    return kx.A,ky.A # turn back into arrays
    

def xcor2(inspec,refspec,maxshift,dx=1.,quiet='false'):#       tshift=xcor(tflux,tay0,tay2,maxshift)
    # do crude discrete cross-correlation on each slit (between y0,y1) to find best offset:
    # make running vector from -maxshift,+maxhsift
    ivec=np.arange(-maxshift,maxshift,dx) # doh! - much easier!
        
    refvec=refspec
    
    xcorval=np.zeros(len(ivec))
    for ii in range(len(ivec)):
        shifrefvec=np.roll(refvec,int(ivec[ii]))
        
#        print np.shape(shifrefvec)
#        print np.shape(inspec)
        xcorval[ii]=np.sum(shifrefvec * inspec)
#    if(quiet=='false'):
#        print ii,ivec[ii],xcorval[ii]
        
    pk=np.argmax(xcorval)
#    stop()
    if(quiet=='false'):
        if ((pk==0) | (pk==(len(xcorval)-1))):
            print
            print 'WARNING: xcor2 peak at edge of window!'
            print 'CHECK RESULTS'
            print
    return ivec[pk],xcorval[pk]



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
#    if(quiet=='false'):
#        print ii,ivec[ii],xcorval[ii]

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






def cenlines(speclam,specflux,templam,tempflux,lines,chunksize=100.,linewidth=5.,maxshift=50.,iplot=False,idebug=False): # all numbers in A.
    # Take an approximately-aligned template and arc and a list of clean lines and "fine-tune":

    Nsig=3.
    # -- First clean template using linelist. only keep lines +/-N*linewidth from line positions
    goodmask=np.zeros(np.size(templam))
    for line in lines: # (wavelength of selected lines)
        ok=np.where(np.abs(templam-line) <= Nsig*linewidth)
        # might be off spectrum!
        if (np.size(ok)>0):
            goodmask[ok]=1
    # apply mask:
    bad=np.where(goodmask==0)
##    tempflux[bad]=0. # not masking seems to work better!

    if (iplot): plt.plot(templam,tempflux,'k-')
    if (iplot): plt.title('template spectrum and line list')
    if (iplot): plt.show()

    # -- we have spectra quantised in pixels but wth approx WLC attached to wavelength arrays.
    #    let's work out approx. shifts necessary for better alignment of sub-chunks:

    alines=[]
    alamshifts=[]
    for line in lines:
        # take data at line wavelength +/- chunksize. ////if this line has already been considered as part of another chunk, then skip
        tempok=np.reshape(np.where(np.abs(templam-line) <= Nsig*linewidth),-1)
        if(np.size(tempok)==0): continue
        specok=np.reshape(np.where(np.abs(speclam-line) <= Nsig*linewidth),-1)
        if(np.size(specok)==0): continue
        # xcor chunks:
        if(np.size(specok) != np.size(tempok)): continue # may need to resample corrected wavelength scale for d(pixelsize)


        # assumes pixels are aligned! sort out angstrom shifts later
        refspec=tempflux[tempok]
        inspec=specflux[specok]
        pixshift,corrval = xcor2(inspec,refspec,maxshift,dx=1.,quiet='true')
#        print np.shape(tempok),np.shape(specok)
        
        dw=templam[tempok[1]]-templam[tempok[0]]
        lamshift = pixshift * dw #+templam[tempok[0]]-speclam[spec[ok]] # or summat... ***
        #pl.plot(line+lamshift,line,'D',mec='r',mfc='none') # this should indicate where the line *was* in the input spectrum
        tline=line+lamshift
        if (iplot): plt.axvline(tline,color='gray' )
        # we need to *add* this shift to correct the spectrum
        
        if (idebug): print line,pixshift
        #-- don't include shifts at edge of xcor window!:
        ###if( (np.abs(pixshift)>0) & (np.abs(pixshift)<(maxshift-1)) ):
        if(  (np.abs(pixshift)<(maxshift-1)) ): # WHY WOULD WE REJECT 0 shifts??
            alines=np.append(alines,line)
            alamshifts=np.append(alamshifts,lamshift)
    return alines,alamshifts
        

def robust_poly_fit(xList,yList,order=1,idebug=False):
    
    '''fit the data using a least squares and polynomial'''
    A = np.vander(xList,order+1)
    b = np.matrix(yList).T
    (w,residuals,rank,sing_vals) = linalg.lstsq(A,b)
    if (idebug): 
        print w
    #--
        print A,np.shape(A)
        print np.shape(b)
    # DGG: My own modification to allow biweight rejection of outliers:
    X=A
    #X=np.c_[A,np.ones(np.shape(A))]
    #X=np.c_[xList,len(xList)]
    y2=b
    # resrlm2 = sm.RLM(y2, X,).fit()
    resrlm2 = sm.RLM(y2, X,M=sm.robust.norms.TukeyBiweight()).fit()
    coeff=resrlm2.params
    if (idebug): print coeff
    w=coeff
    #--
    return w.T.tolist()##[0]




###def zeroshift(filename,reflam,refspec,maxshift=50.):
def zeroshift(inspec,reflam,refspec,maxshift=50.,idebug=False,globshift=0.,smooth=0.):

    # represent pix-->wavelength in reflam as cubic
    # this has been checked to be ana accurate representation for this arc.
    xpix=np.arange(len(reflam))
    coeff=np.polyfit(xpix,reflam,3)
    yf=np.polyval(coeff,xpix)
    
###    inspec,hdr=pyfits.getdata(filename,header='true')
###    dx=-float(hdr['xp']) + 1700. # shift to match refarc.txt pixels
###    if (idebug==1): print dx
##    dx=0.
    dx=globshift # global shift on input (to all RSS to work)
#    if(idebug==1):
#        plt.figure()
#        plt.title(filename+'  '+str(dx))
###        plt.plot(reflam+dx,inspec/10.,'b--')
#        rline=plt.plot(reflam,refspec,'k-')    
    #pl.show()
    
    if (smooth>1.):
        #nn=5
        nn=smooth
        kern=np.ones(nn)/float(nn)
        srefspec=np.convolve(refspec,kern,mode='same')
        refspec=srefspec
###
#    refspec=srefspec
###
                        
###    foo=np.roll(np.reshape(inspec,-1),int(dx))/10.
#    nn=10
#    kern=np.ones(nn)/float(nn)
#    foo=np.convolve(inspec,kern,mode='same')
    
##    foo=inspec
    foo=np.roll(np.reshape(inspec,-1),int(dx))#/10.
    
    ###
#    if (idebug==1): 
#        print np.shape(inspec)
###        print np.shape(dx)
        ###plt.plot(inspec/1.,'g:')
    
    pixshift,corrval=xcor2(foo,refspec,maxshift=maxshift)
    #if (idebug==1): print pixshift
    foo2=np.roll(np.reshape(foo,-1),int(-pixshift))
    #if (idebug==1): plt.plot(foo2/1.,'r-')
    # works in pixels, now need to shift in A accordingly
    totdx=dx-pixshift
    #**** this isn't quite right yet:
    #midx=int(round(float(len(reflam)/2.)))
    #totdl=totdx*(reflam[midx]-reflam[midx-1])
    #speclam=reflam-totdl
    speclam=np.polyval(coeff,(xpix+totdx)) # shift cubic fit for wl soln.
    if (idebug):
#        plt.figure()
        nomline=plt.plot(xpix,foo,'b--')
        rline=plt.plot(xpix,refspec,'k-')  
        fline=plt.plot(xpix+totdx,inspec,'r--')
#        plt.legend((rline,nomline,fline),('ref','nominal','fitted'))
        plt.legend(('nominal','ref','fitted'))
    #    plt.plot(speclam,inspec,'r-')
#        plt.show()
    
    
    
    # -- write spectrum and wavelength solution to file:
    #outfile=re.sub('.ms.fits','_wlc0.txt',filename)
    #np.savetxt(outfile,np.transpose((speclam,inspec)))
    return (speclam,inspec)
#    stop()

"""
def gaussfit(x,y,xPeakGuess,sigmaGuess=5.0,nSigWidth=3.0):
    trimmed=np.where(np.abs(x-xPeakGuess)<(nSigWidth*sigmaGuess))[0]
    # this is for histogram of data? not x,y fit
    (mu1,sigma1) = norm.fit(x[trimmed])
    ###(mu1,sigma1) = norm.fit(p1)
    y1 = mlab.normpdf( bins, mu1, sigma1)
    plt.figure()
    a1=plt.hist(p1,bins=np.arange(-1,1,0.04),normed=1)
    l1 = plt.plot(bins, y1, 'r--', linewidth=2)    
"""


def interactiveWLC(inlam,inspec,lines,pkwin=10.0):#,reftemp=None):
    # click on lines in spectrum and lines in linelist to interactively centroid and fit
    """
    plt.figure()
    plt.plot(inlam,inspec,'r-')
    for i in range(len(lines)):
        if ( (lines[i]<np.min(inlam)) |  (lines[i]>np.max(inlam)) ): continue
        plt.axvline(lines[i],color='gray')
    plt.show()    
    """
    
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    
    textsize = 9
    left, width = 0.1, 0.8
    #rect1 = [left, 0.7, width, 0.2]
    #rect2 = [left, 0.3, width, 0.4]
    #rect3 = [left, 0.1, width, 0.2]
    
    lamdiff=[]
    lamline=[]
    
    rect1 = [left, 0.4, width, 0.55]
    rect2 = [left, 0.1, width, 0.3]

    fig = plt.figure(facecolor='white')
    axescolor  = '#f6f6f6'  # the axes background color
 
  
    ax1 = fig.add_axes(rect1, axisbg=axescolor)  #left, bottom, width, height
    ax2 = fig.add_axes(rect2, axisbg=axescolor, sharex=ax1)

    
    ax2.set_ylim([-10,10])
    ax2.set_xlim([np.min(inlam),np.max(inlam)])
#    ax2.set_xlabel('log(M_* Brinchmann)')
    #ax2.set_ylabel('log(M_* magphys,'+model+')-log(M_* Brinchmann)')
    ax2.set_ylabel('difference')
#    ax2.plot(
    plt.show()

    
    #ax2t = ax2.twinx()
    #ax3  = fig.add_axes(rect3, axisbg=axescolor, sharex=ax1)
    ax1.plot(inlam,inspec,'r-')
    #**
#    matplotlib.widgets.Cursor(ax1)
    #**
    for i in range(len(lines)):
        if ( (lines[i]<np.min(inlam)) |  (lines[i]>np.max(inlam)) ): continue
        ax1.axvline(lines[i],color='gray')
    plt.show()    

    



    
    
    
        
    key=''
    print 'left click line in spectrum followed by expected line position'
    print 'middle click to stop'
    print
    y0=-10.0
    y1=10.0
    while(key<>'n'):
#        key=raw_input('Adjust points (y/n)?')
    
#        if (key=='n'): break

        #--
        #plt.show()
        #--
        
        x=plt.ginput(2,show_clicks=True) # wait for two clicks
        # centroid on 1st coord (spec) and tag second (known line lambda)
        if(np.size(x)==0): break
        roughObsLambda=x[0][0]
        texpectLambda=x[1][0]
        #**

        inwin=np.where( np.abs(inlam-roughObsLambda) <= (pkwin/2.0) )[0]
        cenObsLambda=inlam[inwin][np.argmax(inspec[inwin])]
#        cenObsLambda= roughObsLambda #*** for now
        ax1.axvline(cenObsLambda,color='yellow',alpha=0.7)

        expectLambda=lines[np.argmin(np.abs(lines-texpectLambda))]
        
        #**
        lamdiff=np.append(lamdiff,cenObsLambda-expectLambda)
        lamline=np.append(lamline,expectLambda)
        print "%.1f --> %.1f"%(cenObsLambda,expectLambda)
        ax1.axvline(expectLambda,color='green')
        ax2.plot(lamline,lamdiff,'ok')
        y0=np.min([y0,lamdiff.min()])
        y1=np.max([y1,lamdiff.max()])
        ax2.set_ylim([y0,y1])
        plt.show()
        
        
        
    # fit results:
    coeff = robust_poly_fit(lamline,lamdiff,order=3)
    tweaks=np.polyval(coeff,inlam)
    outlam=inlam-tweaks
    #xxx=np.arange(inlam.min(),inlam.max(),1)
    ax2.plot(inlam,np.polyval(coeff,inlam),'g--')
    ax1.plot(outlam,inspec,'g-')
    plt.show()
    return outlam
    
        

# linelist='arc78_wlc.lis'
###def wlc_arc(arcroot,linelist='_tmparc.lis',chunksize=200.,linewidth=5.,maxshift=20.,tweakorder=2,no_tweak='false'):
def wlc_arc(arcspec,reflam,refspec,linelist='_tmparc.lis',chunksize=200.,linewidth=5.,maxCoarseShift=20.,maxFineShift=10.,\
                tweakorder=2,no_tweak='false',globshift=0.,idebug=False, lenfrac=0.50,initorder=2,iplot=True):
#    files=glob.glob(arcroot+'_??.ms.fits')
#    reflam,refspec=np.loadtxt('/home/gilbank/Proj/imacs/RSS/20120318/product/TEST/refarc.txt',unpack='t')
    
    # lenfrac -- use this percentage of wavelength range for initial fit of wavelength solution

    try:
        xp,lam0,lam1=np.loadtxt(linelist,usecols=(0,1,2),unpack='t')
        lines=lam0
    except:
        lines = np.loadtxt(linelist) #*** just list of wavelengths of lines
    #**
#    lines= np.append([5502.,5529.],lam0)# manually add extra lines where sparse
#    evenlines=np.arange(4500,7500,chunksize/2.)
#    lines=np.append(lines,evenlines)
    #**
    if (iplot): fig=plt.figure(figsize=(15,10))
#    for filename in files:
    # -- first, check for small global shift (relative to rsmt posn.):
    #print 'finding zeropoint shift from reference arc...'
    #speclam,specflux=zeroshift(filename,reflam,refspec,maxshift=50)
    speclam,specflux=zeroshift(arcspec,reflam,refspec,maxshift=maxCoarseShift,globshift=globshift)
    # speclam,specflux now contain solution after only a [small] global shift
    # -- now fit more refined solution:
    if (idebug & iplot): fig=plt.figure()
    #plt.title='wavelength calibration (fine tuning)' 
    if (iplot): a = fig.add_subplot(2,2,1)

#    alines,alamshifts=cenlines(speclam,specflux,reflam,refspec,lines,chunksize=chunksize,linewidth=linewidth,maxFineShift=10.)
###    alines,alamshifts=cenlines(speclam,specflux,reflam,refspec,lines,chunksize=chunksize,linewidth=linewidth,maxshift=10.)
    print iplot
    alines,alamshifts=cenlines(speclam,specflux,reflam,refspec,lines,chunksize=chunksize,\
                               linewidth=linewidth,maxshift=maxFineShift,idebug=idebug,iplot=False)
    
    #print flibbe
    #**
    #plt.quiver(alines+alamshifts,np.copy(alines)*0.+15000.,-alamshifts/10.,np.copy(alines)*0.,units='dots')
    #**
    if (iplot): a = fig.add_subplot(2,2,3)

    #-- store
    alines0=np.copy(alines)
    alamshifts0=np.copy(alamshifts)
    #--

    if (iplot):
        pp1=plt.plot(alines,alamshifts,'ok')
        plt.xlabel='r$\lambda(\AA)$'
        plt.ylabel='dr$\lambda$'
    w=alines
    z=alamshifts
    polycoeffs = np.polyfit(w,z, 2)
    
    # **** Might want a little more error-checking in here
    #      but currently, tweaking seems not too bad! ****

    #-- fit centre regions first and iterate:
    usecen=np.median(speclam)
    maxl=np.max(speclam)
    minl=np.min(speclam)
    lamrng=maxl-minl
    ###userng=np.where(np.abs(w-usecen)<=1000.)[0]
    # use only central 50% of wavelnegth range in fit (+/-25%)
    ##userng=np.where(np.abs(w-usecen)<=lamrng/4.)[0]

    if (lenfrac<1.0):
        halflamrng=lamrng*(1.-lenfrac)/2.
        userng=np.where(np.abs(w-usecen)<=halflamrng)[0]
        print 'USING: ',np.min(w[userng]),np.max(w[userng])
    else: 
        userng=range(len(w))  

#   pcoeffs=robust_poly_fit(w,z,order=2)
#    pcoeffs=robust_poly_fit(w,z,order=tweakorder)
##    pcoeffs=robust_poly_fit(w[userng],z[userng],order=tweakorder)
#    try:
##    pcoeffs=robust_poly_fit(w[userng],z[userng],order=2)#tweakorder)
    pcoeffs=robust_poly_fit(w[userng],z[userng],order=initorder)#tweakorder)
#    except:
#        pass
    print pcoeffs
    # iterate and use full range:
    lamcal0=speclam-np.polyval(pcoeffs,speclam)
#    alines2,alamshifts2=cenlines(lamcal0,specflux,reflam,refspec,lines,chunksize=chunksize,linewidth=linewidth,maxFineShift=maxFineShift,iplot=False)
    alines2,alamshifts2=cenlines(lamcal0,specflux,reflam,refspec,lines,chunksize=chunksize,\
                                 linewidth=linewidth,maxshift=maxFineShift,iplot=False)
    #alines2,alamshifts2=cenlines(speclam,specflux,reflam,refspec,lines,chunksize=200.,linewidth=5.,maxshift=10.,iplot=False)
    w2=alines2
    z2=alamshifts2
#    print 'w2,z2:',w2,z2
#    polycoeffs = np.polyfit(w,z, 2)
    pcoeffs0=np.copy(pcoeffs)
#--***
    try:
        pcoeffs=robust_poly_fit(w2,z2,order=tweakorder)
    except:
        print '2nd robust_poly_fit with order=%s failed! Using initial n=%s fit instead'%(initorder,tweakorder)
        
        pcoeffs=np.copy(pcoeffs0) 
#--***        
    # 2nd iteration shifts must all be small! if not, reject:
    tweaks=np.polyval(pcoeffs,speclam)
    print "rms fit %.2f A"%(np.max(np.abs(tweaks)))
##    if ( np.max(np.abs(tweaks))>10.0 ) : 
    if ( np.max(np.abs(tweaks))>25.0 ) : 
        print 'TWEAKS TOO LARGE. Rejecting iterated solution. Beware of results!'
        pcoeffs=pcoeffs0
    else:
        pass
    
    if (iplot):
        ll1=plt.plot(speclam,np.polyval(pcoeffs0,speclam), 'g:')
        pp2=plt.plot(alines2,alamshifts2,'ob')
    print alines2,alamshifts2
    #--
    
    ###pl.plot(speclam,np.polyval(polycoeffs,speclam), 'g-')
    if (iplot):
        ll2=plt.plot(speclam,np.polyval(pcoeffs,speclam), 'g--')
        plt.legend(('initial nrst line offsets (after best global shift)','initial quadratic fit',\
                    'nrst line offsets after n=2 polyfit','final (used) fit'))
    #plt.legend((pp1,pp2),('initial nrst line offsets (after best global shift)'\
    #                              'nrst line offsets after n=%s polyfit'))
        plt.show()

        a = fig.add_subplot(2,2,2)  
##    lamcal=speclam-np.polyval(pcoeffs,speclam)
    lamcal=lamcal0-np.polyval(pcoeffs,speclam)
    # lamcal,specflux now contain the WLC after final cubic correction
    if (iplot):
        l1=plt.plot(reflam,refspec,'k-')
        l2=plt.plot(speclam,specflux,'b:')
        l3=plt.plot(lamcal,specflux,'r-')
    #plt.legend((l1,l2,l3),('template','arc after coarse shift','final tweaked arc'))
        plt.legend(('template','arc after coarse shift','final tweaked arc'))
    
    for ii in range(len(alamshifts0)):
        lok=np.argmin(np.abs(reflam-alines0[ii]))
        if (iplot):
            plt.plot((alines0[ii],alines0[ii]-alamshifts0[ii]),(refspec[lok],refspec[lok]),'g-')              
        else:
            pass
    
    if (iplot):
        a = fig.add_subplot(2,2,4)  
        # this is a calibrated arc!
    
        plt.plot(lamcal,specflux,'r-')
        plt.title('final calibrated arc')
        lmin=np.min(lamcal)
        lmax=np.max(lamcal)
        for iii in range(len(lines)):
            tline=lines[iii]
            #if ( (tline<np.min(lamcal)) | (tline>np.max(lamcal)) ) : continue
            plt.axvline(tline,color='gray')
            plt.xlim([lmin,lmax])
            
            fig.add_subplot(2,2,1)
            plt.xlim([lmin,lmax])
            fig.add_subplot(2,2,2)
            plt.title('best global shift')
            plt.xlim([lmin,lmax])
            fig.add_subplot(2,2,3)
            plt.title('tweak to individual lines and polyfit')
            plt.xlim([lmin,lmax])
            

    #stop()
    ## Write wavelength-calibrated results to .ms.fits file:
###    outfile=re.sub('.ms.fits','_wlc.txt',filename)
    #np.savetxt('_tmp_spec.txt',np.transpose( (lamcal,specflux) ))
###    np.savetxt(outfile,lamcal)
    #dcfile=re.sub('.ms.fits','_dc.ms.fits',filename)
    #iraf.rspectext ('_tmp_spec.txt',\
    #                dcfile, title="", flux="no", dtype="nonlinear", crval1=1., cdelt1=1., fd1="", fd2="")
    #os.system('rm -f _tmp_spec.txt')
    #*** copy useful header keywords from 2D to 1D files??
        
    return lamcal

#https://gist.github.com/250860

#import sys
#from numpy import NaN, Inf, arange, isscalar, asarray
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Currently returns two lists of tuples, but maybe arrays would be better
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
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



def centroid(xvec,yvec):
    data=yvec
    
    # http://www.scipy.org/Cookbook/FittingData#head-a44b49d57cf0165300f765e8f1b011876776502f
    #X = np.arange(data.size)
    X=xvec
    x = np.sum(X*data)/np.sum(data)
    #width = np.sqrt(np.abs(np.sum((X-x)**2*data)/np.sum(data)))
    #max = data.max()
    #fit = lambda t : max*exp(-(t-x)**2/(2*width**2))
    #plot(fit(X))
    return x

def tracelines(image,linexcoords,yref,ystep=2,lineFWHM=5.0,ybord=1,idebug=False):
    # trace lines in the spatial direction, given a list of input coordinates 
    # and a starting y-value (normally centre of image working outwards. 
    # Proceed until trace is lost

    if (idebug):
        plt.figure()
        plt.imshow(image)
        plt.plot(linexcoords,np.repeat(yref,len(linexcoords)),'or')

    xi=[]
    yi=[]
    xo=[]
    yo=[]
    
    ny,nx=np.shape(image)
#    cens=np.zeros((np.round((ny/ystep)),len(linexcoords)))
    # centroid each line:
    for ii,linex in enumerate(linexcoords):
        minx=int(linex-2.0*lineFWHM)
        maxx=int(linex+2.0*lineFWHM)

        # trace centre to top of image:
        for jj,yval in enumerate(np.arange(yref,ny-ybord,ystep)):
            xdata=np.arange(minx,maxx)
            ydata=image[int(yval),minx:maxx]
            xcen=centroid(xdata,ydata)

            minx=int(xcen-2.0*lineFWHM)
            maxx=int(xcen+2.0*lineFWHM)
            
            xi=np.append(xi,linex)
            yi=np.append(yi,yval)
            xo=np.append(xo,xcen)
            yo=np.append(yo,yval)
#            cens[jj,ii]=ycen
        #trace centre to bottom of image:
        minx=int(linex-2.0*lineFWHM)
        maxx=int(linex+2.0*lineFWHM)
      
        for jj,yval in enumerate(np.arange(yref,ybord,-ystep)):
            xdata=np.arange(minx,maxx)
            ydata=image[int(yval),minx:maxx]
            xcen=centroid(xdata,ydata)

            minx=int(xcen-2.0*lineFWHM)
            maxx=int(xcen+2.0*lineFWHM)


#            xcen=centroid(image[yval,minx:maxx])
#            cens[jj,ii]=ycen
            xi=np.append(xi,linex)
            yi=np.append(yi,yval)
            xo=np.append(xo,xcen)
            yo=np.append(yo,yval)
            
        if (idebug): 
            plt.plot(xo,yo,'ok')
            plt.show()
            
    # map transformation:
    kx,ky=polywarp(xi,yi,xo,yo,3)
    
    if (idebug):
        # check:
        xf,yf=applywarp(xo,yo,kx,ky)
        plt.plot(xf,yf,'Dr',mec='r',mfc='none')
        plt.show()
    
    return kx,ky,xo,yo

#def transformByRow(image,kx,ky):
#    ny,nx=np.shape(image)
#    # transform line-by-line without resampling:

"""
#   https://github.com/tiagopereira/python_tips/wiki/Scipy%3A-curve-fitting
def gaussian(B,x):
    ''' Returns the gaussian function for B=m,stdev,max,offset '''
#    return B[3]+B[2]/(B[1]*np.sqrt(2*np.pi))*np.exp(-((x-B[0])**2/(2*B[1]**2)))
    return B[3]+B[2]*np.exp(-((x-B[0])**2/(B[1]**2)))
def errfunc(p,x,y):
   return y-gaussian(p,x)
"""


def gaussian(x, a, b, c, d):
    val = a * np.exp(-(x - b)**2 / c**2) + d
    return val

    
def extract1D(spec2d,width=3.0,idebug=False,rdnoise=6.0,sideBandSub=False):    
    xsec=np.median(spec2d[:,:],1)
    if(idebug):
        plt.figure()
        plt.plot(xsec)
        plt.show()
    pkmax,pkmin=peakdet(xsec,100)
    npks=np.shape(pkmax)[0]
    print npks
    if(npks>1): 
        print "only works for one peak currently!"
        stop()
    
    # fit Gaussian to position
    p0=[pkmax[0][1],pkmax[0][0],width,np.min(xsec)]
    x=np.arange(len(xsec))
    y=xsec

    
    e=np.ones(len(xsec))
    popt, pcov = curve_fit(gaussian, x, y, sigma=e)
    #print flibble
    ###fit = optimize.leastsq(errfunc,p0,args=(x,y))
    #pky,sigma
    print popt
    xx=np.arange(spec2d.shape[0])
    if(idebug):
#        plt.plot(xx,gaussian(np.array(fit[0]),xx))
        plt.plot(xx,gaussian(xx,popt[0],popt[1],popt[2],popt[3]))
#    print flibbe
    # check for tilt in trace:
    # [looks like near enough horizontal]
    fitpos=popt[1]
    fitsigma=np.min([3.0*np.abs(popt[2]),2.0]) # ** careful, sigma can be -ve.!

    if(idebug):
        plt.figure()
        plt.imshow(spec2d)
        plt.axhline(fitpos+fitsigma,color='r')
        plt.axhline(fitpos-fitsigma,color='g')
        plt.show()

    print fitpos,fitsigma
        
    # crude optimal (Horne) extraction (really just Gaussian weighting):
    
    # (object flux = raw flux - sky)
    # sum in spatial dirn  object flux*profile/noise**2 / sum ( profile**2/noise**2 ) : 
    
    #noise=np.sqrt(spec2d + rdnoise**2)
    noise2 = spec2d + rdnoise**2
    
    #** need to do proper sky sub
#    spec2d=spec2d-popt[3] # subtract const. bgd level
#    popt[3]=0.0
    #**
    #-- simple side band sky sub. :
    
    #mask object +/-3 sigma and make side bands 10 pixels wide starting there:
    b1=int(fitpos-3.0*fitsigma)
    b0=b1-10
    b2=int(fitpos+3.0*fitsigma)
    b3=b2+10
    
    sky01=spec2d[b0:b1,:]
    sky23=spec2d[b2:b3,:]
    sky03=np.append(sky01,sky23,axis=0)
    medsky1d=np.median(sky03,0)
    medsky2d=np.transpose(np.repeat(np.reshape(medsky1d,(-1,1)),(spec2d.shape[0]),axis=1))
    
    profile1d=gaussian(xx,popt[0],popt[1],popt[2],popt[3])
    profile2d=np.repeat(np.reshape(profile1d,(-1,1)),(spec2d.shape[1]),axis=1)
    
    #**** STILL NEED TO CHECK WHETHER I'VE NORMALISED THESE TO 1 PIXEL COUNTS OR SUMMED OVER THE APERTURE OR WHAT
    #     should probably produce counts/sec for final 1d spectra or similar
    
    if (idebug):
        plt.figure()
        plt.imshow(profile2d)
        plt.show()
    
    s0=int(np.round(fitpos-fitsigma))
    s1=int(np.round(fitpos+fitsigma))
    numerbrack = (spec2d-medsky2d) * profile2d/noise2
    numer= np.sum(numerbrack[s0:s1,:],0)
    denombrack = profile2d**2 / noise2
    denom = np.sum(denombrack[s0:s1,:],0)
    
    if (idebug):
        plt.figure()
        plt.imshow(numerbrack/denombrack)
        plt.show()    
    
    flux1d=numer/denom
    return flux1d



    
    
    
        

