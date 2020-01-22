import numpy as np

def zscale(image, contrast=1.0):
   """Implementation of the IRAF zscale algorithm to find vmin and vmax parameters for the dynamic range of a display. It finds the image values near the median image value without the time consuming process of computing a full image histogram."""

   from scipy import optimize
   #import matplotlib.pyplot as plt

   # Get ordered list of points
   I=np.sort(image.flatten())

   # Get number of points
   npoints=len(I)

   # Find the midpoint (median)
   midpoint=(npoints-1)/2

   # Fit a linear function
   # I(i) = intercept + slope * (i - midpoint)

   fitfunc = lambda p, x: p[0]*x+p[1]
   errfunc = lambda p, x, y: fitfunc(p, x) - y

   # Initial guess for the parameters
   p0 = [(I[-1]-I[0])/npoints,I[midpoint]] 

   # Fit
   i=np.arange(len(I))
   p1, success = optimize.leastsq(errfunc, p0[:], args=(i, I))

#    plt.plot(i,I,'r+')
#    plt.plot(i,fitfunc(p1,i))
#    plt.show()

   if success in [1,2,3,4]:
       slope=p1[0]
       z1=I[midpoint]+(slope/contrast)*(1-midpoint)
       z2=I[midpoint]+(slope/contrast)*(npoints-midpoint)
   else:
       z1=np.min(image)
       z2=np.max(image)

   return z1, z2

