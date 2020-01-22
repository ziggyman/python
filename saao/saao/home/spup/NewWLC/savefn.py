import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import warnings
from astropy.io import ascii,fits

from glob import glob
import json

# This assumes 2D Legendre 

def savefn(p,savfile):
    names = [x.encode('ascii') for x in p.param_names]
    pars = p.parameters
    xord = p.x_degree
    yord = p.y_degree
    #return xord,yord,names,pars
    np.savez(savfile,names=names,pars=pars,xord=xord,yord=yord)

#def restorefn(xord,yord,names,pars):
def restorefn(savfile):
    f=np.load(savfile)
    names=f.f.names
    pars=f.f.pars
    xord=f.f.xord#[0] # arrays to ints
    yord=f.f.yord#[0]
    print yord.shape
    d = {}
    for k,v in zip(names,pars):
        d[k]=v
    f = models.Legendre2D(xord,yord,**d)
    return f

    
