#! /usr/bin/env python
import sys
import numpy as np
from scipy.special import erf,erfc
from fractions import Fraction

def lCen_Z05(x,Mmin,sig): # Eq. from Zheng, Z et al. 2005
    y = np.zeros(shape=(len(x)))
    if (Mmin>0. and sig>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_Z05, WARNING: Negative masses as input'
        y[ind] = 0.5*(1.+ erf((np.log10(x[ind])-np.log10(Mmin))/sig))
    return y
def lSat_Z05(x,Msat,asat,Mcut):
    y = np.zeros(shape=(len(x)))
    if (Msat>0.):
        ind = np.where(x>Mcut) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lSat_Z05, WARNING: Negative masses as input'
        y[ind] = ((x[ind]-Mcut)/Msat)**asat
    return y 

def lCen_Z05y(x,fc,Mc,sigc): # Eq. from Zheng, Z et al. 2005 for young galaxies
    y = np.zeros(shape=(len(x)))
    if (Mc>0. and sigc>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_Z05y, WARNING: Negative masses as input'
        y[ind] = fc/(1.+np.exp((np.log10(x[ind])-np.log10(Mc))/sigc))
    return y
def lSat_Z05y(x,fs,M0,sigs):
    y = np.zeros(shape=(len(x)))
    if (M0>0. and sigs>0.):
        ind = np.where(x>M0) 
        y[ind] = fs/np.exp((np.log10(x[ind])-np.log10(M0))/sigs) # The slope seems wrong
    return y

def lCen_Z07(x,Mmin,sig): # Eq. from Zheng, Z et al. 2007
    return lCen_Z05(x,Mmin,sig)
def lSat_Z07(x,ncen,Msat,asat,Mcut):
    y = np.zeros(shape=(len(x)))
    if (Msat>0. and Mcut>0.):
        ind = np.where(x>Mcut)
        y[ind] = ncen[ind]*((x[ind]-Mcut)/Msat)**asat
    return y 


def lCen_T12(x,params): # Eq. from Tinker, J et al. 2012
    return lCen_Z05(x,params)
def lSat_T12(x,params):
    Mmin = params[0]
    sig = params[1]
    Msat = params[2]
    asat = params[3]
    Mcut = params[4]

    ncen = lCen_T12(x,[Mmin,sig])

    y = np.zeros(shape=(len(x)))
    ind = np.where(x>0.)
    y[ind] = ncen[ind]*(x[ind]/Msat)**asat*np.exp(-Mcut/x[ind])
    return y


def lCen_G12(x,fa,fb,Mc,sigc): # Eq. from Geach, J et al. 2012
    y = np.zeros(shape=(len(x)))
    if (Mc>0. and sigc>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_G12, WARNING: Negative masses as input'
        gau = np.exp(-(x[ind]-Mc)**2./(2.*sigc**2.))
        step = 1.+erf((x[ind]-Mc)/sigc)
        y[ind] = fb*(1.-fa)*gau + fa*step

    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y
def lSat_G12(x,fs,Mmin,sigs,asat):
    y = np.zeros(shape=(len(x)))
    if (Mmin>0. and sigs>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_G12, WARNING: Negative masses as input'
        y[ind] = fs*(1.+erf((x[ind]-Mmin)/sigs))*\
            (np.np.exp(x[ind])/np.np.exp(Mmin))**asat

    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y


def lCen_C13(x,fa,fb,siga,sigb,Mc): # Eq. from Contreras, S et al. 2013
    gau = np.zeros(shape=(len(x)))
    step = np.zeros(shape=(len(x)))

    ind = np.where(x<Mc)
    gau[ind] = fb*(1.-fa)/np.exp((x[ind]-Mc)**2./(2.*siga**2.))
    step[ind] = fa*(1.+erf((x[ind]-Mc)/siga))

    ind = np.where(x>=Mc)
    gau[ind] = fb*(1.-fa)/np.exp((x[ind]-Mc)**2./(2.*sigb**2.))
    step[ind] = fa*(1.+erf((x[ind]-Mc)/sigb))

    y = gau + step

    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y
def lSat_C13(x,params):
    return lSat_G12(x,params)


def lCen_gasym(x,p0): 
    # Using an np.exponentially modified Gaussian
    # The masses are np.expected as np.log10(M)
    fb,mb,sigb,fd,md,sigd = p0

    y = np.zeros(shape=(len(x))) ; y.fill(-np.inf)
    if (mb>0. and sigb>0. and md>0. and sigd>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_gasym, WARNING: Very small values as input'

        # Bulges ~ Step function
        step = (1.+ erf((x[ind]-mb)/sigb))

        # Disks ~ Asymmetric Gaussian
        inexp = (2*md + fd*sigd**2 - 2*x[ind])*fd/2. 
        mgau = np.exp(inexp)*erfc((md+fd*sigd**2-x[ind])/(np.sqrt(2.)*sigd))

        # Normalised total
        y[ind] = step*fb/2. + mgau*fd/2.

    # Take log10
    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])

    return y

def lCen_g(lx,p0): 
    # Using a step function plus a Gaussian
    # The masses are np.expected as np.log10(M)
    fb,mb,sigb,fd,md,sigd = p0

    y = np.zeros(shape=(len(lx))) ; y.fill(-999.)
    if (mb>0. and sigb>0. and md>0 and sigd>0.):
        ind = np.where(lx>0.) ; negs = np.shape(ind)[1]-np.shape(lx)[0]
        if (negs>0):
            print 'lCen_g, WARNING: Negative values as input'

        # Bulges ~ Step function
        step = (1.+ erf((lx[ind]-mb)/sigb))

        # Disks ~ Gaussian
        inexp = ((lx-md)**2)/(-2.*sigd**2)
        gau = np.exp(inexp)/np.sqrt(2.*np.pi*sigd**2)

        # Normalised total
        y[ind] = step*fb/2. + gau*fd/2.
        
    ind = np.where(y<=0.)
    if (np.shape(ind)[1]>1):
        y[ind] = -np.inf
    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y

def lCen_gasym2(x,p0): 
    # Using an np.exponentially modified Gaussian
    # The masses are np.expected as np.log10(M)
    fb,mb,fd,ad,md,sigd = p0

    y = np.zeros(shape=(len(x))) ; y.fill(-np.inf)
    if (mb>0. and md>0. and sigd>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_gasym2, WARNING: Very small values as input'

        # Bulges ~ Step function
        step = (1.+ erf((x[ind]-mb)/sigd))

        # Disks ~ Asymmetric Gaussian
        inexp = (2*md + ad*sigd**2 - 2*x[ind])*ad/2. 
        mgau = np.exp(inexp)*erfc((md+ad*sigd**2-x[ind])/(np.sqrt(2.)*sigd))

        # Normalised total
        y[ind] = step*fb/2. + mgau*fd/2.

    # Take log10
    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])

    return y

def lCen_g2(lx,p0): 
    # Using a step function plus a Gaussian
    # The masses are np.expected as np.log10(M)
    fb,mb,fd,md,sigd = p0

    y = np.zeros(shape=(len(lx))) ; y.fill(-999.)
    if (mb>0. and md>0 and sigd>0.):
        ind = np.where(lx>0.) ; negs = np.shape(ind)[1]-np.shape(lx)[0]
        if (negs>0):
            print 'lCen_g, WARNING: Negative values as input'

        # Bulges ~ Step function
        step = (1.+ erf((lx[ind]-mb)/sigd))

        # Disks ~ Gaussian
        inexp = ((lx-md)**2)/(-2.*sigd**2)
        gau = np.exp(inexp)/np.sqrt(2.*np.pi*sigd**2)

        # Normalised total
        y[ind] = step*fb/2. + gau*fd/2.
        
    ind = np.where(y<=0.)
    if (np.shape(ind)[1]>1):
        y[ind] = -np.inf
    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y

def lSat_gasym(x,Mcut,Msat,asat):
    y = np.zeros(shape=(len(x)))  ; y.fill(-999.)
    if (Mcut>0. and Msat>0.):
        ind = np.where(x>Mcut) 
        y[ind] = ((x[ind]-Mcut)/Msat)**asat

    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y 

def lg(x,p0): 
    fg,mg,sigg = p0

    y = np.zeros(shape=(len(x))) ; y.fill(-999.)
    if (mg>0 and sigg>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_gasym, WARNING: Negative values as input'

        # Gaussian
        inexp = -0.5*((x[ind]-mg)/sigg)**2
        y[ind] = fg*np.exp(inexp)
        
    ind = np.where(y<=0.)
    if (np.shape(ind)[1]>1):
        y[ind] = -np.inf
    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y


def lgasym(x,p0): 
    # Using an np.exponentially modified Gaussian
    # The masses are np.expected as np.log10(M)
    fd,md,sigd = p0

    y = np.zeros(shape=(len(x))) ; y.fill(-999.)
    if (md>0 and sigd>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_gasym, WARNING: Negative values as input'

        # Disks ~ Asymmetric Gaussian
        inexp = (2*md + fd*sigd**2 - 2*x[ind])*fd/2.
        mgau = np.exp(inexp)*erfc((md+fd*sigd**2-x[ind])/(np.sqrt(2.)*sigd))

        # Normalised total
        y[ind] = mgau*fd/2.

    ind = np.where(y<=0.)
    if (np.shape(ind)[1]>1):
        y[ind] = -np.inf
    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y


def lstep(x,p0): 
    fb,mb,sigb = p0

    y = np.zeros(shape=(len(x))) ; y.fill(-999.)
    if (mb>0. and sigb>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print 'lCen_gasym, WARNING: Negative values as input'

        # Bulges ~ Step function
        step = (1.+ erf((x[ind]-mb)/sigb))

        # Normalised
        y[ind] = step*fb/2. 
        
    ind = np.where(y<=0.)
    if (np.shape(ind)[1]>1):
        y[ind] = -np.inf
    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return y
