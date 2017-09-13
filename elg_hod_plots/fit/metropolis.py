#! /usr/bin/env python

"""
Metropolis algorithm and symmetric jump functions

NOTE: this module requires the numpy and random libraries
"""

import sys
import numpy as np
import random as random

# Likelihood of a normal distribution
def dnorm(val,mu,verr=None,vmin=None,vmax=None,logd=True):
    nv = len(val)

    if(verr is None):
        verr = np.zeros(shape=(nv)) ; verr.fill(1)
    s2 = 2*np.sqrt(verr)

    if (vmin is None and vmax is None):
        if(logd):
            dnorm = np.log(s2*np.pi)/2. -(val-mu)**2/s2
        else:
            dnorm = exp(-(val-mu)**2/s2)/np.sqrt(s2*np.pi)
    else:
        print 'dnorm: not ready to handle Gaussians with limits'
    return dnorm

# Likelihood of a uniform function
def duni(val,vmin,vmax,logd=True):
    nv = len(val)

    duni = np.zeros(shape=(nv))    
    duni.fill(-np.inf)

    for i in range(nv):
        if(val[i]>=vmin[i] and val[i]<=vmax[i]):
            duni[i] = 1./(vmax[i]-vmin[i])            
            if(logd):
                duni[i] = np.log(duni[i])
    return duni

# Likelihood of a chi^2 function
def chi2_likelihood(params, func, x, y, yerr=None, logd=True):
    if(yerr is None):
        yerr = np.zeros(shape=(len(y))) ; yerr.fill(1)

    model = func(x,params)
    chi2 = sum((y-model)**2/(yerr**2))

    if (logd):
        fchi2 = -chi2/2.
    else:
        fchi2 = np.exp(-chi2/2.)

    return fchi2

# Get the unnormalised posterior probability, i.e.
# the likelihood(d|h)*prior(h), h= hypothesis, d=data
def posterior(x,y,params,lb,ub,funcp,funcl,model,logd=True):
    likelihood = funcl(params,model,x,y,logd=logd)

    priors = funcp(params,lb,ub,logd=logd)
    if (logd):
        prior = sum(priors)
        posterior = likelihood+prior
    else:
        prior = np.prod(priors)
        posterior = likelihood*prior

    return posterior


def metropolis(x,y,params,lb,ub,funcp,funcl,model,stepsize=None,nit=1e5):
    # The Metropolis algorithm needs the proporsal or jumping density
    # to be symmetrical (such a Gaussian, without limits)
    random.seed()    
    npass = 3

    A = [params] #; print A
    L = [posterior(x,y,params,lb,ub,funcp,funcl,model)]  

    if(stepsize is None):
        stepsize = (ub-lb)/1000.

    should_restart = True ; ipass = 0 
    while should_restart:
        accepted = 0.0
        should_restart = False
        if(ipass==0):
            # Find a reasonable starting point
            iit = 1000
        else:
            # Refine the starting point
            iit = int(nit/npass)
 
        for i in range(iit):
            oldA = A[len(A)-1]
            oldL = posterior(x,y,oldA,lb,ub,funcp,funcl,model)
        
            # Assuming a gaussian distribution
            newA = np.zeros(len(params))
            for ip in range(len(params)):
                newA[ip] = random.gauss(oldA[ip],stepsize[ip])
            newL = posterior(x,y,newA,lb,ub,funcp,funcl,model)
        
            #Check if to accept or not the new value:
            if(np.isinf(newL)):
                A.append(oldA) ; L.append(oldL)
            else:
                if(np.isnan(newL)):
                    sys.exit('STOP: newL is NAN')
                vv = min(newL-oldL,100.) #To avoid too large numbers
                v = min(np.exp(vv),1.)
                u = random.uniform(0.0,1.0) #u within [low, high)
                if(v>u):
                    A.append(newA) ; L.append(newL)
                    accepted = accepted +1.0            
                else:
                    A.append(oldA) ; L.append(oldL)

        acc = accepted/float(iit)
        if (ipass<npass):
            should_restart = True
            ipass = ipass + 1
            ind = np.argmax(L)
            A = [A[ind]] ; L = [L[ind]]
            print 'Acceptance = ',acc,' ; Restarting from:'
            print '      ',A
            print '      ',L

    A = np.array(A) ; L = np.array(L)
    
    print 'Metropolis, Acceptance=',accepted/float(nit)
    print '            Best value=',A[np.argmax(L),:]

    return A,L 

