import os.path, sys
import numpy as np
import hod_functions as hf
#import metropolis as m
import emcee
from matplotlib import pyplot as plt
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

# log(likelihood)
def lnprob(params, x, y, erry, lb, ub, funcm):
    priors = m.duni(params,lb,ub)

    lprior = sum(priors)
    if not np.isfinite(lprior):
        # if the parameters are outside the prior range
        return -np.inf

    # Likelihood of a chi^2 function
    likelihood = m.chi2_likelihood(params, funcm, x, y, yerr=erry)

    posterior = lprior+likelihood
    return posterior

#############################
funcm = hf.lCen_gasym # hf.lCen_g

# Starting point : fb,mb,sigb,fd,md,sigd
p0 = np.array([0.07, 13.5, 1.9, 2.9, 11., 0.09])

# Limits for priors
lb = np.array([0.007, 8.5, 0.01,0.007, 8.5,0.01])
ub = np.array([ 3., 14.,  50.,  3.,14., 50.])

surveys = ['DEEP2']#, 'VVDS-DEEP', 'VVDS-WIDE', 'eBOSS', 'DESI']
snap_list = [44]#[44, 42, 41, 40, 37] #MillGas
zleg = 'z=0.62'

cols = get_distinct(len(surveys))

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
model = 'MillGas/gp15newmg/'

outdir = '/gpfs/data/violeta/lines/desi_hod_o2/'
plotfile = outdir+'plots/hodfits/'+model+'hod_OII3727.pdf'

##############################

# Loop over redshifts and surveys
for survey in surveys:
    for iz,zsnap in enumerate(snap_list):
        # Input generated with:
        infile = outdir+'hod/'+model+'z'+str(snap_list[iz])\
            +'_OII3727_'+survey+'_bot_central.txt'

        if(os.path.isfile(infile)):
            ffile = open(infile,'r')
            data = ffile.readlines()
            ffile.close()

            # Count number of lines that are not the header
            nl = 0
            for line in data:        
                if line[0].isdigit():
                    nl = nl + 1
            print nl,' read lines in ',infile
    
            # Initialize vectors
            logM, hod, cen, sat, cs, cd  = \
                [np.zeros(shape=(nl)) for i in range(6)]
            # Store data
            nl = 0
            for line in data:
                if line[0].isdigit():
                    logM[nl] = float(line.split(',')[0])
                    hod[nl] = float(line.split(',')[1])
                    cen[nl] = float(line.split(',')[2])
                    cs[nl] = float(line.split(',')[4])
                    cd[nl] = float(line.split(',')[5])
                    
                    nl = nl + 1            

            # Data to be fit
            ind = np.where((cen > 0.001))
            x = logM[ind] 
            y = np.log10(cen[ind])
            erry = np.sqrt(cen[ind])

            # Set up the MCMC sampler
            ndim = len(p0) # number of parameters being fit
            # number of walkers to use (should be > 2*ndim)
            nwalkers = 500 #100
            # start each walker in a small ball around this position
            pos = [p0 + 1e-2 * np.random.randn(ndim) \
                       for i in range(nwalkers)]

            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, \
                                                args=(x, y, erry, lb, ub, funcm))
            sampler.run_mcmc(pos, 1000)
            # The output, sampler.chain has shape [nwalkers, nsteps, ndim]

            # Plot 
            plotmcmc = False
            if (plotmcmc):
                # fb,mb,sigb,fd,md,sigd
                plt.figure(1, figsize=(12, 8))
                plt.subplot(3, 2, 1)
                plt.plot(sampler.chain[:, :, 0].T, alpha=0.05, color='k')
                plt.ylabel(r'$fb$') ; plt.xlabel('step')
                plt.subplot(3, 2, 2)
                plt.plot(sampler.chain[:, :, 3].T, alpha=0.05, color='k')
                plt.ylabel(r'$fd$') ; plt.xlabel('step')
                plt.subplot(3, 2, 3)
                plt.plot(sampler.chain[:, :, 1].T, alpha=0.05, color='k')
                plt.ylabel(r'$mb$') ; plt.xlabel('step')
                plt.subplot(3, 2, 4)
                plt.plot(sampler.chain[:, :, 4].T, alpha=0.05, color='k')
                plt.ylabel(r'$md$') ; plt.xlabel('step')
                plt.subplot(3, 2, 5)
                plt.plot(sampler.chain[:, :, 2].T, alpha=0.05, color='k')
                plt.ylabel(r'$sigb$') ; plt.xlabel('step')
                plt.subplot(3, 2, 6)
                plt.plot(sampler.chain[:, :, 5].T, alpha=0.05, color='k')
                plt.ylabel(r'$sigd$') ; plt.xlabel('step')
                plt.show()

            # Remove the burn-in (need to understand this better)
            samples = sampler.chain[:, 400:, :].reshape((-1, ndim))

            # Find the 1-sigma confidence interval ()
            fits = np.percentile(samples, [16, 50, 84], axis=0)
            fits_pm = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*fits))
            best = fits[1][:]
            print 'Best parameters =',best

            # pick 40 random parameter sets from the final sample
            idx = np.random.randint(samples.shape[0], size=40)
            #fig = plt.figure(2, figsize=(8., 9.))
            fig = plt.figure(figsize=(8., 9.))
            xtit = "${\\rm log}_{10}({\\rm M_{halo}}/M_{\odot}h^{-1})$"
            ytit = "${\\rm log}_{10}(\\langle N\\rangle _{\\rm [OII]})$"
            xmin = 10. ; xmax = 15.
            ymin = -3. ; ymax = 1.
            ax = plt.subplot()
            ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
            ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
            ax.text(xmin+(xmax-xmin)*0.03, ymax-(ymax-ymin)*0.04,\
                        survey+' cuts, '+zleg)

            model = funcm(x,p0)
            #plt.plot(x,model,'r',label='Initial guess')
            model = funcm(x,best)
            plotbest = True
            if plotbest:
                #plt.plot(x,model,'r',label='Best fit')
                for i in idx:
                    # plot a sample of best fit lines
                    f = samples[i]
                    model = funcm(x,f)
                    plt.plot(x, model, color='grey', alpha=0.1)

            #Plot input
            ind = np.where(cen>0.)
            print logM[ind],np.log10(cen[ind])
            plt.plot(logM[ind],np.log10(cen[ind]),color=cols[0],\
                         label='All centrals', linewidth=5)

            ind = np.where(cs>0.)
            plt.plot(logM[ind],np.log10(cs[ind]),color=cols[0],\
                         label='Spheroid centrals', linestyle='-')

            ind = np.where(cd>0.)
            plt.plot(logM[ind],np.log10(cd[ind]),color=cols[0],\
                         label='Disk centrals', linestyle='--')
           
                
            # Legend
            leg = plt.legend(loc=1)
            leg.draw_frame(False)
            plt.show()
    ##########################
    # Save figure
    fig.savefig(plotfile)
    print plotfile
    print infile

