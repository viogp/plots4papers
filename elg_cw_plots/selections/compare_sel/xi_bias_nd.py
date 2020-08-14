#! /usr/bin/env python

import numpy as np
import os.path, sys
from scipy import interpolate
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

def chi2(obs,model,err,nparams):
    val = 0.
    for i,iobs in enumerate(obs):
        val = val + ((iobs-model[i])/err[i])**2

    dof = len(obs) - nparams

    return val/dof
##################################

space = 'r' #'r' or 'z'

model = 'gp19/'

# Paths 
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
plotpath = inpath+'plots/'+model+'selections/xi/'
filepath = inpath+'selections/'

snapnum = ['39','41']
surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']

# Min. number of galaxy pairs to be considered
npairs = 1.
##################################

# Initialize the parameters for the figures
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15

xmin = -2. ; xmax = 2.
ytitb = "$\sqrt{\\xi_{gg}/\\xi_{DM}}$"
yminb = -1. ; ymaxb = 3.
ymin = -2.9 ; ymax = 5.

if (space == 'r'):
    xtit = "${\\rm log}_{10}(\\rm{r/Mpc}\, h^{-1})$"
    ytit = "${\\rm log}_{10}\\xi (\\rm{r})$" 
else:
    xtit = "${\\rm log}_{10}(\\rm{s/Mpc}\, h^{-1})$"
    ytit = "${\\rm log}_{10}\\xi (\\rm{s})$" 

cols = ['darkred','dodgerblue','palegreen']
linestyle = ['-','--',':','.-']
##################################
# Initialize arrays to minimise chi2 for the bias
abias = np.linspace(0.1,10.,10000)
chis = np.zeros((len(abias))) ; chis.fill(999.)
rl = 8.
rh = 50.

# Loop over the different files
for ii,iz in enumerate(snapnum):
    epath = filepath+model+'iz'+iz+'/'+space+'/'
    # File to output bias values
    bfile = epath+'bias_sn'+str(iz)+'.dat'

    with open(bfile, 'w') as outf: 
        outf.write('# bias({}-{} Mpc/h), log(bias), log(nd/Mpc-3h3), Cut, Survey \n'.format(rl,rh))

    if (space == 'r'):
        # DM
        pathDM = '/cosma5/data/durham/violeta/mr7corr/'
        dmfile = pathDM+'xi_real_sn0'+iz+'.txt' 
        if (not os.path.isfile(dmfile)):
            print('STOP: {} not found'.format(dmfile)) ; sys.exit()
        # CUTE output: r,xi(r), error(r), DD(r)
        # the third column is the Poisson error calculated from DD.
        dm_r, dm_xi, dm_inerr, dm_dd = np.loadtxt(dmfile,unpack=True)
        
        ind = np.where((dm_r>0.) & (dm_xi>0.))
        lr_dm = np.log10(dm_r[ind])
        y_dm  = np.log10(dm_xi[ind])
        eg_dm = np.log10(np.exp(1.))*dm_inerr[ind]/dm_xi[ind]

    survey1 = surveys1[ii]
    survey2 = surveys2[ii]

    inleg = ['Mass cut, All','Mass cut, '+survey1,'Mass cut, '+survey2,
             'SFR cut, All','SFR cut, '+survey1,'SFR cut, '+survey2,
             'L[OII] cut, All','L[OII] cut, '+survey1,'L[OII] cut, '+survey2]

    for nd in ['-2.0','-3.0','-4.2']:
        # Initialize the parameters for the figures
        fig = plt.figure(figsize=(8.,9.))
        gs = gridspec.GridSpec(4,1)
        gs.update(wspace=0., hspace=0.)

        # Bias
        axb = plt.subplot(gs[3,:])
        axb.set_autoscale_on(False) ; axb.minorticks_on()
        axb.set_xlim(xmin,xmax) ; axb.set_ylim(yminb,ymaxb)
        axb.set_xlabel(xtit, fontsize=fs)
        axb.set_ylabel(ytitb, fontsize=fs)

        # Plot 2PCF r-space
        ax = plt.subplot(gs[:-1,:],sharex=axb)
        ax.set_ylabel(ytit, fontsize=fs)
        ax.set_autoscale_on(False) ; ax.minorticks_on()
        ax.set_ylim(ymin,ymax) ; start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end, 1.))
        if(iz == '41'):
            ztext = 'z = 0.83; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
        elif(iz == '39'):
            ztext = 'z = 0.99; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
        ax.text(xmax-0.5*(xmax-xmin),ymax-0.07*(ymax-ymin), ztext)
        plt.setp(ax.get_xticklabels(), visible=False)

        # DM
        ax.errorbar(lr_dm,y_dm,yerr=eg_dm,color='k')

        ii = -1
        for ic, cut in enumerate(['m','sfr','lo2']):
            surv = ['All',survey1,survey2]

            for iif, ss in enumerate(surv):
                ii += 1

                efile = epath+cut+'cut_'+ss+\
                        '_nd'+nd+'_sn'+iz+'_CUTExi_'+space+'.dat'
                if (not os.path.isfile(efile)):
                    continue

                # r   xi(r)   error(r)   DD(r)
                in_r,in_xi,in_err,in_dd = np.loadtxt(efile, unpack=True)

                ind = np.where((in_dd >npairs) & (in_r>=10**xmin))
                if (np.shape(ind)[1] > 0):
                    xg = in_r[ind] ; yg = in_xi[ind]
                    eg = in_err[ind]

                    # 2PCF
                    ie = np.where((xg>0.) & (yg>0.))
                    x  = np.log10(xg[ie])
                    y  = np.log10(yg[ie])
                    lerr =  np.log10(np.exp(1.))*eg[ie]/yg[ie]

                    ax.errorbar(x,y,yerr=lerr,color=cols[ic],
                                    linestyle=linestyle[iif],label=inleg[ii])

                    # Interpolate DM to x-galaxies
                    ydm = np.interp(xg,dm_r,dm_xi)
                    ddm = np.interp(xg,dm_r,dm_dd)
                    derr = (1.+ydm)/np.sqrt(ddm)

                    # Bias
                    ie = np.where((ydm > 0.)  & (yg > 0.))
                    rr = xg[ie]
                    xig = yg[ie]  ; exig = eg[ie]
                    xid = ydm[ie] ; exid = derr[ie]

                    x_bias = np.log10(rr)
                    y_bias = np.sqrt(xig/xid)
                    berr = 0.5*y_bias*np.sqrt(exig**2/xig**2 +
                                              exid**2/xid**2)

                    axb.errorbar(x_bias,y_bias, yerr=berr, color=cols[ic],
                                 linestyle=linestyle[iif])
                    
                    # Estimate bias in the range rl-rh
                    rind = np.where((rr>=rl) & (rr<=rh))
                    oxi = xig[rind]
                    oerr = exig[rind]
                    for ib,bb in enumerate(abias):
                        dmh = xid[rind]*bb*bb
                        chis[ib] = chi2(oxi,dmh,oerr,1)

                    chi2min = np.nanmin(chis)
                    ib = np.where(chis == chi2min)
                    bias = abias[ib][0]

                    # Error
                    ind = np.argsort(chis)
                    b1 = np.interp(chi2min+1,chis[ind],abias[ind])
                    ebb = abs(b1-bias)

                    # Write bias to file
                    #b8 = np.interp(8.,xg[ie],y_bias)
                    lbb = np.log10(bias)
                    with open(bfile, 'a') as outf:
                        outf.write('{:.3f} $\pm${:.3f}  {:.2f}  {}  {}  {} \n'.format(bias,ebb,lbb,nd,cut,ss))

        # Legends
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] for h in handles]
        leg = ax.legend(handles, labels, loc=3, fontsize = fs-2)
        leg.draw_frame(False)
    
        # Save figure
        newnd = nd.replace('.','p')
        plotfile = plotpath+'xi_'+space+\
                   '_'+survey1+'_nd'+newnd+'_sn'+iz+'_env2.pdf'
        fig.savefig(plotfile)
        print('Output: {}'.format(plotfile))
    # Bias file
    print('Bias file: {}'.format(bfile))
