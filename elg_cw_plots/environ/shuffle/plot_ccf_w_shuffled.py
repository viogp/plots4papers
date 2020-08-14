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

Testing = False

space = 'r' #'r' or 'z'

model = 'gp19/'

# Min. number of galaxy pairs to be considered
npairs = 1.

# Paths 
inpath = '/cosma6/data/dp004/dc-gonz3/lines/cosmicweb/'
plotpath = inpath+'plots/'+model+'selections/crosscorr/'
filepath = inpath+'selections/'

snapnum = ['39','41']
zzs = ['0.99','0.83']
surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']

cw_list = ['Vweb','Pweb']  

if Testing:
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI'] 
    snapnum = ['39'] ; zzs = ['0.99']
    cw_list = ['Vweb'] 

##################################

# Initialize the parameters for the figures
xmin = -2. ; xmax = 1.8
ytitb = "$\\xi_{\\rm All,X}/\\xi_{\\rm DM}$"
yminb = -3.5 ; ymaxb = 3.5
ymin = -2.9 ; ymax = 5.

if (space == 'r'):
    xtit = "${\\rm log}_{10}(r/h^{-1}{\\rm Mpc})$"
    ytit = "${\\rm log}_{10}\\xi_{\\rm All,X} (r)$" 
else:
    xtit = "${\\rm log}_{10}(s/h^{-1}{\\rm Mpc})$"
    ytit = "${\\rm log}_{10}\\xi_{\\rm All,X} (s)$" 

cut = 'elg'

elabels = ['All','Voids','Sheets','Filaments','Knots','All shuffled'] 
cols = ['k','#e7d4e8','#7fbf7b','#af8dc3','#1b7837','grey']
lstyle = ['-','-','-','-','-','-.'] 
lwidth = [4,2.5,2.5,2.5,2.5,2.5] 

##################################
# Initialize arrays to minimise chi2 for the bias
abias = np.linspace(0.1,10.,10000)
chis = np.zeros((len(abias))) ; chis.fill(999.)
rl = 8.
rh = 50.

# Loop over the different files
for cw in cw_list:
    for iiz,sn in enumerate(snapnum):
        zz = zzs[iiz] 
        epath = filepath+model+'iz'+sn+'/'+space+'/'

        if (space == 'r'):
            # DM
            pathDM = '/cosma6/data/dp004/dc-gonz3/mr7corr/'
            dmfile = pathDM+'xi_real_sn0'+sn+'.txt'
            if (not os.path.isfile(dmfile)):
                print('STOP: {} not found'.format(dmfile)) ; sys.exit()
            # CUTE output: r,xi(r), error(r), DD(r)
            # the third column is the Poisson error calculated from DD.
            dm_r, dm_xi, dm_inerr, dm_dd = np.loadtxt(dmfile,unpack=True)
            
            ind = np.where((dm_r>0.) & (dm_xi>0.))
            lr_dm = np.log10(dm_r[ind])
            y_dm  = np.log10(dm_xi[ind])
            eg_dm = np.log10(np.exp(1.))*dm_inerr[ind]/dm_xi[ind]
    
        surveys = [surveys1[iiz],surveys2[iiz]] 
        for survey in surveys:
            # Initialize the parameters for the figures
            fig = plt.figure(figsize=(6.5,9.))
            gs = gridspec.GridSpec(4,1)
            gs.update(wspace=0., hspace=0.)

            # Bias
            axb = plt.subplot(gs[3,:])
            axb.set_autoscale_on(False) ; axb.minorticks_on()
            axb.set_xlim(xmin,xmax) ; axb.set_ylim(yminb,ymaxb)
            axb.set_xlabel(xtit) ; axb.set_ylabel(ytitb)
    
            # Plot 2PCF r-space
            ax = plt.subplot(gs[:-1,:],sharex=axb)
            ax.set_ylabel(ytit)
            ax.set_autoscale_on(False) ; ax.minorticks_on()
            ax.set_ylim(ymin,ymax) ; start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 1.))
            #ax.text(xmax-0.5*(xmax-xmin),ymax-0.07*(ymax-ymin), ztext)
            plt.setp(ax.get_xticklabels(), visible=False)
    
            # DM
            ax.errorbar(lr_dm,y_dm,yerr=eg_dm,color='k',
                        linewidth=lwidth[1],linestyle=':',label='Dark Matter')
            axb.plot(lr_dm,y_dm/y_dm, color='k',
                         linewidth=lwidth[1],linestyle=':')
    
            # Loop over environments
            for iie,elabel  in enumerate(elabels):
                if (elabel == 'All shuffled'):
                    efile = epath+'shuffled_'+cut+'cut_'+survey+\
                        '_sn'+sn+'_CUTExi_'+space+'.dat'
                elif (elabel == 'All'):
                    efile = epath+cut+'cut_'+survey+\
                        '_sn'+sn+'_CUTExi_'+space+'.dat'
                    elabel = elabel+' '+survey+', z = '+zz
                else:
                    efile = epath+'ccf/'+elabel+'_'+cw+'_'+cut+'cut_'+survey+\
                        '_sn'+sn+'_ccf_'+space+'.dat'

                if (not os.path.isfile(efile)): continue
                # r   xi(r)   error(r)   DD(r)                
                in_r,in_xi,in_err,in_dd = np.loadtxt(efile, usecols=(0,1,2,3),unpack=True)

                ind = np.where((in_dd >npairs) & (in_r>=10**xmin))
                if (np.shape(ind)[1] > 0):
                    xg = in_r[ind] ; yg = in_xi[ind]
                    eg = in_err[ind]
    
                    # 2PCF
                    ie = np.where((xg>0.) & (yg>0.))
                    x  = np.log10(xg[ie])
                    y  = np.log10(yg[ie])
                    lerr =  np.log10(np.exp(1.))*eg[ie]/yg[ie]
                    ax.errorbar(x,y,yerr=lerr,color=cols[iie],
                                linewidth=lwidth[iie],linestyle=lstyle[iie],
                                label=elabel)

                    ie = np.where((xg>0.) & (yg<0.))
                    x  = np.log10(xg[ie])
                    y  = np.log10(abs(yg[ie]))
                    lerr =  np.log10(np.exp(1.))*eg[ie]/abs(yg[ie])
                    ax.errorbar(x,y,yerr=lerr,color=cols[iie],
                                linewidth=lwidth[iie],linestyle='--')
    
                    if (space == 'r'):
                        # Interpolate DM to x-galaxies
                        ydm = np.interp(xg,dm_r,dm_xi)
                        ddm = np.interp(xg,dm_r,dm_dd)
                        derr = (1.+ydm)/np.sqrt(ddm)

                        # Ratios
                        ir = np.where(ydm>0.)
                        rr = xg[ir]
                        xig = yg[ir]  ; exig = eg[ir]
                        xid = ydm[ir]  ; exid = derr[ir]

                        x_ratio = np.log10(rr)
                        y_ratio = xig/xid
                        berr = y_ratio*np.sqrt(exig**2/xig**2 +
                                               exid**2/xid**2)
                        axb.errorbar(x_ratio,y_ratio, yerr=berr, color=cols[iie],
                                     linewidth=lwidth[iie],
                                     linestyle=lstyle[iie])                    

            # Legends
            handles, labels = ax.get_legend_handles_labels()
            handles = [h[0] for h in handles]
            leg = ax.legend(handles, labels, loc=1)
            leg.draw_frame(False)
    
            # Save figure
            plotfile = plotpath+'ccf_w_shuffled_'+space+'_env'+cw+\
                       '_'+survey+'_sn'+sn+'.pdf'

            fig.savefig(plotfile)
            print('Output: {}'.format(plotfile))
