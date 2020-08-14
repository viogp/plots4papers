#! /usr/bin/env python

import numpy as np
import os.path, sys
from scipy import interpolate
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

#############################
#
#  Input ARGUMENTS
#
narg = len(sys.argv)
if(narg == 2):
    survey = sys.argv[1]
else:
    survey = 'DEEP2'
#############################

#model = 'MillGas/gp15newmg/'
#pathgal = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'+model

#model = 'MillGas/gp17/'
model = 'MillGas/gp17.spin.ramt0.01/'
pathgal = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/'+model

npairs = 1.
plotdir = '/cosma5/data/durham/violeta/lines/desi_hod_o2/xi/'
plotfile = plotdir+model+ 'xi_bias_'+survey+'.pdf'

#############################

# Surveys
surveys = ['DEEP2','VVDSDEEP','VVDSWIDE','eBOSS','DESI']
sn =['44','42','40','37']#,'34']
zz = ['0.62','0.76','0.91','1.17','1.5']

iss = surveys.index(survey) 

cols = get_distinct(len(sn))

# Plots
fig = plt.figure(figsize=(8.,9.))
gs = gridspec.GridSpec(4,1)
gs.update(wspace=0., hspace=0.)
xtit = "${\\rm log}_{10}(\\rm{r/Mpc}\, h^{-1})$"
xmin = -2. ; xmax = 2.

# Plot bias
ymin = 0.2 ; ymax = 1.9
ytit = "$\sqrt{\\xi_{gg}/\\xi_{DM}}$"
axb = plt.subplot(gs[3,:])
axb.set_autoscale_on(False) ; axb.minorticks_on()
axb.set_xlim(xmin,xmax) ; axb.set_ylim(ymin,ymax)
axb.set_xlabel(xtit) ; axb.set_ylabel(ytit)

# Plot 2PCF r-space
ymin = -2.9 ; ymax = 4.
ytit = "${\\rm log}_{10}\\xi (\\rm{r})$" 
ax = plt.subplot(gs[:-1,:],sharex=axb) 
ax.set_ylabel(ytit)
ax.set_autoscale_on(False) ; ax.minorticks_on()
ax.set_ylim(ymin,ymax) ; start, end = ax.get_xlim()
ax.xaxis.set_ticks(np.arange(start, end, 1.))
ax.text(xmax-0.25*(xmax-xmin),ymax-0.07*(ymax-ymin), survey)
plt.setp(ax.get_xticklabels(), visible=False)

# Galaxy predictions
for i,isn in enumerate(sn):
    # DM
    pathDM = '/cosma5/data/durham/violeta/lines/desi_hod_o2/dm_w7_xi/' 
    dmfile = pathDM+'xi_real_sn0'+sn[i]+'.txt' 

    if os.path.isfile(dmfile):
        # CUTE output: r,xi(r), error(r), DD(r)
        # the third column is the Poisson error calculated from DD.
        dm_r, dm_xi, dm_inerr, dm_dd = np.loadtxt(dmfile,unpack=True)

        ind = np.where((dm_r>0.) & (dm_xi>0.))
        lr_dm  = np.log10(dm_r[ind])
        y_dm   = np.log10(dm_xi[ind])
        eg_dm = ((dm_xi[ind]+1.)/np.sqrt(dm_dd[ind]))\
            /(np.log(10.)*dm_xi[ind])
        #ax.errorbar(lr_dm,y_dm,yerr=eg_dm,color='k')

    infil = pathgal+'iz'+sn[i]+'/r/OII3727_'+survey+'_CUTExi.dat'
    if(os.path.isfile(infil)):
        # r   xi(r)   error(r)   DD(r)
        in_r,in_xi,in_err,in_dd = np.loadtxt(infil, unpack=True)

        ind = np.where((in_dd >npairs) & (in_r>=10**xmin))
        if (np.shape(ind)[1] > 0):
            xg = in_r[ind] ; yg = in_xi[ind]
            eg = (yg+1.)/np.sqrt(in_dd[ind])

            #2PCF
            ie = np.where((xg>0.) & (yg>0.))
            x  = np.log10(xg[ie])
            y  = np.log10(yg[ie])
            err = eg[ie]/(np.log(10.)*yg[ie])

            col = cols[i]      
            szz = 'z='+zz[i]
            ax.errorbar(x,y,yerr=err,color=col,label =szz)

            # Bias
            ydm = np.interp(xg,dm_r,dm_xi)
            ddm = np.interp(xg,dm_r,dm_dd)

            ie = np.where((ydm > 0.)  & (yg > 0.)) 
            x_bias = np.log10(xg[ie])
            y_bias = np.sqrt(yg[ie]/ydm[ie])
            ed = (ydm[ie]+1.)/np.sqrt(ddm[ie])
            err = np.sqrt(eg[ie]*eg[ie]/(4.*ydm[ie]*yg[ie])\
                              + ed*ed*yg[ie]/(4.*ydm[ie]**2)) 

            b = np.interp(10.,xg[ie],y_bias)
            print survey,zz[i],b,np.log10(b)

            axb.errorbar(x_bias,y_bias, yerr=err, color=col)
        else:
            print ' * Too few data for ',infil
    else:
        print ' * File not found: ', infil

# Legends
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] for h in handles]
leg = ax.legend(handles, labels, loc=3, handlelength=0, handletextpad=0)
for item in leg.legendHandles:
    item.set_visible(False)
for color,text in zip(cols,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)

# Save figure
fig.savefig(plotfile)
print 'Output: ',plotfile
