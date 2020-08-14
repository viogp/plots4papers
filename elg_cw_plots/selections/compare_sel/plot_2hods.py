#! /usr/bin/env python

import numpy as np
import os.path, sys
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import read_jc_obs as jc
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

model = 'gp19/'

sn = '41' ; zz = '0.83'
surveys = ['All','VVDS-DEEP']
nd = '-2.0'
cuts = ['m','sfr','lo2']

#############################
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
hodpath = inpath+'hod/'

plotfile = inpath+'plots/'+model+'selections/hod/hods_All_VVDS.pdf'
############################################
# Initialize the parameters for the figures
fig = plt.figure(figsize=(6.5,14.))

xtit = "${\\rm log}_{10}(M_{\\rm halo}/M_{\odot}h^{-1})$"
ytit = "${\\rm log}_{10}\\langle N_M\\rangle$"

xmin = 11. ; xmax = 15.
ymin = -2.9 ; ymax = 2.

cols = ['darkred','dodgerblue','palegreen']
cutlabel = ['Mass cut','SFR cut','L[OII] cut']

# Loop over the redshifts of interest
jj = 410
for iis, survey in enumerate(surveys):
    jj = jj + 1
    #Plot
    stext = survey+'; z='+zz+'; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
    if (iis == 0):
        ax1 = fig.add_subplot(jj) ; ax1.set_autoscale_on(False)
        ax1.set_xlim([xmin,xmax]) ; ax1.set_ylim([ymin,ymax])
        ax1.set_autoscale_on(False) ;  ax1.minorticks_on()
        ax1.set_ylabel(ytit) 
        ax1.tick_params(labelbottom=False)
        ax1.text(xmin+0.49*(xmax-xmin), ymin+0.04*(ymax-ymin), stext, fontsize='small')
    else:
        ax = fig.add_subplot(jj,sharex=ax1,sharey=ax1)
        ax.set_autoscale_on(False) ;  ax.minorticks_on()
        ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
        ax.text(xmin+0.32*(xmax-xmin), ymin+0.04*(ymax-ymin), stext,fontsize='small')

    for ic, cut in enumerate(cuts):
        hfile = hodpath+model+cut+'cut_'+survey+'_nd'+nd+'_sn'+sn+'.dat'
        # Check if the file exists
        if (not os.path.isfile(hfile)):
            if verbose:
                print('WARNING: {} not found'.format(infile))
            continue
        mh,nall,ncen,nsat = np.loadtxt(hfile,unpack='True')

        # All
        ind = np.where(nall>-999.)
        if (iis == 0):
            ax1.plot(mh[ind],nall[ind],color=cols[ic],linestyle='-',
                               linewidth=2.5,label=cutlabel[ic])
        else:
            ax.plot(mh[ind],nall[ind],color=cols[ic],linestyle='-',
                              linewidth=2.5,label=cutlabel[ic])

        # Centrals
        ind = np.where(ncen>-999.)
        if (iis == 0):
            ax1.plot(mh[ind],ncen[ind],color=cols[ic],linestyle='--')
        else:
            ax.plot(mh[ind],ncen[ind],color=cols[ic],linestyle='--')

        # Satellites
        ind = np.where(nsat>-999.)
        if (iis == 0):
            ax1.plot(mh[ind],nsat[ind],color=cols[ic],linestyle=':')
        else:
            ax.plot(mh[ind],nsat[ind],color=cols[ic],linestyle=':')

plt.subplots_adjust(wspace=0, hspace=0)
leg = ax1.legend(loc=2, fontsize='small',
                          handlelength=0, handletextpad=0)
for item in leg.legendHandles:
    item.set_visible(False)
for color,text in zip(cols,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)
            
# Save figures
fig.savefig(plotfile)
print('Output: {}'.format(plotfile))
