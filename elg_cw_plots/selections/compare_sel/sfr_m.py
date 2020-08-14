#! /usr/bin/env python

import numpy as np
import os.path, sys
import subprocess
from Cosmology import *
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import ndimage
from stats import *
import mpl_style
plt.style.use(mpl_style.style1)

# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
#-----------------------------------------------------

Testing = False

volume = (500.**3.) ; verbose = True

model = 'gp19/'

sn_list = ['39','41']
zz_list = ['0.99','0.83']
surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']
nds = np.array([-2.0,-3.0,-4.2])
cuts = ['lo2','sfr','m']

if Testing:
    sn_list = ['39']

#################################################
# Prepare the plot
#al = np.sort(nds - 0.5) #; print(al)
al = np.sort(nds) 
cols=['navy','royalblue','lightsteelblue']
lsty=['-','--',':']

# Paths
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
file_path = inpath+'selections/'+model+'ascii_files/'
plt_path = inpath+'plots/'+model+'selections/sfr_m/'

# Initialize GSMF
mmin = 8.5
mmax = 15.
dm = 0.1
mbins = np.arange(mmin,mmax,dm)
mhist = mbins + dm*0.5
# Initialize SSFR
smin = 3.
smax = 13.
ds = 0.1
sbins = np.arange(smin,smax,ds)
shist = sbins + ds*0.5

nlines = len(nds)*len(cuts)

# Read the information from the different files
for iiz,sn in enumerate(sn_list): 
    surveys = ['All',surveys1[iiz], surveys2[iiz]]
    for survey in surveys:
        # Figure http://matplotlib.org/users/gridspec.html
        mtext = survey+'_sn'+sn
        plt_file = plt_path+mtext+'.pdf'

        fig = plt.figure(figsize=(8.5,9.))
        gs = gridspec.GridSpec(3, 3)
        gs.update(wspace=0., hspace=0.)
        ax = plt.subplot(gs[1:,:-1])

        # Fig. SFRF vs M
        xtit="$log_{10}(\\rm M_{*}/M_{\odot}h^{-1})$"
        ytit="$log_{10}(\\rm SFR/M_{\odot}h^{-1}Gyr^{-1})$"
        xmin=mmin ; xmax=11.9 ; ymin=6. ; ymax=smax
        ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
        ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)

        # GSMF 
        axm = plt.subplot(gs[0, :-1],sharex=ax)
        ytit="$log_{10}(\Phi)$" ; axm.set_ylabel(ytit)
        axm.set_autoscale_on(False) ;  axm.minorticks_on()
        axm.set_ylim(-5.5,-1.) 
        plt.setp(axm.get_xticklabels(), visible=False)


        # SFRF
        axs = plt.subplot(gs[1:, 2],sharey=ax)
        xtit="$log_{10}(\Phi)$" ; axs.set_xlabel(xtit)
        axs.set_autoscale_on(False) ;  axs.minorticks_on()
        axs.set_xlim(-4.4,0.0)
        start, end = axs.get_xlim()
        axs.xaxis.set_ticks(np.arange(-4., end, 1.))
        plt.setp(axs.get_yticklabels(), visible=False)

        # Read the files
        nfiles = 0 
        for iic,cut in enumerate(cuts):
            for iin,nd in enumerate(nds):
                # Check the existance of the files
                mtext = cut+'cut_'+survey+'_nd'+str(nd)+'_sn'+sn
                infile = file_path+mtext+'.dat'

                if (not os.path.isfile(infile)):
                    if verbose:
                        print('WARNING: {} not found'.format(infile))
                    continue
                # Jump files with only a one line header
                wcl_line = subprocess.check_output(["wc", "-l",infile])
                wcl = int(wcl_line.split()[0])
                if (wcl <= 1):
                    if verbose:
                        print('WARNING: {} has too few lines'.format(infile))
                    continue
                nfiles += 1 

                lmass,lsfr = np.loadtxt(infile,usecols=(7,8),
                                        unpack=True)

                # GSMF
                H, bins_edges = np.histogram(lmass,bins=np.append(mbins,mmax))
                gsmf = H/volume/dm  # In Mpc^3/h^3

                # sSFR
                H, bins_edges = np.histogram(lsfr,bins=np.append(sbins,smax))
                sfrf = H/volume/ds

                # sSFR-GSMF
                H, xedges, yedges = np.histogram2d(lsfr,lmass,
                                                   bins=[np.append(sbins,smax),
                                                         np.append(mbins,mmax)])
                smf = H/volume/dm/ds  

                # Plot SMF vs SFR
                matplotlib.rcParams['contour.negative_linestyle'] = lsty[iic]
                zz = np.zeros(shape=(len(shist),len(mhist))) ; zz.fill(-999.) 
                ind = np.where(smf>0.)
                zz[ind] = np.log10(smf[ind]) 

                ind = np.where(zz>-999.)
                if (np.shape(ind)[1]>1):
                    # Contours
                    xx,yy = np.meshgrid(mbins,sbins)
                    #al = nds[iin] ; print(al)
                    cs = ax.contour(xx, yy, zz, levels=al, \
                                    colors=cols[iin])
                    #cs.levels = [nf(val) for val in cs.levels]
                    #ax.clabel(cs, cs.levels, inline=1,inline_spacing=0,\
                        #  fontsize=10,fmt='%r')#fmt='%r %%')

                # GSMF 
                py = gsmf ; ind = np.where(py>0.)
                x = mhist[ind] ; y = np.log10(py[ind])
                ind = np.where(y < 0.)
                axm.plot(x[ind],y[ind],color=cols[iin],\
                         linestyle=lsty[iic])

                # SFRF
                px = sfrf ; ind = np.where(px>0.)
                y = shist[ind] ; x = np.log10(px[ind])
                ind = np.where(x < 0.)
                if (iic == 1):
                    inleg = '$n_{\\rm gal}=10^{'+str(nds[iin])+'}{\\rm Mpc}^{-3}h^{-3}$'
                    axs.plot(x[ind],y[ind],color=cols[iin],\
                             linestyle=lsty[iic],label=inleg)
                else:
                    if (iin == 2 and iic == 0):
                        axs.plot([],[],' ',
                                 label=survey+', z='+zz_list[iiz])

                    axs.plot(x[ind],y[ind],color=cols[iin],\
                             linestyle=lsty[iic])

        if (nfiles>0):
            # Legend
            leg = axs.legend(bbox_to_anchor=(1.05, 1.4),fontsize='small', \
                             handlelength=0, handletextpad=0)
            for item in leg.legendHandles:
                item.set_visible(False)
            allcols = ['k'] + cols
            for color,text in zip(allcols,leg.get_texts()):
                text.set_color(color)
                leg.draw_frame(False)

            # Save figures
            fig.savefig(plt_file)
            print('Output: ',plt_file)

