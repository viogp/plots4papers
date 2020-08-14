#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import string
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from Cosmology import * 
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

#############################
#
#  Input ARGUMENTS
#
narg = len(sys.argv)
if(narg == 2):
    iz = sys.argv[1]
else:
    iz = '0.83' #'0.99' '0.83' 
#############################
#selection = 'fm16' # fm18, g22.2, g22.8
selection = 'fm18' 
#selection = 'g22.2'
#selection = 'g22.8'

models = ['gp19','gp19.font','gp19.starvation'] 
#inleg = ['1% gas stripping','10% gas stripping ','100% gas stripping'] 
inleg = ['This work','10% stripping ','Starvation'] 
lsty=['-','--',':'] ; lwt = [3., 1.5, 1.5]

pathgal = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'

npairs = 1.

plotdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/xi/ram/'
plotfile = plotdir+selection+'_z'+iz+'_favole.pdf'

#############################

sn = ['41','39']
zz = ['0.83','0.99']

isn = zz.index(iz) 
szz = 'z='+zz[isn]+', '+selection
#szz = 'Universe Age = 6.7Gyr'

cols = get_distinct(len(models)+2)
cols.insert(0,'grey') #; cols.insert(0,'grey') 

#############################

rule = string.maketrans('D', 'E')

#############################

# Plots
fig = plt.figure(figsize=(8.,9.))

xtit = "$s(h^{-1}\\rm{Mpc})$"
ytit = "$s\\xi (\\rm{s})$"

xmin = 2. ; xmax = 30.
ymin = 0. ; ymax = 10.

ax = plt.subplot() 
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
ax.set_autoscale_on(False) ; ax.minorticks_on()
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax)
#start, end = ax.get_xlim()
#ax.xaxis.set_ticks(np.arange(start, end, 1.))
ax.text(xmin+0.05*(xmax-xmin),ymin+0.05*(ymax-ymin), szz)

# Plot Favole+15
obspath = '/cosma5/data/durham/violeta/lines/desi_hod_o2/xi_obs_data/'
ho = 0.677

fil = obspath+'SHAM_norm-mean1e+12-sig500000000000.0-fsat0.225_ELGdensity_radecz_monopole.dat'
#  s, xi(s), poisson errors, DD, DR, RR
so, xio, in_err, dd = np.loadtxt(fil,usecols=(0,1,2,3),unpack=True)
ind = np.where((dd>0.) & (xio >0.))
x   = so[ind] #*ho
y   = x*xio[ind]
err = x*in_err[ind] #err = x*(1.+xio[ind])/np.sqrt(dd[ind])
#ax.errorbar(x,y,yerr=err, color=cols[1],\
#                label ='Favole et al. 2015: model')

#  s, xis, poisson errors
fil = obspath+'Combined_W134_xi.CORR.dat'
so, xio, erro = np.loadtxt(fil, unpack=True)
ind = np.where(xio >0.)
ox   = so[ind] #*ho
oy   = so[ind]*xio[ind]
oerr = so[ind]*erro[ind] 
obslabel = 'Favole et al. 2016'
#obslabel = 'Observations'
ax.errorbar(ox,oy, yerr=oerr, color=cols[0], ecolor=cols[0],\
                label =obslabel, fmt = 'o')

# Galaxy predictionsmv 
for i,model in enumerate(models):
    infil = pathgal+model+'/iz'+sn[isn]+'/z/OII3727_'+\
        selection+'gflux_CUTExi_z.dat'

    if(not os.path.isfile(infil)):
        print 'NOT found ',infil
    else:
        print 'Reading ',infil
        # r   xi(r)   error(r)   DD(r)
        in_r,in_xi,in_err,in_dd = np.loadtxt(infil, unpack=True)
    
        ind = np.where((in_dd >npairs) & (in_r>=xmin))
    
        if (np.shape(ind)[1] > 0):
            x = in_r[ind] 
            y = x*in_xi[ind]
            err = x*(1.+in_xi[ind])/np.sqrt(in_dd[ind])
        
            ax.errorbar(x,y,yerr=err,color=cols[i+1],\
                            linestyle=lsty[i],linewidth=lwt[i], \
                            label =inleg[i])

            my = np.interp(ox,x,y)
            chi = (my-oy)**2/oerr**2
            rchi2 = np.sum(chi)/len(chi)
            print '   Reduced Chi2',rchi2

# Legends
leg = ax.legend(loc=1)
for color,text in zip(cols,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)

# Save figures
fig.savefig(plotfile)
print 'Output: ',plotfile

