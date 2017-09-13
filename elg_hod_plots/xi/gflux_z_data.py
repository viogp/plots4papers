#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import string
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
    iz = '0.91' #'0.76' 
#############################

model = 'MillGas/gp15newmg/'

npairs = 1.
plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/xi/'
plotfile = plotdir+model+ 'zspace_z'+iz+'_gflux.pdf'

#############################
# Surveys
#surveys = ['4cute_g22.0_z.dat','4cute_g22.2_z.dat','4cute_g22.5_z.dat','gcut_4cute_gflux_z.dat']
#surveys = ['fm16_4cute_gflux_z.dat','5fm16_4cute_gflux_z.dat','fm17_4cute_gflux_z.dat','5fm17_4cute_gflux_z.dat','fm18_4cute_gflux_z.dat'] 
#surveys = ['fm17_CUTExi_gflux_z.dat','f5m17_CUTExi_gflux_z.dat','fm16_CUTExi_gflux_z.dat','f5m16_CUTExi_gflux_z.dat','gcut_CUTExi_gflux_z.dat']
#inleg = surveys

surveys=['fm16_CUTExi_gflux_z.dat','fm18_CUTExi_gflux_z.dat','gcut_CUTExi_gflux_z.dat','4cute_g22.2_z.dat']
inleg = ['$F_{[OII]}>10^{-16}erg\,s^{-1}cm^{-2}$','$F_{[OII]}>10^{-18}erg\,s^{-1}cm^{-2}$','20<g<22.8 (Favole cut)','20<g<22.2']

sn = ['44','42','40','37','34']
zz = ['0.62','0.76','0.91','1.17','1.50']

isn = zz.index(iz) 
szz = 'z='+zz[isn]

cols = get_distinct(len(surveys))
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
obspath = '/gpfs/data/violeta/lines/desi_hod_o2/xi_obs_data/'
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
ax.errorbar(ox,oy, yerr=oerr, color=cols[0], ecolor=cols[0],\
                label ='Favole et al. 2016', fmt = 'o')

# Galaxy predictionsmv 
pathgal = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'+model
for i,s in enumerate(surveys):
    infil = pathgal+'iz'+sn[isn]+'/z/OII3727_'+s
    print 'Reading ',infil
    if(os.path.isfile(infil)):
        # r   xi(r)   error(r)   DD(r)
        in_r,in_xi,in_err,in_dd = np.loadtxt(infil, unpack=True)
    
        ind = np.where((in_dd >npairs) & (in_r>=xmin))
    
        if (np.shape(ind)[1] > 0):
            x = in_r[ind] 
            y = x*in_xi[ind]
            err = x*(1.+in_xi[ind])/np.sqrt(in_dd[ind])
        
            ax.errorbar(x,y,yerr=err,color=cols[i+1],label =inleg[i])

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

#7_fm16_4cute_gflux_z.dat
#d Chi2 4.69739647227
#7_fm18_4cute_gflux_z.dat
#d Chi2 3.45560440144
#7_g22.2_4cute_gflux_z.dat
#d Chi2 2.67008089966
#7_gcut_4cute_gflux_z.dat
#d Chi2 3.1559949062
#/gpfs/data/violeta/lines/desi_hod_o2/xi/MillGas/gp15newmg/zspace_z0.9_gflux.pdf
