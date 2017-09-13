#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from Cosmology import *
import scipy.ndimage
from distinct_colours import get_distinct
import mpl_style
from astropy.io import fits
plt.style.use(mpl_style.style1)

dir = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'

model = 'MillGas/gp15newmg/' #'MillGas/gp14/'

nvol = 64

#zlist = [61, 44, 42, 40, 37] #z=0, 0.6, 0.75, 0.9, 1.18
zlist = [44, 42, 40, 37, 34] #z=0.6, 0.75, 0.9, 1.18, 1.5
cols = get_distinct(len(zlist)) 
cols.insert(0,'k')

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/colours/'+model
plotfile = outdir+line+'_decam_desi_mag.pdf'

fluxcut = 8.*10.**-17.

#############################
# Prepare plots
fig = plt.figure(figsize=(8.5,9.))

plt.ylabel("$(r-z)_{\\rm DECam}$")
plt.xlabel("$(g-r)_{\\rm DECam}$")

# Prepare histograms for contours
dxy = 0.05
xmin = -0.7 ; xmax = 1.7
ymin = 0.0 ; ymax = 1.6
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)

xbins = np.arange(xmin,xmax,dxy)
xhist = xbins +dxy*0.5
ybins = np.arange(ymin,ymax,dxy)
yhist = ybins +dxy*0.5

zhist = np.zeros(shape=(len(zlist),len(xhist),len(yhist)))
############################
vols = np.zeros(shape=(len(zlist)))
zleg = [] ; lines=[]

# Observations

# Stars
obspath = '/gpfs/data/violeta/lines/desi_hod_o2/xi_obs_data/'

fil = obspath+'lensing14_stars_decam_match.fits' 
hdulist = fits.open(fil) ; tbdata = hdulist[1].data
x = tbdata['g_decam'] - tbdata['r_decam']
y = tbdata['r_decam'] - tbdata['z_decam']
hdulist.close()
plt.scatter(x,y,color='grey',marker=(5, 1),label='Stars')


# DESI
x = [0.2,0.70,-0.3,-0.3,0.2]
y = [0.3,0.75, 1.6, 0.3,0.3]
cs, = plt.plot(x,y,color=cols[0],linewidth=3.5)
zleg.append('DESI cuts')
lines.append(cs)

# Model
for ii,iz in enumerate(zlist):
    firstpass = True
    for ivol in range(nvol):
        fil = dir+model+'/iz'+str(iz)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (os.path.isfile(fil)):
            f = h5py.File(fil,'r')
            zz = f['Output001/redshift'].value
            h0 = f['Parameters/h0'].value 
            omega0 = f['Parameters/omega0'].value
            omegab = f['Parameters/omegab'].value
            lambda0 = f['Parameters/lambda0'].value 
            vol1 = f['Parameters/volume'].value 
            f.close()

            set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                          h0=h0, universe="Flat",include_radiation=False)
            lcut = emission_line_luminosity(fluxcut,zz)
            tomag = band_corrected_distance_modulus(zz) 

            if(firstpass):
                zstring = "{:.2f}".format(zz) 
                zleg.append('z='+zstring)
                firstpass = False
            vols[ii] = vols[ii] + vol1
        else:
            print ' Not found: ',fil

        efil = dir+model+'/iz'+str(iz)+'/ivol'+str(ivol)+'/elgs.hdf5'
        if (os.path.isfile(efil)):
            f = h5py.File(efil,'r') ; g = f['Output001']

            #lum = g['L_tot_'+line].value
            #gr = g['magDgo_tot'].value - g['magDro_tot'].value
            #rz = g['magDro_tot'].value - g['magDzo_tot'].value 
            #ind = np.where(lum >lcut)

            lum_ext = g['L_tot_'+line+'_ext'].value
            r = g['magDro_tot_ext'].value + tomag
            gr = g['magDgo_tot_ext'].value - g['magDro_tot_ext'].value
            rz = g['magDro_tot_ext'].value - g['magDzo_tot_ext'].value 
            ind = np.where((lum_ext >lcut) & (r<23.4))
            
            f.close()

            xcol = gr[ind] ; ycol = rz[ind]
            #xcol = gr ; ycol = rz

            H,xedges,yedges = np.histogram2d(xcol,ycol,bins=[np.append(xbins,xmax),np.append(ybins,ymax)])
            zhist[ii,:,:] = zhist[ii,:,:] + H

print 'Start plotting'
# Plot the model predictions
for ii in range(len(zlist)):
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = np.zeros(shape=(len(X),len(Y))) ; Z.fill(-999.)
    al = [-4.0,-1.5,1.] 
    if (vols[ii]>0.):  
        zzhist = zhist[ii,:,:]
        ind = np.where(zzhist >0.)
        Z[ind] = np.log10(zzhist[ind]/vols[ii]/dxy/dxy)  # In Mpc^3/h^3
        cs = plt.contour(X, Y, Z.T, levels=al, \
                            linestyles='-',colors=cols[ii+1],\
                             label=zleg[ii])
        lines.append(cs.collections[0])


# Legend
leg = plt.legend(lines,zleg,loc=2, handlelength=0, handletextpad=0)
for item in leg.legendHandles:
    item.set_visible(False)
for color,text in zip(cols,leg.get_texts()):
    print text
    text.set_color(color)
    leg.draw_frame(False)        


#plt.show()
print 'Plot: ',plotfile
fig.savefig(plotfile)
