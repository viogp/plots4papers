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

dir = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/' 
nvol = 64
model = 'gp19/'

zlist = [61, 44, 41, 39, 34] #z=0.6, 0.83, 1., 1.5
cols = get_distinct(len(zlist)) 

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/'+model
plotfile = outdir+line+'_decam.pdf'

fluxcut = 8.*10.**-17.

#############################
# Prepare plots
fig = plt.figure(figsize=(8.5,9.))

plt.xlabel("$(r-z)_{\\rm DECam}$")
plt.ylabel("$(g-r)_{\\rm DECam}$")

# Prepare histograms for contours
dxy = 0.05
xmin = 0. ; xmax = 2.
ymin = -0.3 ; ymax = 1.45
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
obspath = '/cosma5/data/durham/violeta/lines/desi_hod_o2/xi_obs_data/'

fil = obspath+'lensing14_stars_decam_match.fits' 
hdulist = fits.open(fil) ; tbdata = hdulist[1].data
x = tbdata['r_decam'] - tbdata['z_decam']
y = tbdata['g_decam'] - tbdata['r_decam']
hdulist.close()
plt.scatter(x,y,color='grey',marker=(5, 1))#,label='Stars')

# eBOSS SGC
rz1 = (0.571+0.218*0.457)/(1.+0.218*0.068)
gr1 = -0.068*rz1 +0.457
rz2 = (0.571+0.218*0.773)/(1.-0.218*0.112)
gr2 = 0.112*rz2 + 0.773
rz3 = (1.901 - 0.555*0.773)/(1. + 0.555*0.112)
gr3 = 0.112*rz3 + 0.773
rz4 = (1.901 - 0.555*0.457)/(1. - 0.555*0.068)
gr4 = -0.068*rz4 + 0.457
x = [rz1,rz2,rz3,rz4,rz1]
y = [gr1,gr2,gr3,gr4,gr1]
cs, = plt.plot(x,y,'k',linewidth=3.5,label='eBOSS SGC')
lines.append(cs)

# DESI+16
rz1 = 0.3
gr1 = -0.3
rz2 = 0.3
gr2 = 0.2
rz3 = 0.75
gr3 = 0.7
rz4 = 1.6
gr4 = -0.3
x = [rz1,rz2,rz3,rz4,rz1]
y = [gr1,gr2,gr3,gr4,gr1]
cs, = plt.plot(x,y,'k',linewidth=3.5,linestyle='--',label='DESI')
lines.append(cs)

# Model
for ii,iz in enumerate(zlist):
    firstpass = True 
    for ivol in range(nvol):
        fil = dir+model+'/iz'+str(iz)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (not os.path.isfile(fil)):
            print('Not found {}'.format(fil))
        else:
            f = h5py.File(fil,'r')
            try:
                zz = f['Output001/redshift'].value
                h0 = f['Parameters/h0'].value 
                omega0 = f['Parameters/omega0'].value
                omegab = f['Parameters/omegab'].value
                lambda0 = f['Parameters/lambda0'].value 
                vol1 = f['Parameters/volume'].value 
            except:
                print('Problem with {}'.format(fil))
                continue
            f.close()

            set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                          h0=h0, universe="Flat",include_radiation=False)
            lcut = emission_line_luminosity(fluxcut,zz)

            if(firstpass):
                zstring = "{:.2f}".format(zz) 
                zleg.append('z='+zstring)
                firstpass = False
            vols[ii] = vols[ii] + vol1
            set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                              h0=h0, universe="Flat",include_radiation=False)
            inleg = 'z='+str(zz)


            efil = dir+model+'/iz'+str(iz)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efil)):
                f = h5py.File(efil,'r') ; g = f['Output001']

                gr = g['mag_DES-g_o_tot_ext'].value -\
                     g['mag_DES-r_o_tot_ext'].value
                rz = g['mag_DES-r_o_tot_ext'].value -\
                     g['mag_DES-z_o_tot_ext'].value 

                lum_ext = g['L_tot_'+line+'_ext'].value
                ind = np.where(lum_ext >lcut)
                #gmag = g['mag_DES-g_o_tot_ext'].value 
                #ind = np.where((lum_ext >lcut) & \
                    #                   (gmag>21.825) & \
                    #                   (gmag<22.825))            
                f.close()

            ycol = gr[ind] ; xcol = rz[ind]

            H,xedges,yedges = np.histogram2d(xcol,ycol,bins=[np.append(xbins,xmax),np.append(ybins,ymax)])
            zhist[ii,:,:] = zhist[ii,:,:] + H

# Plot the model predictions
for ii in range(len(zlist)):
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = np.zeros(shape=(len(X),len(Y))) ; Z.fill(-999.)
    al = [-4.5,-1.5,1.] 
    if (vols[ii]>0.):  
        zzhist = zhist[ii,:,:]
        ind = np.where(zzhist >0.)
        Z[ind] = np.log10(zzhist[ind]/vols[ii]/dxy/dxy)  # In Mpc^3/h^3
        cs = plt.contour(X, Y, Z.T, levels=al, \
                         linestyles='-',colors=cols[ii])
        lines.append(cs.collections[0])

# Legend
leg = plt.legend(lines,zleg,loc=4, handlelength=0, handletextpad=0)
for item in leg.legendHandles:
    item.set_visible(False)
for color,text in zip(cols,leg.get_texts()):
    #print text
    text.set_color(color)
    leg.draw_frame(False)        

#plt.show()
fig.savefig(plotfile)
print 'Plot: ',plotfile
