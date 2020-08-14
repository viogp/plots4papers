#! /usr/bin/env python

import numpy as np
import h5py, os.path, sys
from Cosmology import *
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from distinct_colours import get_distinct,pault_cmap
import matplotlib.gridspec as gridspec
from scipy import ndimage
from stats import *
import mpl_style
plt.style.use(mpl_style.style1) ; ptmap=pault_cmap(1)

model = 'gp19'
path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
nvol = 64

snap_list = [39,41]

plotdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/'+model+'/lo2/' 
line = 'OII3727' ; lline = '[OII]'

# Initialize prop
pmin = -0.25
pmax = 6.
dp = 0.1
pbins = np.arange(pmin,pmax,dp)
phist = pbins + dp*0.5

# Initialize LOII
lmin = 38.
lmax = 46.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

nlevel = 10

lo2prop = np.zeros(shape=(len(pbins),len(lbins)))
lo2prop_ext = np.zeros(shape=(len(pbins),len(lbins)))

###########################################
# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()


for zsnap in snap_list:
    volume = 0. ; first = True
    for ivol in range(nvol):
        gfile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if(os.path.isfile(gfile)):
            # Read the relevant information from galaxies.hdf5
            f = h5py.File(gfile,'r') #; print gfile
            zz   = f['Output001/redshift'].value
            vol1 = f['Parameters/volume'].value ; volume = volume + vol1
            h0 =   f['Parameters/h0'].value 
            omega0 = f['Parameters/omega0'].value
            omegab = f['Parameters/omegab'].value
            lambda0 =f['Parameters/lambda0'].value
            f.close()
    
            set_cosmology(omega0=omega0,omegab=omegab, \
                          lambda0=lambda0,h0=h0, \
                          universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz)
    
            efile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r')
                lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                lum = f['Output001/L_tot_'+line].value

                g_ext = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                r_ext = f['Output001/mag_DES-r_o_tot_ext'].value + tomag
                z_ext = f['Output001/mag_DES-z_o_tot_ext'].value + tomag
                rz_ext = r_ext-z_ext ; gr_ext = g_ext-r_ext
    
                g = f['Output001/mag_DES-g_o_tot'].value + tomag
                r = f['Output001/mag_DES-r_o_tot'].value + tomag
                z = f['Output001/mag_DES-z_o_tot'].value + tomag
                rz = r-z ; gr = g-r
    
                # Intrinsic
                ind = np.where(lum>0.)
                if (np.shape(ind)[1] > 0.):
                    ll = np.log10(lum[ind]) + 40.
                    H, xedges, yedges = np.histogram2d(rz[ind],ll,
                                                       bins=[np.append(pbins,pmax),np.append(lbins,lmax)])
                    lo2prop = lo2prop + H
    
                # Extincted
                ind = np.where(lum_ext>0.)
                if (np.shape(ind)[1] > 0.):
                    ll = np.log10(lum_ext[ind]) + 40.                    
                    H, xedges, yedges = np.histogram2d(rz_ext[ind],ll,
                                                       bins=[np.append(pbins,pmax),np.append(lbins,lmax)])
                    lo2prop_ext = lo2prop_ext + H
    
                f.close()
            else:
                print 'NOT found:',gfile
                
    if (volume>0.):
        lo2prop = lo2prop/volume/dl/dp  
        lo2prop_ext = lo2prop_ext/volume/dl/dp  
        print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Figure http://matplotlib.org/users/gridspec.html
    fig = plt.figure(figsize=(8.5,9.))
    ax = fig.add_subplot(1,1,1)

    # Prop vs LOII
    xtit="${\\rm log}_{10}(L\\rm{"+lline+"}/h^{-2}erg\, s^{-1})$"
    ytit="(r-z)"
    xmin=lmin ; xmax=43 ; ymin=pmin ; ymax=2.5
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
    
    # Extincted
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    zp = np.log10(lo2prop_ext)
    xx,yy = np.meshgrid(lbins,pbins)
    al = [-3.5,-2.5,-1.5] 
    cs = ax.contour(xx, yy, zp, levels=al,
                    colors='black',linewidths=3.)
    cs.levels = [nf(val) for val in cs.levels]
    ax.clabel(cs, cs.levels, inline=1,inline_spacing=0,\
              fontsize=10,fmt='%r')

    # Intrinsic
    matplotlib.rcParams['contour.negative_linestyle'] = 'dotted'
    zp = np.log10(lo2prop) 
    xx,yy = np.meshgrid(lbins,pbins)
    al = [-3.5,-2.5,-1.5] 
    cs = ax.contour(xx, yy, zp, levels=al,
                    colors='grey',linewidths=1.5)
    cs.levels = [nf(val) for val in cs.levels]
    ax.clabel(cs, cs.levels, inline=1,inline_spacing=0,\
              fontsize=10,fmt='%r')
 
    zleg = 'z={:.2f}'.format(zz)
    ax.text(xmax-0.2*(xmax-xmin),ymax-0.05*(ymax-ymin),zleg)

    # Save figures
    plotfile = plotdir + 'rz_z'+str(zsnap)+'.pdf'
    fig.savefig(plotfile)
    print 'Output: ',plotfile


#print(min(r),max(r),min(r_ext),max(r_ext)): (18.228985, 35.46558, 19.548662, 40.938416)
#print(min(g),max(g),min(g_ext),max(g_ext)): (18.029503, 36.001396, 19.714731, 41.651962)
#print(min(rz),max(rz),min(rz_ext),max(rz_ext)): (-0.25122643, 2.0978699, -0.19802475, 5.591667)
#print(min(gr),max(gr),min(gr_ext),max(gr_ext)): (-0.6182518, 2.4435616, -0.32359123, 2.4566479)
