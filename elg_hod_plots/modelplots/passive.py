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

path6 = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/'
path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/'
nvol = 64

#plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/'
#models = ['gp14','gp15newmg']
#inleg = ['GP14','This work']

plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/anders.'
models = ['gp15newmg','gp15newmg.anders'] ; inleg = models

# Initialize GSMF
mmin = 8.5
mmax = 15.
dm = 0.1
mbins = np.arange(mmin,mmax,dm)
mhist = mbins + dm*0.5

ntot  = np.zeros(shape=(len(models),len(mbins)))
npas1 = np.zeros(shape=(len(models),len(mbins)))
npas2 = np.zeros(shape=(len(models),len(mbins)))

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


# Loop over redshifts
for index,model in enumerate(models):
    print model
    volume = 0. ; first = True
    for ivol in range(nvol):
        gfile = path+model+'/iz61/ivol'+str(ivol)+'/galaxies.hdf5'
        if(os.path.isfile(gfile)):
            # Read the relevant information from galaxies.hdf5
            f = h5py.File(gfile,'r') #; print gfile
            vol1 = f['Parameters/volume'].value ; volume = volume + vol1
            h0 =   f['Parameters/h0'].value 
            omega0 = f['Parameters/omega0'].value
            omegab = f['Parameters/omegab'].value
            lambda0 =f['Parameters/lambda0'].value


            group = f['Output001']
            zz   = group['redshift'].value
            set_cosmology(omega0=omega0,omegab=omegab, \
                              lambda0=lambda0,h0=h0, \
                              universe="Flat",include_radiation=False)
            slim = 1./tHubble(zz) 
            if first:
                print 'zz, slim, 0.3*slim=',zz, slim, 0.3*slim
                first = False

            mdisk = group['mstars_disk'].value
            mbulge = group['mstars_bulge'].value
            mass1 = mdisk + mbulge
            sdisk = group['mstardot'].value # Msolar/h/Gyr
            sbulge = group['mstardot_burst'].value
            sfr1 = sdisk + sbulge            

            ind = np.where((mass1>0.) & (sfr1>0.))
            ssfr = np.zeros(shape=(len(sfr1)))
            ssfr[ind] = sfr1[ind]/mass1[ind]
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            ntot[index,:] = ntot[index,:] +H

            ind = np.where((mass1>0.) & (ssfr<0.3*slim))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            npas1[index,:] = npas1[index,:] +H

            ind = np.where((mass1>0.) & (ssfr<slim))
            mass = np.log10(mass1[ind])
            H, bins_edges =\
                np.histogram(mass,bins=np.append(mbins,mmax))
            npas2[index,:] = npas2[index,:] +H

            f.close()
        else:
            print 'NOT found:',gfile
            
    if (volume>0.):
        print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)
        for i in range(len(mbins)):
            if ntot[index,i]>0.:
                npas1[index,i] = npas1[index,i]/ntot[index,i]
                npas2[index,i] = npas2[index,i]/ntot[index,i]
            else:
                npas1[index,i] = -999.
                npas2[index,i] = -999.


# Figure http://matplotlib.org/users/gridspec.html
fig = plt.figure(figsize=(8.5,9.))
ax = plt.subplot()
cols = get_distinct(len(models)) 
colors = cols ; g = ['grey'] ; colors.extend(g) ; colors.extend(g)
xtit = "$log(\\rm{M_{*}/M_{\odot} h^{-1}})$"
ytit="Passive fraction"
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
ax.set_autoscale_on(False) ;  ax.minorticks_on()
ax.set_xlim(mmin,12) ; ax.set_ylim(0.,1.) 

# Observations
col = colors[len(colors)-1]

fobs = 'Obs_Data/gilbank10.txt'
#log(M_*)   passive_fraction error
lm, p, erro = np.loadtxt(fobs, unpack=True)
xo = lm + np.log10(0.7)
ax.errorbar(xo,p, yerr=erro, color=col, ecolor=col,\
                label ='Gilbank+10', fmt = 'o')

fobs = 'Obs_Data/bauer13.txt'
#log(M_*)   passive_fraction error
lm, p = np.loadtxt(fobs, unpack=True)
xo = lm + np.log10(0.7)
erro = lm*0.
ax.errorbar(xo,p, yerr=erro, color=col, ecolor=col,\
                label ='Bauer+10', fmt = '^')

# Models
for ii in range(len(models)):
    py = npas1[ii,:] ; ind = np.where(py>0.)
    x = mhist[ind] ; y = py[ind]
    ax.plot(x[ind],y[ind],color=colors[ii],label=inleg[ii])

    #py = npas2[ii,:] ; ind = np.where(py>0.)
    #y = py[ind]
    #ax.plot(x[ind],y[ind],color=colors[ii],\
    #            linestyle='--')

## Legend
leg = ax.legend(loc=2,fontsize='small')
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)


# Save figures
plotfile = plotdir + 'passive.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile

