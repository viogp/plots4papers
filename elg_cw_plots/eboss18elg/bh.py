#! /usr/bin/env python

import numpy as np
import os.path,sys
import h5py
from matplotlib import pyplot as plt
import stats as s
from distinct_colours import get_distinct 
import mpl_style
plt.style.use(mpl_style.style1)

path6 = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/'
plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/'

nvol = 64
# 
models = ['gp14','gp15newmg']
inleg = ['GP14','This work']

# Limits
lmin = 9.
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

med, p1, p9 = [np.zeros(shape=(len(models),len(lhist))) for i in range(3)]

for index,model in enumerate(models):
    volh = 0. ; first = True
    for ivol in range(nvol):
        infile = path6+model+'/iz61/ivol'+str(ivol)+'/galaxies.hdf5'

        if (os.path.isfile(infile)):
            f = h5py.File(infile,'r')
            volh = volh + f['Parameters/volume'].value

            group = f['Output001']
            if (first):
                h0 = f['Parameters/h0'].value
                mstars_bulge = group['mstars_bulge'].value
                mbh = group['M_SMBH'].value
                first = False
            else:
                mstars_bulge = np.append(mstars_bulge,group['mstars_bulge'].value)
                mbh = np.append(mbh,group['M_SMBH'].value)

            f.close()
        else:
            print 'Not found: ',infile

    # Take logs and normalize
    print 'Total volume considered =(', np.power(volh,1./3.), ' Mpc/h)^3'

    ind = np.where((mstars_bulge>0) & (mbh>0))
    y = np.log10(mbh[ind]/h0) ; x = np.log10(mstars_bulge[ind]/h0)
    w = np.zeros(shape=(len(x))) ; w.fill(1./volh/dl)
    nmin=100
    med[index,:] = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.5)
    p1[index,:] = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.1)
    p9[index,:] = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.9)

######################
# Make plots
xtit = "$log(\\rm{M_{*}/M_{\odot} h^{-1}})$"
ytit="$log(\Phi/dlog{\\rm M_{*}}/{\\rm Mpc}^{-3} )$"

cols = get_distinct(len(models)) 
colors = cols ; g = ['grey'] ; colors.extend(g)

fig = plt.figure(figsize=(8.5,9.))
plt.xlim(8.,13.)
plt.ylim(5.,10.) 
plt.xlabel(xtit) ; plt.ylabel(ytit)

# Plot observations from Baldry+2012
dobs = '/cosma/home/violeta/Galform2/galform-2.6.0/Obs_Data2/'
file1 = np.genfromtxt(dobs+'MBH_Mbulge_HR04.data')
#MBH(M_sun) +err -err sigma(km/s) Mbulge(M_sun)
mbh_obs = file1[:,0]
log_mbh_obs = np.log10(mbh_obs)

mbh_obs_upper = file1[:,1]
log_mbh_obs_upper = np.log10(1 + mbh_obs_upper/mbh_obs)

mbh_obs_lower = file1[:,2]
log_mbh_obs_lower = -1*np.log10(1 - mbh_obs_lower/mbh_obs)

mbulge_obs = file1[:,4]
log_mbulge_obs = np.log10(mbulge_obs)

oh = 0.7 
plt.errorbar(log_mbulge_obs,log_mbh_obs, xerr=0.18,\
                 yerr=[log_mbh_obs_upper,log_mbh_obs_lower], fmt='o',\
                 ecolor='grey',color='grey', mec='grey', \
                 label='Haering & Rix 2004')

# Plot predictions
for im in range(len(models)):
    py = med[im,:] ; pl = p1[im,:] ; ph = p9[im,:] 
    ind = np.where(py>-999.)
    plt.fill_between(lhist[ind], pl[ind], ph[ind], color=cols[im], alpha=0.2)

for im in range(len(models)):
    py = med[im,:] 
    ind = np.where(py>-999.)
    plt.plot(lhist[ind],py[ind],cols[im],label=inleg[im])

# Legend
leg = plt.legend(loc=2)
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)

# Save figures
plotfile = plotdir + 'bh.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile
