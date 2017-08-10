#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import read_jc_obs as jc
from Cosmology import * 
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)
cols = get_distinct(5)

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
model = 'MillGas/gp15newmg.anders/' #model = 'MillGas/gp14/'

zsnap = '44' ; inz = '0.62' 
#zsnap = '42' ; inz = '0.76' 
#zsnap = '40' ; inz = '0.91' 
#zsnap = '37' ; inz = '1.17' 
#zsnap = '34' ; inz = '1.50' 
nvol = 64

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/cuts/'+model+'dust_deep2vvds_z'
############################# Obs
obsh0 = 0.677
obs_dir = '/gpfs/data/violeta/lines/desi_hod_o2/lf_obs_data/lf_may16_comparat/individual_LF/'
#############################

plotfile = outdir+inz+'_'+line+'.pdf'
zleg = 'z='+inz

obsnom = ['DEEP2','VVDSWIDE','VVDSDEEP']
obands = ['R24.2','I22.5','I24']

bands = ['RK','m2','m2']
mcuts = [24.1,22.5,24] 
fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17.]

inleg = ['DEEP2 cuts','VVDS-Wide cuts','VVDS-Deep cuts']


# Initialize histogram
lmin = 38.
lmax = 46.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
fig = plt.figure(figsize=(8.,9.))

xtit = "${\\rm log}_{10}(L\\rm{"+lline+"}/h^{-2}erg\, s^{-1})$"
ytit = "${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dlog}_{10}L)$"

xmin = 40. ; xmax = 44.
ymin = -5.5 ; ymax = -1.

lf_ext = np.zeros(shape=(len(bands),len(lhist)))
lf = np.zeros(shape=(len(bands),len(lhist)))

volume = 0.
for ivol in range(nvol):
    gfile = path+model+'/iz'+zsnap+'/ivol'+str(ivol)+'/galaxies.hdf5'
    if (os.path.isfile(gfile)):
        # Get some of the model constants
        f = h5py.File(gfile,'r')
        zz   = f['Output001/redshift'].value
        group = f['Parameters']
        vol1 = group['volume'].value ; volume = volume + vol1
        h0 = group['h0'].value 
        omega0 = f['Parameters/omega0'].value
        omegab = f['Parameters/omegab'].value
        lambda0 = f['Parameters/lambda0'].value
        tsf    = 3*f['Output001/tsfburst'].value 
        tburst = f['Output001/tburst'].value
        f.close()

        set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                          h0=h0, universe="Flat",include_radiation=False)
        tomag = band_corrected_distance_modulus(zz) 

        efile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
        if (os.path.isfile(efile)):
            f = h5py.File(efile,'r')
            lum_ext = f['Output001/L_tot_'+line+'_ext'].value
            lum = f['Output001/L_tot_'+line].value
            BoT = f['Output001/BoT'].value

            if(len(lum)==len(tsf)):
                for index,ib in enumerate(bands):
                    mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                    icut = mcuts[index] ;  fluxcut = fcuts[index]
                    lcut = emission_line_luminosity(fluxcut,zz)

                    # All w dust
                    ind  = np.where((mag<icut) & (lum_ext>lcut))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H

                    # All no dust
                    ind  = np.where((mag<icut) & (lum>lcut))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H


            else:
                print 'WARNING: arrays of different size for ivol=',ivol

            f.close()

lf_ext = lf_ext/dl/volume
lf = lf/dl/volume
print zz,', Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

# Plot
jj = 111
ax = fig.add_subplot(jj) ; ax.set_autoscale_on(False)
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
start, end = ax.get_xlim()
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.text(40.1, -1.3, zleg)

# Plot the observations  O2_3728-VVDSWIDEI21.5-z1.234.txt
for i,isurvey in enumerate(obsnom):        
    ox, oy, el, eh = jc.read_jc_indlf(obs_dir,zz,h0=obsh0,line='O2_3728',\
                                          survey=isurvey,band=obands[i])
    if(isinstance(ox, (np.ndarray))):
        ax.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                    ecolor=cols[i],color=cols[i],mec=cols[i])

# Plot the model predictions
for index,ib in enumerate(bands):
    py = 0. ; py = lf_ext[index,:]
    ind = np.where(py > 0)
    x = lhist[ind]
    y = np.log10(py[ind])
    ind = np.where(y < 0.)
    ax.plot(x[ind],y[ind],color=cols[index],linestyle='-',\
                label=inleg[index])

    py = 0. ; py = lf[index,:]
    ind = np.where(py > 0)
    x = lhist[ind]
    y = np.log10(py[ind])
    ind = np.where(y < 0.)
    ax.plot(x[ind],y[ind],color=cols[index],linestyle='--')

                
# Legend
leg = plt.legend(loc=1, handlelength=0, handletextpad=0)
for item in leg.legendHandles:
    item.set_visible(False)
for color,text in zip(cols,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)

# Save figures
fig.savefig(plotfile)
print 'Output: ',plotfile
