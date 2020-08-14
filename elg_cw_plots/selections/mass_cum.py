#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import read_jc_obs as jc
from Cosmology import * 
from stats import n_gt_x
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp19/'

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/'+model+'selections/mass_cum_'
plotfile = outdir+line+'.pdf'
#############################

snap_list = [41,39] #MillGas
nvol = 64

#####
obsnom = ['DEEP2','VVDSDEEP','VVDSWIDE']
obands = ['R24.2','I24','I22.5']

bands = ['DEIMOS-R','MegaCam-i-atmos','MegaCam-i-atmos','eBOSS-SGC','DESI']
mcuts = [24.1, 24, 22.5]
fcuts = [2.7*10.**-17., 1.9*10.**-17., 3.5*10.**-17.,10.**-16.,8.*10.**-17.]

inleg = ['All','DEEP2','VVDS-DEEP','VVDS-Wide','eBOSS-SGC','DESI']
##########
# No VVDSWIDE
obsnom = ['DEEP2','VVDSDEEP']
obands = ['R24.2','I24']

bands = ['DEIMOS-R','MegaCam-i-atmos','eBOSS-SGC','DESI']
mcuts = [24.1, 24]
fcuts = [2.7*10.**-17., 1.9*10.**-17., 10.**-16.,8.*10.**-17.]

inleg = ['All','DEEP2','VVDS-DEEP','eBOSS-SGC','DESI']
##########

ntypes = len(inleg)
zleg = []
cols = get_distinct(ntypes-1)
cols.insert(0,'k')

# Initialize bins
pmin = 8.
pmax = 13.
dp = 0.1
pbins = np.arange(pmin,pmax,dp)

############################################
# Initialize the parameters for the figures
fig = plt.figure(figsize=(6.,12.))
ax = fig.add_subplot(111) 
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_xlabel('${\\rm log}_{10}({\\rm M_*}/M_{\odot}h^{-1})$')
ax.set_ylabel('${\\rm log}_{10}(n_{\\rm gal}(>X)/Mpc^{-3}h^3)$')

xmin = 8.5 ; xmax = 12.
ymin = -5.9 ; ymax = 0.5

# Loop over the redshifts of interest
jj = 210
for iz,zsnap in enumerate(snap_list):
    jj = jj + 1

    ncum = np.zeros(shape=(ntypes,len(pbins)))
    ncum_ext = np.zeros(shape=(ntypes,len(pbins)))

    volume = 0. ; firstpass = True
    for ivol in range(nvol):
        gfile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (os.path.isfile(gfile)):
            # Get some of the model constants
            f = h5py.File(gfile,'r')
            group = f['Parameters']
            vol1 = group['volume'].value ; volume = volume + vol1
            h0 = group['h0'].value 
            omega0 = group['omega0'].value
            omegab = group['omegab'].value
            lambda0 =group['lambda0'].value

            group = f['Output001']
            zz   = group['redshift'].value
            mdisk = group['mstars_disk'].value 
            mbulge = group['mstars_bulge'].value
            mass1 = mdisk + mbulge 

            lmass = np.zeros(len(mass1)) ; lmass.fill(-999.)
            ind = np.where(mass1>0.) 
            lmass[ind] = np.log10(mass1[ind])
            f.close()

            if(firstpass):
                zstring = "{:.2f}".format(zz) 
                zleg.append('z='+zstring)
                firstpass = False
            set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                              h0=h0, universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz) 

            efile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r')
                lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                lum = f['Output001/L_tot_'+line].value

                for index in range(ntypes):
                    if index==0:
                        ind  = np.where(lum_ext>0.)
                        indi = np.where(lum>0.)
                    elif ((inleg[index] == 'eBOSS-SGC') or (inleg[index] == 'DESI')):
                        fluxcut = fcuts[index-1]
                        lcut = emission_line_luminosity(fluxcut,zz)

                        g = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                        r = f['Output001/mag_DES-r_o_tot_ext'].value + tomag
                        z = f['Output001/mag_DES-z_o_tot_ext'].value + tomag
                        rz = r-z ; gr = g-r

                        if (inleg[index] == 'eBOSS-SGC'): 
                            ind = np.where((lum_ext>lcut) & \
                                               (g>21.825) & (g<22.825) & \
                                               (gr>-0.068*rz + 0.457) & \
                                               (gr<0.112*rz + 0.773) & \
                                               (rz>0.218*gr + 0.571) & \
                                               (rz<-0.555*gr + 1.901))
                            indi = np.where((lum>lcut) & \
                                               (g>21.825) & (g<22.825) & \
                                               (gr>-0.068*rz + 0.457) & \
                                               (gr<0.112*rz + 0.773) & \
                                               (rz>0.218*gr + 0.571) & \
                                               (rz<-0.555*gr + 1.901))

                        elif (inleg[index] == 'DESI'): 
                            ind  = np.where((r<23.4) & \
                                                (rz>0.3) & (gr>-0.3) & \
                                                (gr<1.1*rz-0.13) & \
                                                (gr<1.6-1.18*rz) & \
                                                (lum_ext>lcut))

                            indi  = np.where((r<23.4) & \
                                                (rz>0.3) & (gr>-0.3) & \
                                                (gr<1.1*rz-0.13) & \
                                                (gr<1.6-1.18*rz) & \
                                                (lum>lcut))

                    else:
                        ib = bands[index-1]
                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value\
                            + tomag
                        icut = mcuts[index-1]
                        fluxcut = fcuts[index-1]
                        lcut = emission_line_luminosity(fluxcut,zz)

                        ind  = np.where((mag<icut) & (lum_ext>lcut))
                        indi = np.where((mag<icut) & (lum>lcut))
                        
                    if (np.shape(ind)[1] > 0.):
                        ll = lmass[ind]
                        H = n_gt_x(pbins,ll)
                        ncum_ext[index,:] = ncum_ext[index,:] + H

                    if (np.shape(indi)[1] > 0.):
                        ll = lmass[indi]
                        H = n_gt_x(pbins,ll)
                        ncum[index,:] = ncum[index,:] + H
                f.close()

    ncum = ncum/dp/volume 
    ncum_ext = ncum_ext/dp/volume
    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Write output
    outfile = '/cosma5/data/durham/violeta/lines/cosmicweb/selections/'+model+'mass_cum_sn'+\
        str(zsnap)+'.dat'
    with open(outfile,'w') as outf:
        outf.write('# log10(M/Msun/h) log10(ngal for ')
        outf.write(' '.join(str(i) for i in inleg))
        outf.write(') \n')

        tofile1 = np.copy(np.transpose(ncum_ext))
        ind = np.where(tofile1<1e-8)
        tofile1[ind] = -999.
        
        ind = np.where(tofile1>0.)
        tofile1[ind] = np.log10(tofile1[ind])

        tofile = np.column_stack((pbins,tofile1))

        np.savetxt(outf, tofile, fmt='%.5f')
    print('Output: ',outfile)

    # Plot
    if (iz == 0):
        ax1 = fig.add_subplot(jj) ; ax1.set_autoscale_on(False)
        ax1.set_xlim([xmin,xmax]) ; ax1.set_ylim([ymin,ymax]) 
        ax1.set_autoscale_on(False) ;  ax1.minorticks_on()
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.text(xmin+0.1, ymin+0.2, zleg[iz])
    else:
        ax2 = fig.add_subplot(jj,sharex=ax1,sharey=ax1)
        ax2.set_autoscale_on(False) ;  ax2.minorticks_on()
        ax2.text(xmin+0.1, ymin+0.2 , zleg[iz])

    # Plot the model predictions
    for index in range(ntypes):
        # Attenuated
        py = np.zeros(len(pbins)) ; py.fill(-999.)
        a = np.copy(ncum_ext[index,:])
        ind = np.where(a > 0) ; py[ind] = np.log10(a[ind])
        if (iz == 0):
            ax1.step(pbins,py,color=cols[index],linestyle='-',\
                         label=inleg[index])
        else:
            ax2.step(pbins,py,color=cols[index],linestyle='-',\
                        label=inleg[index])

        ## Intrinsic
        #py = np.zeros(len(pbins)) ; py.fill(-999.)
        #a = np.copy(ncum[index,:])
        #ind = np.where(a > 0) ; py[ind] = np.log10(a[ind])
        #if (iz == 0):
        #    ax1.step(pbins,py,color=cols[index],linestyle='-',\
        #                 label=inleg[index])
        #else:
        #    ax2.step(pbins,py,color=cols[index],linestyle='-',\
        #                label=inleg[index])

        # Legend
        if (iz == len(snap_list)-1):
            leg = plt.legend(loc=1)
            renderer = fig.canvas.get_renderer()
            shift1 = max([t.get_window_extent(renderer).width for t in leg.get_texts()])
            shift2 = min([t.get_window_extent(renderer).width for t in leg.get_texts()])
            for item in leg.legendHandles:
                item.set_visible(False)
            for color,text in zip(cols,leg.get_texts()):
                text.set_color(color)
                text.set_ha('right') ; text.set_position((shift1-shift2,0))
                leg.draw_frame(False)           

# Plot number densities
x = [xmin-xmin*0.5,xmax+xmax*0.5]
for nd in [-2.,-3.,-4.2]:
    y = [nd,nd]
    ndtext = '$10^{'+str(nd)+'}$'

    ax1.plot(x,y,color='grey',linestyle=':')
    ax1.text(xmin+0.1, nd-0.3, ndtext, color='grey',fontsize=11)
    
    ax2.plot(x,y,color='grey',linestyle=':')
    ax2.text(xmin+0.1, nd-0.3, ndtext, color='grey',fontsize=11)

# Save figures
fig.subplots_adjust(hspace=0)
fig.savefig(plotfile)
print 'Output: ',plotfile
