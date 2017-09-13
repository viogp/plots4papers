#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from Cosmology import * 
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
model = 'MillGas/gp15newmg/' #model = 'MillGas/gp14/'

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/desi_hod_o2/'
plotfile = outdir+'plots/cuts/'+model+'1z_hod_'+line+'_z0.75_bot.pdf'
#############################

snap_list = [42] #MillGas
zleg = ['z=0.75']

nvol = 64

bands = ['RK','m2','m2']
mcuts = [24.1,22.5,24]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17.,10.**-17., 10.**-17.]

surveys = ['DEEP2','VVDS-WIDE','VVDS-DEEP','eBOSS','DESI']   

# Initialize histogram
lmin = 8.5
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
cols = get_distinct(len(surveys))

xtit = "${\\rm log}_{10}({\\rm M_{halo}}/M_{\odot}h^{-1})$"
ytit = "$\\langle N\\rangle _{\\rm [OII] emitters}$"

xmin = 10. ; xmax = 15.
ymin = -4. ; ymax = 1.

# Loop over the redshifts of interest
for iz,zsnap in enumerate(snap_list):
    nm = np.zeros(shape=(len(surveys),len(lhist)))
    nc = np.zeros(shape=(len(surveys),len(lhist)))
    cs = np.zeros(shape=(len(surveys),len(lhist)))
    cd = np.zeros(shape=(len(surveys),len(lhist)))
    nh = np.zeros(shape=(len(lhist)))

    volume = 0.
    for ivol in range(nvol):
        gfile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (os.path.isfile(gfile)):
            # Get some of the model constants
            f = h5py.File(gfile,'r') #; print gfile
            zz   = f['Output001/redshift'].value
            vol1 = f['Parameters/volume'].value ; volume = volume + vol1
            h0 =   f['Parameters/h0'].value 
            omega0 = f['Parameters/omega0'].value
            omegab = f['Parameters/omegab'].value
            lambda0 =f['Parameters/lambda0'].value

            set_cosmology(omega0=omega0,omegab=omegab,\
                              lambda0=lambda0, h0=h0,\
                              universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz)

            efile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r') #; print efile
                # Number of haloes per bin
                gtype  = f['Output001/type'].value
                mhhalo = f['Output001/mhhalo'].value
                lum_ext =f['Output001/L_tot_'+line+'_ext'].value
                BoT = f['Output001/BoT'].value
                ind = np.where(gtype == 0)
                ll = np.log10(mhhalo[ind])
                H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                nh = nh + H
            
                for index,ib in enumerate(bands):
                    mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                    icut = mcuts[index] ; fluxcut = fcuts[index]
                    lcut = emission_line_luminosity(fluxcut,zz)

                    # All
                    ind = np.where((mag<icut) & (lum_ext>lcut))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        nm[index,:] = nm[index,:] + H

                    # Centrals
                    ind = np.where((mag<icut) & \
                                       (lum_ext>lcut) & (gtype==0))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        nc[index,:] = nc[index,:] + H

                    # Centrals Spheroids
                    ind = np.where((mag<icut) & (lum_ext>lcut) &\
                                       (gtype==0) & (BoT>0.5))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cs[index,:] = cs[index,:] + H

                    # Centrals Disk
                    ind = np.where((mag<icut) & (lum_ext>lcut) &\
                                       (gtype==0) & (BoT<=0.5))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cd[index,:] = cd[index,:] + H

                #----------------------
                #eBOSS & DESI
                g = f['Output001/magDgo_tot_ext'].value + tomag
                r = f['Output001/magDro_tot_ext'].value + tomag
                z = f['Output001/magDzo_tot_ext'].value + tomag
                rz = r-z ; gr = g-r

                #eBOSS decam180 selection ----------------------
                index = len(surveys)-2  #; print inleg[index]
                lcut = emission_line_luminosity(fcuts[index],zz)

                # All
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & \
                                   (lum_ext>lcut))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    nm[index,:] = nm[index,:] + H

                # Centrals
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & \
                                   (lum_ext>lcut) & (gtype==0))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    nc[index,:] = nc[index,:] + H

                # Censtral spheroids
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & \
                                   (lum_ext>lcut) & \
                                   (gtype==0) & (BoT>0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cs[index,:] = cs[index,:] + H

                # Central disks
                ind = np.where((g>22.1) & (g<22.8) & \
                                   (gr>0.3) & (gr<0.7) & \
                                   (rz>0.25) & (rz<1.4) & \
                                   (rz>0.5*gr+0.4) & \
                                   (rz<0.5*gr+0.8) & \
                                   (lum_ext>lcut) & \
                                   (gtype==0) & (BoT<=0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cd[index,:] = cd[index,:] + H


                #DESI -------------------------------------------------
                index = len(surveys)-1  #; print inleg[index],fcuts[index]
                lcut = emission_line_luminosity(fcuts[index],zz)

                # All
                ind = np.where((g>20.45) & (r<22.79) & \
                                   (rz>0.285) & (rz<1.585) & \
                                   (gr<1.1458*rz-0.209) & \
                                   (gr<1.4551-1.1458*rz) & \
                                   (lum_ext>lcut))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    nm[index,:] = nm[index,:] + H

                # Centrals
                ind = np.where((g>20.45) & (r<22.79) & \
                                   (rz>0.285) & (rz<1.585) & \
                                   (gr<1.1458*rz-0.209) & \
                                   (gr<1.4551-1.1458*rz) & \
                                   (lum_ext>lcut) & (gtype==0))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    nc[index,:] = nc[index,:] + H

                # Central spheroids
                ind = np.where((g>20.45) & (r<22.79) & \
                                   (rz>0.285) & (rz<1.585) & \
                                   (gr<1.1458*rz-0.209) & \
                                   (gr<1.4551-1.1458*rz) & \
                                   (lum_ext>lcut) &\
                                   (gtype==0) & (BoT>0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cs[index,:] = cs[index,:] + H

                # Central disks
                ind = np.where((g>20.45) & (r<22.79) & \
                                   (rz>0.285) & (rz<1.585) & \
                                   (gr<1.1458*rz-0.209) & \
                                   (gr<1.4551-1.1458*rz) & \
                                   (lum_ext>lcut) & \
                                   (gtype==0) & (BoT<=0.5))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    cd[index,:] = cd[index,:] + H


                f.close()

    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Plot
    plt.xlim(xmin,xmax) ; plt.ylim(ymin,ymax) 
    plt.xlabel(xtit) ; plt.ylabel(ytit)
    axes = plt.gca() 
    start, end = axes.get_xlim()
    axes.xaxis.set_ticks(np.arange(start, end, 1.))
    axes.xaxis.set_ticks(np.arange(start, end, 0.2),minor=True)
    axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    axes.text(xmin+(xmax-xmin)*0.05, ymax-(ymax-ymin)*0.05, zleg[iz])

    # Plot the model predictions
    for index,survey in enumerate(surveys):
        ileg = survey+' cuts'
        indh = np.where(nh > 0) ; nhalos = sum(nh)
        if (np.shape(indh)[1]<1): continue
        print index,ileg,np.shape(indh)[1]

        py = 0. ; py = nm[index,:] ; nall = sum(py)
        x = lhist[indh] ; y = np.zeros(len(x))
        yall = py[indh]/nh[indh] ; ind = np.where(yall>0.) 
        y[ind] = np.log10(yall[ind]) ; ind = np.where(y < 0.)
        #plt.plot(x[ind],y[ind],color=cols[index],linestyle='-',label=inleg[index])
        #plt.plot(x[ind],y[ind],color=cols[index],linestyle='-')

        #py = 0. ; py = nc[index,:] ; ncen = sum(py)
        #x = lhist[indh]  ; y = np.zeros(len(x))
        #ycen = py[indh]/nh[indh] ; ind = np.where(ycen>0.) 
        #y[ind] = np.log10(ycen[ind]) ;ind = np.where(y < 0.)
        #if (index==0):
        #    plt.plot(x[ind],y[ind],color=cols[index],\
        #                 linestyle='-', label='Centrals')
        #else:
        #    plt.plot(x[ind],y[ind],color=cols[index],linestyle='-')
        
        py = 0. ; py = cs[index,:]  ; nsat = sum(py)
        x = lhist[indh] ; y = np.zeros(len(x))
        ysat = py[indh]/nh[indh] ; ind = np.where(ysat>0.) 
        y[ind] = np.log10(ysat[ind]) ; ind = np.where(y < 0.)
        if (index==0):
            plt.plot(x[ind],y[ind],color=cols[index],linestyle='-',label='Central spheroids')
        else:
            plt.plot(x[ind],y[ind],color=cols[index],linestyle='-')

        py = 0. ; py = cd[index,:]  ; nsat2 = sum(py)
        x = lhist[indh] ; y = np.zeros(len(x))
        y2 = py[indh]/nh[indh] ; ind = np.where(y2>0.) 
        y[ind] = np.log10(y2[ind]) ; ind = np.where(y < 0.)
        if (index==0):
            plt.plot(x[ind],y[ind],color=cols[index],linestyle='--',label='Central disks')
        else:
            plt.plot(x[ind],y[ind],color=cols[index],linestyle='--')
                
        # Legend
        leg = plt.legend(loc=1)
        leg.draw_frame(False)
    # Save figures
    plt.savefig(plotfile)
    print 'Output: ',plotfile
