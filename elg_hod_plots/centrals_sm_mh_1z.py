#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from distinct_colours import get_distinct
from Cosmology import *
import stats as s 
import mpl_style
plt.style.use(mpl_style.style1)

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
model = 'MillGas/gp15newmg.anders/' #model = 'MillGas/gp14/'

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/desi_hod_o2/'
#############################

snap_list =[44, 42, 40, 37, 34]

nvol = 64

zleg = []
bands = ['RK','m2','m2']
mcuts = [24.1,22.5,24]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17.,10.**-16.,8.*10.**-17.]
inleg = ['DEEP2 cuts','VVDS-Wide cuts','VVDS-Deep cuts','eBOSS cuts','DESI cuts']

cols = ['gray'] ; cols.extend(get_distinct(len(inleg)))
# Initialize histogram
lmin = 8.5
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
ytit="$log_{10} (\\rm M_{*}/M_{\odot}h^{-1})$"
xtit="$log_{10} (\\rm M_{\\rm halo}/M_{\odot}h^{-1})$"

ymin = 8. ; ymax = 11.5
xmin = 10.5 ; xmax = 13.5

# Plot
jj = 111 ; first = True
for iz,zsnap in enumerate(snap_list):
    nall_med, nall_p1, nall_p9 = [np.zeros(shape=(len(lhist))) for i in range(3)]
    nelgs_med, nelgs_p1, nelgs_p9 = [np.zeros(shape=(len(bands),len(lhist))) for i in range(3)]
    nboss_med, nboss_p1, nboss_p9 = [np.zeros(shape=(len(lhist))) for i in range(3)]
    ndesi_med, ndesi_p1, ndesi_p9 = [np.zeros(shape=(len(lhist))) for i in range(3)]

    volume = 0.  ; firstpass = True
    for ivol in range(nvol):
        gfile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
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

            if(firstpass):
                zstring = "{:.2f}".format(zz) 
                zleg.append('z='+zstring)
                plotfile = outdir+'plots/cuts/'+model+\
                    '/sm_mh/centrals_sm_mh_'+\
                    line+'_z'+zstring+'.pdf'
                firstpass = False
            set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                              h0=h0, universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz)

            efile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r')
                # Number of haloes per bin
                if (first):
                    mhhalo = f['Output001/mchalo'].value
                    mass = f['Output001/mstars_tot'].value
                    types = f['Output001/type'].value
                    for index,ib in enumerate(bands):
                        if (index==0):
                            mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                            lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                        else:
                            arr = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                            mag = np.column_stack((mag,arr))

                            arr = f['Output001/L_tot_'+line+'_ext'].value
                            lum_ext = np.column_stack((lum_ext,arr))
                            
                    #eBOSS & DESI
                    g = f['Output001/magDgo_tot_ext'].value + tomag
                    r = f['Output001/magDro_tot_ext'].value + tomag
                    z = f['Output001/magDzo_tot_ext'].value + tomag
                    rz = r-z ; gr = g-r

                    #eBOSS decam180 selection
                    fluxcut = 10.**-16.
                    lcut = emission_line_luminosity(fluxcut,zz)
                    ind = np.where((g>22.1) & (g<22.8) & \
                                       (gr>0.3) & (gr<0.7) & \
                                       (rz>0.25) & (rz<1.4) & \
                                       (rz>0.5*gr+0.4) & \
                                       (rz<0.5*gr+0.8) & \
                                       (f['Output001/type'].value==0) & \
                                       (f['Output001/L_tot_'+line+'_ext'].value>lcut))
                    mh_eboss = mhhalo[ind]
                    mass_eboss = mass[ind]
                    #DESI
                    fluxcut = 8.*10.**-17.
                    lcut = emission_line_luminosity(fluxcut,zz)
                    ind = np.where((r<23.4) & \
                                       (rz>0.3) & (gr>-0.3) & \
                                       (rz>0.9*gr+0.12) & \
                                       (rz<1.345-0.85*gr) & \
                                       (f['Output001/type'].value==0) & \
                                       (f['Output001/L_tot_'+line+'_ext'].value>lcut))
                    mh_desi = mhhalo[ind]
                    mass_desi = mass[ind]

                    first = False
                else:
                    in_mh = f['Output001/mchalo'].value ; mhhalo = np.append(mhhalo,in_mh)
                    in_mass= f['Output001/mstars_tot'].value ; mass = np.append(mass,in_mass)
                    types = np.append(types,f['Output001/type'].value)
                    for index,ib in enumerate(bands):
                        if (index==0):
                            inmag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                            inlum_ext = f['Output001/L_tot_'+line+'_ext'].value
                        else:
                            arr = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                            inmag = np.column_stack((inmag,arr))

                            arr = f['Output001/L_tot_'+line+'_ext'].value
                            inlum_ext = np.column_stack((inlum_ext,arr))
                    mag = np.vstack((mag,inmag))
                    lum_ext = np.vstack((lum_ext,inlum_ext))

                    #eBOSS & DESI
                    g = f['Output001/magDgo_tot_ext'].value + tomag
                    r = f['Output001/magDro_tot_ext'].value + tomag
                    z = f['Output001/magDzo_tot_ext'].value + tomag
                    rz = r-z ; gr = g-r

                    #eBOSS decam180 selection
                    fluxcut = 10.**-16.
                    lcut = emission_line_luminosity(fluxcut,zz)
                    ind = np.where((g>22.1) & (g<22.8) & \
                                       (gr>0.3) & (gr<0.7) & \
                                       (rz>0.25) & (rz<1.4) & \
                                       (rz>0.5*gr+0.4) & \
                                       (rz<0.5*gr+0.8) & \
                                       (f['Output001/type'].value==0) & \
                                       (f['Output001/L_tot_'+line+'_ext'].value>lcut))
                    mh_eboss = np.append(mh_eboss,in_mh[ind])
                    mass_eboss = np.append(mass_eboss,in_mass[ind])
                    #DESI
                    fluxcut = 8.*10.**-17.
                    lcut = emission_line_luminosity(fluxcut,zz)
                    ind = np.where((r<23.4) & \
                                       (rz>0.3) & (gr>-0.3) & \
                                       (rz>0.9*gr+0.12) & \
                                       (rz<1.345-0.85*gr) & \
                                       (f['Output001/type'].value==0) & \
                                       (f['Output001/L_tot_'+line+'_ext'].value>lcut))
                    mh_desi = np.append(mh_desi,in_mh[ind])
                    mass_desi = np.append(mass_desi,in_mass[ind])

                f.close()

    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # All 
    ind = np.where((mass>0) & (mhhalo>0) & (types==0))
    y = np.log10(mass[ind]) ; x = np.log10(mhhalo[ind])
    w = np.zeros(shape=(len(x))) ; w.fill(1./volume/dl)
    nmin=100
    nall_med = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.5)
    nall_p1 = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.1)
    nall_p9 = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.9)
    
    # Galaxies selected using the survey cuts
    for index,ib in enumerate(bands):
        icut = mcuts[index] ; fluxcut = fcuts[index]
        lcut = emission_line_luminosity(fluxcut,zz)

        # ELGs
        inmag = mag[:,index] ; inlum_ext = lum_ext[:,index]
        ind = np.where((inmag<icut) & (inlum_ext>lcut)  & (types==0))
        if (np.shape(ind)[1] > 0):
            y = np.log10(mass[ind])
            x = np.log10(mhhalo[ind])
            w = np.zeros(shape=(len(x))) ; w.fill(1./volume/dl)
            nelgs_med[index,:] = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.5)
            nelgs_p1[index,:] = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.1)
            nelgs_p9[index,:] = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.9)
            print inleg[index]
            print nelgs_med,nelgs_p1,nelgs_p9

    # eBOSS 
    ind = np.where((mass_eboss>0) & (mh_eboss>0))
    y = np.log10(mass_eboss[ind]) ; x = np.log10(mh_eboss[ind])
    w = np.zeros(shape=(len(x))) ; w.fill(1./volume/dl)
    nmin=100
    neboss_med = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.5)
    neboss_p1 = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.1)
    neboss_p9 = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.9)
    print 'eBOSS'
    print neboss_med,neboss_p1,neboss_p9

    # DESI
    ind = np.where((mass_desi>0) & (mh_desi>0))
    y = np.log10(mass_desi[ind]) ; x = np.log10(mh_desi[ind])
    w = np.zeros(shape=(len(x))) ; w.fill(1./volume/dl)
    nmin=100
    ndesi_med = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.5)
    ndesi_p1 = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.1)
    ndesi_p9 = s.perc_2arrays(np.append(lbins,lmax),x,y,w,nmin,0.9)
    print 'DESI'
    print ndesi_med,ndesi_p1,ndesi_p9

    # Plot
    fig = plt.figure(figsize=(8.5,9.))

    ax = fig.add_subplot(jj)
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
    #ax.tick_params(labelsize=fs-2)   
    #start, end = ax.get_xlim()
    #ax.xaxis.set_ticks(np.arange(start, end, 1.))
    #ax.xaxis.set_ticks(np.arange(start, end, 0.2),minor=True)
    #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    ax.text(xmax-(xmax-xmin)*0.2, ymin+(ymax-ymin)*0.05, zleg[iz])

    # Plot the model predictions
    py = nall_med ; ind = np.where(py>-999.)
    ax.plot(lhist[ind],py[ind],color=cols[0],linestyle='-',\
                label='All centrals')
    ax.plot(lhist[ind],nall_p1[ind],color=cols[0],linestyle='--')
    ax.plot(lhist[ind],nall_p9[ind],color=cols[0],linestyle='--')

    for index,ib in enumerate(bands):        
        py = nelgs_med[index,:] ; ind = np.where(py>-999.)
        ax.plot(lhist[ind],py[ind],color=cols[index+1],linestyle='-',label=inleg[index])
        if(index==0):
            y=nelgs_p1[index,:] ; ax.plot(lhist[ind],y[ind],color=cols[index+1],linestyle='--')
            y=nelgs_p9[index,:] ; ax.plot(lhist[ind],y[ind],color=cols[index+1],linestyle='--')
                
    # Plot eBOSS
    ii = 3
    py = neboss_med ; ind = np.where(py>-999.)
    ax.plot(lhist[ind],py[ind],color=cols[ii+1],linestyle='-',label=inleg[ii])
    #ax.plot(lhist[ind],neboss_p1[ind],color=cols[ii],linestyle='--')
    #ax.plot(lhist[ind],neboss_p9[ind],color=cols[ii],linestyle='--')

    # Plot DESI
    ii = 4
    py = ndesi_med ; ind = np.where(py>-999.)
    ax.plot(lhist[ind],py[ind],color=cols[ii+1],linestyle='-',label=inleg[ii])
    #ax.plot(lhist[ind],ndesi_p1[ind],color=cols[ii],linestyle='--')
    #ax.plot(lhist[ind],ndesi_p9[ind],color=cols[ii],linestyle='--')


    # Legend
    leg = plt.legend(loc=2, handlelength=0, handletextpad=0)
    for item in leg.legendHandles:
        item.set_visible(False)
    for color,text in zip(cols,leg.get_texts()):
        text.set_color(color)
        leg.draw_frame(False)

    # Save figures
    #fig.tight_layout()
    fig.savefig(plotfile)
    print 'Output: ',plotfile
