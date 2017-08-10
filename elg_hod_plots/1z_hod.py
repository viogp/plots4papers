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
#############################

snap_list = [42]#[44, 42, 40, 37, 34] #MillGas    

nvol = 64

#bands = ['RK','m2','m2']
#mcuts = [24.1,22.5,24]
#fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17.,10.**-17., 8*10.**-17.]
#
#inleg = ['DEEP2 cuts','VVDS-WIDE cuts','VVDS-DEEP cuts','eBOSS cuts','DESI cuts']

bands = ['RK','m2']
mcuts = [24.1,22.5]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 10.**-16., 8*10.**-17.]

inleg = ['DEEP2 cuts','VVDS-Wide cuts','eBOSS cuts','DESI cuts']

# Initialize histogram
lmin = 8.5
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

zleg = []
############################################
# Initialize the parameters for the figures
ntypes = len(inleg) 
cols = get_distinct(ntypes)

xtit = "${\\rm log}_{10}({\\rm M_{halo}}/M_{\odot}h^{-1})$"
ytit = "${\\rm log}_{10}(\\langle N\\rangle _{\\rm [OII]})$"

xmin = 10. ; xmax = 15.
ymin = -3. ; ymax = 1.

# Loop over the redshifts of interest
for iz,zsnap in enumerate(snap_list):
    nm = np.zeros(shape=(ntypes,len(lhist)))
    nc = np.zeros(shape=(ntypes,len(lhist)))
    ns = np.zeros(shape=(ntypes,len(lhist)))
    n2 = np.zeros(shape=(ntypes,len(lhist)))
    nh = np.zeros(shape=(len(lhist)))

    volume = 0. ; firstpass = True 
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
            if(firstpass):
                szz = "{:.2f}".format(zz) 
                firstpass = False

            efile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r') #; print efile
                # Number of haloes per bin
                gtype  = f['Output001/type'].value
                mhhalo = f['Output001/mhhalo'].value
                lum_ext =f['Output001/L_tot_'+line+'_ext'].value
                ind = np.where(gtype == 0)
                ll = np.log10(mhhalo[ind])
                H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                nh = nh + H
            
                for index,fc in enumerate(fcuts):
                    fluxcut = fc
                    lcut = emission_line_luminosity(fluxcut,zz)

                    g = f['Output001/magDgo_tot_ext'].value + tomag
                    r = f['Output001/magDro_tot_ext'].value + tomag
                    z = f['Output001/magDzo_tot_ext'].value + tomag
                    rz = r-z ; gr = g-r

                    # All
                    if index<2:
                        ib = bands[index]
                        mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                        icut = mcuts[index]
                        ind = np.where((mag<icut) & (lum_ext>lcut))
                    elif index==2: #eBOSS
                        ind = np.where((g>22.1) & (g<22.8) & \
                                           (gr>0.3) & (gr<0.7) & \
                                           (rz>0.25) & (rz<1.4) & \
                                           (rz>0.5*gr+0.4) & \
                                           (rz<0.5*gr+0.8) & \
                                           (lum_ext>lcut))
                    elif index==3: #DESI
                        ind = np.where((r<23.4) &\
                                           (rz>0.3) & (gr>-0.3) &\
                                           (rz>0.9*gr+0.12) &\
                                           (rz<1.345-0.85*gr) &\
                                           (lum_ext>lcut))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        nm[index,:] = nm[index,:] + H

                    # Centrals
                    if index<2:
                        ib = bands[index]
                        mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                        icut = mcuts[index]
                        ind = np.where((mag<icut) &\
                                           (lum_ext>lcut)  & (gtype==0))
                    elif index==2: #eBOSS
                        ind = np.where((g>22.1) & (g<22.8) & \
                                           (gr>0.3) & (gr<0.7) & \
                                           (rz>0.25) & (rz<1.4) & \
                                           (rz>0.5*gr+0.4) & \
                                           (rz<0.5*gr+0.8) & \
                                           (lum_ext>lcut)  & (gtype==0))
                    elif index==3: #DESI
                        ind = np.where((r<23.4) &\
                                           (rz>0.3) & (gr>-0.3) &\
                                           (rz>0.9*gr+0.12) &\
                                           (rz<1.345-0.85*gr) &\
                                           (lum_ext>lcut)  & (gtype==0))

                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        nc[index,:] = nc[index,:] + H

                    # Satellites
                    if index<2:
                        ib = bands[index]
                        mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                        icut = mcuts[index]
                        ind = np.where((mag<icut) &\
                                           (lum_ext>lcut)  & (gtype>0))
                    elif index==2: #eBOSS
                        ind = np.where((g>22.1) & (g<22.8) & \
                                           (gr>0.3) & (gr<0.7) & \
                                           (rz>0.25) & (rz<1.4) & \
                                           (rz>0.5*gr+0.4) & \
                                           (rz<0.5*gr+0.8) & \
                                           (lum_ext>lcut)  & (gtype>0))
                    elif index==3: #DESI
                        ind = np.where((r<23.4) &\
                                           (rz>0.3) & (gr>-0.3) &\
                                           (rz>0.9*gr+0.12) &\
                                           (rz<1.345-0.85*gr) &\
                                           (lum_ext>lcut)  & (gtype>0))

                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        ns[index,:] = ns[index,:] + H

                    # Type 2
                    if index<2:
                        ib = bands[index]
                        mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                        icut = mcuts[index]
                        ind = np.where((mag<icut) &\
                                           (lum_ext>lcut)  & (gtype==2))
                    elif index==2: #eBOSS
                        ind = np.where((g>22.1) & (g<22.8) & \
                                           (gr>0.3) & (gr<0.7) & \
                                           (rz>0.25) & (rz<1.4) & \
                                           (rz>0.5*gr+0.4) & \
                                           (rz<0.5*gr+0.8) & \
                                           (lum_ext>lcut)  & (gtype==2))
                    elif index==3: #DESI
                        ind = np.where((r<23.4) &\
                                           (rz>0.3) & (gr>-0.3) &\
                                           (rz>0.9*gr+0.12) &\
                                           (rz<1.345-0.85*gr) &\
                                           (lum_ext>lcut)  & (gtype==2))

                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        n2[index,:] = n2[index,:] + H

                f.close()

    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Plot
    fig = plt.figure(figsize=(8.,9.))    
    ax = plt.subplot()
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    ax.text(xmax-(xmax-xmin)*0.2, ymin+(ymax-ymin)*0.05, 'z='+szz)

    # Plot the model predictions
    leg2 = []
    for index,ileg in enumerate(inleg):
        indh = np.where(nh > 0) ; nhalos = sum(nh)
        if (np.shape(indh)[1]<1): continue
        print index,ileg,np.shape(indh)[1]

        py = 0. ; py = nm[index,:] ; nall = sum(py)
        x = lhist[indh] ; y = np.zeros(len(x))
        yall = py[indh]/nh[indh] ; ind = np.where(yall>0.) 
        y[ind] = np.log10(yall[ind]) ; ind = np.where(y < 0.)
        ax.plot(x[ind],y[ind],color=cols[index]\
                             ,label=inleg[index])

        if (index==0):
            py = 0. ; py = nc[index,:] ; ncen = sum(py)
            x = lhist[indh]  ; y = np.zeros(len(x))
            ycen = py[indh]/nh[indh] ; ind = np.where(ycen>0.) 
            y[ind] = np.log10(ycen[ind]) ;ind = np.where(y < 0.)
            l1, = plt.plot(x[ind],y[ind],color=cols[index],\
                         linestyle='--')

            py = 0. ; py = ns[index,:]  ; nsat = sum(py)
            x = lhist[indh] ; y = np.zeros(len(x))
            ysat = py[indh]/nh[indh] ; ind = np.where(ysat>0.) 
            y[ind] = np.log10(ysat[ind]) ; ind = np.where(y < 0.)
            l2, = plt.plot(x[ind],y[ind],color=cols[index],\
                                  linestyle=':')

            leg2.append([l1,l2])

            #py = 0. ; py = n2[index,:]  ; nsat2 = sum(py)
            #x = lhist[indh] ; y = np.zeros(len(x))
            #y2 = py[indh]/nh[indh] ; ind = np.where(y2>0.) 
            #y[ind] = np.log10(y2[ind]) ; ind = np.where(y < 0.)
            #l3, = plt.plot(x[ind],y[ind],color=cols[index],\
            #                      linestyle='-.')
            #leg2.append([l3])

    # Legends
    #first_legend = plt.legend(leg2[0],\
    #                              ["Centrals", "Satellites", "Orphans"],\
    #                              loc=1)
    first_legend = plt.legend(leg2[0],\
                                  ["Centrals", "Satellites"],\
                                  loc=1)
    first_legend.draw_frame(False)

    leg = plt.legend(loc=2,\
                         handlelength=0, handletextpad=0)
    for item in leg.legendHandles:
        item.set_visible(False) 
    for color,text in zip(cols,leg.get_texts()):   
        text.set_color(color)
        leg.draw_frame(False)
    plt.gca().add_artist(first_legend)

    # Save figures
    plotfile = outdir+'plots/cuts/'+model+'hod_'+line+\
        '_z'+str(zsnap)+'.pdf'
    plt.savefig(plotfile)
    print 'Output: ',plotfile
