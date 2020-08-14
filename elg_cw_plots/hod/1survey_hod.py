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

#path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
#model = 'MillGas/gp15newmg/' #model = 'MillGas/gp14/'

path = '/gpfs/data/violeta/Galform_Out/v2.7.0/stable/'
model = 'MillGas/gp17/'
#model = 'MillGas/gp17.spin/'

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/desi_hod_o2/'
#############################

snap_list = [44, 42, 40, 37, 34] #MillGas    

nvol = 64

bands = ['DEIMOS-R','MegaCam-i-atmos','MegaCam-i-atmos']
mcuts = [24.1,22.5,24]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17.,10.**-17., 8*10.**-17.]
surveys = ['DEEP2','VVDS-WIDE','VVDS-DEEP','eBOSS','DESI']
inleg = ['DEEP2','VVDS-Wide','VVDS-Deep','eBOSS','DESI']

# Initialize histogram
lmin = 8.5
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
ntypes = len(snap_list) 
cols = get_distinct(ntypes)

xtit = "${\\rm log}_{10}({\\rm M_{halo}}/M_{\odot}h^{-1})$"
ytit = "${\\rm log}_{10}(\\langle N\\rangle _{\\rm [OII]})$"

xmin = 10. ; xmax = 15.
ymin = -3. ; ymax = 2.

# Loop over the redshifts of interest
for index,fc in enumerate(fcuts):
    nm = np.zeros(shape=(ntypes,len(lhist)))
    nc = np.zeros(shape=(ntypes,len(lhist)))
    ns = np.zeros(shape=(ntypes,len(lhist)))
    n2 = np.zeros(shape=(ntypes,len(lhist)))
    nf = np.zeros(shape=(ntypes,len(lhist)))
    nh = np.zeros(shape=(len(lhist)))
    zleg = []
    
    for iz,zsnap in enumerate(snap_list):
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
                f.close()

                set_cosmology(omega0=omega0,omegab=omegab,\
                                  lambda0=lambda0, h0=h0,\
                                  universe="Flat",include_radiation=False)
                tomag = band_corrected_distance_modulus(zz)
                if(firstpass):
                    szz = "{:.2f}".format(zz) 
                    zleg.append('z='+szz)
                    firstpass = False

                efile = path+model+'iz'+\
                    str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
                if (os.path.isfile(efile)):
                    f = h5py.File(efile,'r') #; print efile
                    # Number of haloes per bin
                    gtype  = f['Output001/type'].value
                    mhhalo = f['Output001/mhhalo'].value
                    lum_ext =f['Output001/L_tot_'+line+'_ext'].value
                    ind = np.where(gtype == 0)
                    ll = np.log10(mhhalo[ind])
                    H, bins_edges = \
                        np.histogram(ll,bins=np.append(lbins,lmax))
                    nh = nh + H
            
                    fluxcut = fc
                    lcut = emission_line_luminosity(fluxcut,zz)

                    g = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                    r = f['Output001/mag_DES-r_o_tot_ext'].value + tomag
                    z = f['Output001/mag_DES-z_o_tot_ext'].value + tomag
                    rz = r-z ; gr = g-r

                    # All
                    if index<3:
                        ib = bands[index]
                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value + tomag
                        icut = mcuts[index]
                        ind = np.where((mag<icut) & (lum_ext>lcut))
                    elif index==3: #eBOSS
                        ind = np.where((g>22.1) & (g<22.8) & \
                                           (gr>0.3) & (gr<0.7) & \
                                           (rz>0.25) & (rz<1.4) & \
                                           (rz>0.5*gr+0.4) & \
                                           (rz<0.5*gr+0.8) & \
                                           (lum_ext>lcut))
                    elif index==4: #DESI
                        ind = np.where((r<23.4) &\
                                           (rz>0.3) & (gr>-0.3) &\
                                           (rz>0.9*gr+0.12) &\
                                           (rz<1.345-0.85*gr) &\
                                           (lum_ext>lcut))
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        nm[iz,:] = nm[iz,:] + H

                    # Centrals
                    if index<3:
                        ib = bands[index]
                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value + tomag
                        icut = mcuts[index]
                        ind = np.where((mag<icut) &\
                                           (lum_ext>lcut)  & (gtype==0))
                    elif index==3: #eBOSS
                        ind = np.where((g>22.1) & (g<22.8) & \
                                           (gr>0.3) & (gr<0.7) & \
                                           (rz>0.25) & (rz<1.4) & \
                                           (rz>0.5*gr+0.4) & \
                                           (rz<0.5*gr+0.8) & \
                                           (lum_ext>lcut)  & (gtype==0))
                    elif index==4: #DESI
                        ind = np.where((r<23.4) &\
                                           (rz>0.3) & (gr>-0.3) &\
                                           (rz>0.9*gr+0.12) &\
                                           (rz<1.345-0.85*gr) &\
                                           (lum_ext>lcut)  & (gtype==0))

                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        nc[iz,:] = nc[iz,:] + H

                    # Satellites
                    if index<3:
                        ib = bands[index]
                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value + tomag
                        icut = mcuts[index]
                        ind = np.where((mag<icut) &\
                                           (lum_ext>lcut)  & (gtype>0))
                    elif index==3: #eBOSS
                        ind = np.where((g>22.1) & (g<22.8) & \
                                           (gr>0.3) & (gr<0.7) & \
                                           (rz>0.25) & (rz<1.4) & \
                                           (rz>0.5*gr+0.4) & \
                                           (rz<0.5*gr+0.8) & \
                                           (lum_ext>lcut)  & (gtype>0))
                    elif index==4: #DESI
                        ind = np.where((r<23.4) &\
                                           (rz>0.3) & (gr>-0.3) &\
                                           (rz>0.9*gr+0.12) &\
                                           (rz<1.345-0.85*gr) &\
                                           (lum_ext>lcut)  & (gtype>0))

                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        ns[iz,:] = ns[iz,:] + H

                    # Type 2
                    if index<3:
                        ib = bands[index]
                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value + tomag
                        icut = mcuts[index]
                        ind = np.where((mag<icut) &\
                                           (lum_ext>lcut)  & (gtype==2))
                    elif index==3: #eBOSS
                        ind = np.where((g>22.1) & (g<22.8) & \
                                           (gr>0.3) & (gr<0.7) & \
                                           (rz>0.25) & (rz<1.4) & \
                                           (rz>0.5*gr+0.4) & \
                                           (rz<0.5*gr+0.8) & \
                                           (lum_ext>lcut)  & (gtype==2))
                    elif index==4: #DESI
                        ind = np.where((r<23.4) &\
                                           (rz>0.3) & (gr>-0.3) &\
                                           (rz>0.9*gr+0.12) &\
                                           (rz<1.345-0.85*gr) &\
                                           (lum_ext>lcut)  & (gtype==2))

                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        n2[iz,:] = n2[iz,:] + H

                    # Only flux cut
                    ind = np.where(lum_ext>lcut)
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges =\
                            np.histogram(ll,bins=np.append(lbins,lmax))
                        nf[iz,:] = nf[iz,:] + H

                    f.close()

    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Plot
    fig = plt.figure(figsize=(8.,9.))    
    ax = plt.subplot()
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    #ax.text(xmax-(xmax-xmin)*0.25, ymin+(ymax-ymin)*0.02,\
    #            inleg[index]+' cuts')

    # Plot the model predictions
    leg1 = [] ; leg2 = []
    for iz,ileg in enumerate(zleg):
        #print iz,np.shape(nm)
        indh = np.where(nh > 0) ; nhalos = sum(nh)
        if (np.shape(indh)[1]<1): continue

        py = nm[iz,:] ; nall = sum(py)
        x = lhist[indh] 
        yall = py[indh]/nh[indh] 
        y = np.log10(yall)
        l1, = ax.plot(x,y,color=cols[iz],\
                             label=zleg[iz])

        py = nf[iz,:] ; nfcut = sum(py)
        x = lhist[indh] 
        yall = py[indh]/nh[indh] 
        y = np.log10(yall) 
        l2, = ax.plot(x,y,color=cols[iz],\
                    linestyle='--')

        leg1.append([l1,l2])

        #if (iz==ntypes-1):
        #    py = nc[iz,:] ; ncen = sum(py)
        #    x = lhist[indh] 
        #    ycen = py[indh]/nh[indh]
        #    y = np.log10(ycen)
        #    l1, = plt.plot(x,y,color=cols[iz],\
        #                 linestyle='--')
        #
        #    py = ns[iz,:]  ; nsat = sum(py)
        #    x = lhist[indh]
        #    ysat = py[indh]/nh[indh]
        #    y = np.log10(ysat) 
        #    l2, = plt.plot(x,y,color=cols[iz],\
        #                          linestyle=':')
        #
        #    py = n2[iz,:]  ; nsat2 = sum(py)
        #    x = lhist[indh]
        #    y2 = py[indh]/nh[indh]
        #    y = np.log10(y2)
        #    l3, = plt.plot(x,y,color=cols[iz],\
        #                          linestyle='-.')
        #    leg2.append([l1,l2,l3])

    # Legends
    #first_legend = plt.legend(leg2[0],\
    #                              ["Centrals", "Satellites", "Orphans"],\
    #                              loc=1)

    ltxt = inleg[index]+' cuts'
    first_legend = plt.legend(leg1[0],\
                                  [ltxt, "Flux cut"],\
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
        '_'+str(surveys[index])+'.pdf'
    plt.savefig(plotfile)
    print 'Output: ',plotfile
