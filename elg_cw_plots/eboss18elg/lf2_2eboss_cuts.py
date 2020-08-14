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

path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp19/'

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/'+model+'lfs/lf2_ebossmod'
plotfile = outdir+line+'.pdf'
############################# Obs
obsh0 = 0.677
obs_dir = '/cosma5/data/durham/violeta/lines/desi_hod_o2/lf_obs_data/'
#############################

snap_list = [41,39] #MillGas
nvol = 64

#####
obsnom = ['DEEP2','VVDSDEEP','VVDSWIDE']
obands = ['R24.2','I24','I22.5']

bands = ['DEIMOS-R','MegaCam-i-atmos','MegaCam-i-atmos']
mcuts = [24.1, 24, 22.5]
fcuts = [2.7*10.**-17., 1.9*10.**-17., 3.5*10.**-17.,8.*10.**-17.,10.**-16.,10.**-16.]

inleg = ['All','DEEP2','VVDS-DEEP','VVDS-Wide','DESI','eBOSS-SGC','eBOSSmod']
##########

ntypes = len(inleg)
zleg = []
cols = get_distinct(ntypes-2)
cols.insert(0,'grey')

# Initialize histogram
lmin = 38.
lmax = 46.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
fig = plt.figure(figsize=(6.5,14.))

xtit = "${\\rm log}_{10}(L\\rm{"+lline+"}/h^{-2}erg\, s^{-1})$"
ytit = "${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dex}^{-1})$"

xmin = 40.2 ; xmax = 43.7
ymin = -5.9 ; ymax = -1.

# Loop over the redshifts of interest
jj = 410
for iz,zsnap in enumerate(snap_list):
    print('\n ####### snap={} ####### \n'.format(zsnap))
    jj = jj + 1

    lf = np.zeros(shape=(ntypes,len(lhist)))
    lf_ext = np.zeros(shape=(ntypes,len(lhist)))

    volume = 0. ; firstpass = True
    for ivol in range(nvol):
        gfile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
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
            f.close()

            if(firstpass):
                zstring = "{:.2f}".format(zz) 
                zleg.append('z='+zstring)
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
                    elif ((inleg[index] == 'eBOSS-SGC') or
                          (inleg[index] == 'DESI') or
                          (inleg[index] == 'eBOSSmod')):
                        fluxcut = fcuts[index-1]
                        lcut = emission_line_luminosity(fluxcut,zz)
                        if (firstpass):
                            print('{}: flux_cut={}, lcut={:.3f}'.format(inleg[index],
                                                                    fluxcut,np.log10(lcut)+40.))

                        g = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                        r = f['Output001/mag_DES-r_o_tot_ext'].value + tomag
                        z = f['Output001/mag_DES-z_o_tot_ext'].value + tomag
                        rz = r-z ; gr = g-r

                        g_int = f['Output001/mag_DES-g_o_tot'].value + tomag
                        r_int = f['Output001/mag_DES-r_o_tot'].value + tomag
                        z_int = f['Output001/mag_DES-z_o_tot'].value + tomag
                        rz_int = r_int-z_int ; gr_int = g_int-r_int

                        if (inleg[index] == 'DESI'): 
                            ind  = np.where((r<23.4) & \
                                                (rz>0.3) & (gr>-0.3) & \
                                                (rz>0.9*gr+0.12) & \
                                                (rz<1.345-0.85*gr) & \
                                                (lum_ext>lcut))

                            indi  = np.where((r_int<23.4) & \
                                                (rz_int>0.3) & (gr_int>-0.3) & \
                                                (rz_int>0.9*gr_int+0.12) & \
                                                (rz_int<1.345-0.85*gr_int) & \
                                                (lum>lcut))

                        elif (inleg[index] == 'eBOSS-SGC'): 
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

#                            indi = np.where((lum>lcut) & \
#                                            (g_int>21.825) & (g_int<22.825) & \
#                                            (gr_int>-0.068*rz_int + 0.457) & \
#                                            (gr_int<0.112*rz_int + 0.773) & \
#                                            (rz_int>0.218*gr_int + 0.571) & \
#                                            (rz_int<-0.555*gr_int + 1.901))
#
                        elif (inleg[index] == 'eBOSSmod'): 
                            ind = np.where((lum_ext>lcut) & \
                                           (g>21.825) & (g<22.825) & \
                                           (gr>-0.068*rz + 0.457) & \
                                           (gr<0.112*rz + 0.773) & \
                                           (rz>0.218*gr + 0.85) & \
                                           (rz<-0.555*gr + 1.901))
                            indi = np.where((lum>lcut) & \
                                           (g>21.825) & (g<22.825) & \
                                           (gr>-0.068*rz + 0.457) & \
                                           (gr<0.112*rz + 0.773) & \
                                           (rz>0.218*gr + 0.85) & \
                                           (rz<-0.555*gr + 1.901))

#                            indi = np.where((lum>lcut) & \
#                                            (g_int>21.825) & (g_int<22.825) & \
#                                            (gr_int>-0.068*rz_int + 0.457) & \
#                                            (gr_int<0.112*rz_int + 0.773) & \
#                                            (rz_int>0.218*gr_int + 0.85) & \
#                                            (rz_int<-0.555*gr_int + 1.901))

                    else:
                        fluxcut = fcuts[index-1]
                        lcut = emission_line_luminosity(fluxcut,zz)
                        if (firstpass):
                            print('{}: flux_cut={}, lcut={:.3f}'.format(inleg[index],
                                                                    fluxcut,np.log10(lcut)+40.))

                        icut = mcuts[index-1]
                        ib = bands[index-1]

                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value\
                            + tomag
                        ind  = np.where((mag<icut) & (lum_ext>lcut))

                        mag_int = f['Output001/mag_'+ib+'_o_tot'].value\
                            + tomag
                        indi = np.where((mag_int<icut) & (lum>lcut))

                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H

                    if (np.shape(indi)[1] > 0.):
                        ll = np.log10(lum[indi]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H
                f.close()
                firstpass = False

    for index in range(ntypes):        
        n_lf = sum(lf[index,:])
        n_lf_ext = sum(lf_ext[index,:])
        print('    {}: n_int={}, n_ext={}, n_int/n_ext={:.2f}'.format(inleg[index],
              n_lf,n_lf_ext,n_lf/n_lf_ext))

    lf = lf/dl/volume
    lf_ext = lf_ext/dl/volume
    print('\n Side of the explored box (Mpc/h) = {} \n'.format(pow(volume,1./3.)))

    # Plot
    if (iz == 0):
        ax1 = fig.add_subplot(jj) ; ax1.set_autoscale_on(False)
        ax1.set_xlim([xmin,xmax]) ; ax1.set_ylim([ymin,ymax]) 
        ax1.set_autoscale_on(False) ;  ax1.minorticks_on()
        ax1.set_xlabel(xtit) ; ax1.set_ylabel(ytit)
        #ax1.tick_params(labelsize=fs-2) 
        ax1.text(42.5, -1.7, zleg[iz])
    else:
        ax = fig.add_subplot(jj,sharex=ax1,sharey=ax1)
        ax.set_autoscale_on(False) ;  ax.minorticks_on()
        ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
        ax.text(42.5, -1.7, zleg[iz])

    # Plot all observations
    ox, oy, el, eh = jc.read_jc_lf(obs_dir+'lf_may16_comparat/',zz,h0=obsh0,\
                                       infile=\
                                       'O2_3728-data-summary-Planck15.txt')
    ind = np.where(oy>-5) 
    oxr = ox[ind] ; oyr = oy[ind]
    arrinds = oxr.argsort()
    oxr = oxr[arrinds]
    oyr = oyr[arrinds]

    if(isinstance(ox, (np.ndarray))):
        if (iz == 0):
            ax1.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                             ecolor='grey',color='grey',mec='grey')
        else:
            ax.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                            ecolor='grey',color='grey',mec='grey')
    # Plot the observations  O2_3728-*-z*.txt
    for i,isurvey in enumerate(obsnom):        
        ox, oy, el, eh = jc.read_jc_indlf(obs_dir+'lf_may16_comparat/individual_LF/',\
                                              zz,h0=obsh0,\
                                              line='O2_3728',\
                                              survey=isurvey,band=obands[i])
        if(isinstance(ox, (np.ndarray))):
            col = cols[i+1]
            if (iz == 0):
                ax1.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                                 ecolor=col,color=col,mec=col)
            else:
                ax.errorbar(ox,oy,yerr=[el,eh],fmt='o',\
                                ecolor=col,color=col,mec=col)

    # Plot Prabhakar Tiwari's LF
    if (zsnap == 41):
        oxh, oyh, oeh = np.loadtxt(obs_dir+'tiwari/OII_3728_LF_v11.dat',
                                 usecols=(0,1,2),unpack=True,skiprows=1)

        ox = oxh + 2*np.log10(h0)
        oy = oyh + 3*np.log10(obsh0) - 3*np.log10(h0)
        oe = oeh + 3*np.log10(obsh0) - 3*np.log10(h0)

        i = inleg.index('eBOSS-SGC') 

        if (iz == 0):
            ax1.errorbar(ox,oy,yerr=[oe,oe],fmt='o',\
                             ecolor=cols[i],color=cols[i],mec=cols[i])
        else:
            ax.errorbar(ox,oy,yerr=[oe,oe],fmt='o',\
                            ecolor=cols[i],color=cols[i],mec=cols[i])

    # Plot the model predictions
    for index in range(ntypes):
        # Attenuated
        py = 0. ; py = lf_ext[index,:]
        ind = np.where(py > 0)
        x = lhist[ind]
        y = np.log10(py[ind])
        ind = np.where(y < 0.)       
        if (iz == 0):
            if (index<ntypes-1):
                ax1.plot(x[ind],y[ind],color=cols[index],linestyle='-',\
                         label=inleg[index])
            else:
                ax1.plot(x[ind],y[ind],color=cols[index-1],linestyle='--')
        else:
            if (index<ntypes-1):
                ax.plot(x[ind],y[ind],color=cols[index],linestyle='-',\
                        label=inleg[index])
            else:
                ax.plot(x[ind],y[ind],color=cols[index-1],linestyle='--')

        # Ratios
        if (index == 0):
            my = np.interp(oxr,x,y) 
            diff = abs(my-oyr) ; ratio = 10.**(diff)
            print('Max. LF ratio = {}, Min(diff)={}, Max(diff)={}'.format(max(ratio),
                                                                       min(diff),max(diff)))
            print('x_lf={}, ratios={}'.format(oxr,ratio))

        # Intrinsic
        if (index == 0):
            py = 0. ; py = lf[index,:]
            ind = np.where(py > 0)
            x = lhist[ind]
            y = np.log10(py[ind])
            ind = np.where(y < 0.)
            if (iz == 0):
                ax1.plot(x[ind],y[ind],color=cols[index],linestyle=':')
            else:
                ax.plot(x[ind],y[ind],color=cols[index],linestyle=':')

        # Intrinsic eBOSS
        if (index == 5):
            py = 0. ; py = lf[index,:]
            ind = np.where(py > 0)
            x = lhist[ind]
            y = np.log10(py[ind])
            ind = np.where(y < 0.)
            if (iz == 0):
                ax1.plot(x[ind],y[ind],color=cols[index],linestyle=':')
            else:
                ax.plot(x[ind],y[ind],color=cols[index],linestyle=':')

        # Intrinsic eBOSSmod
        if (index == 6):
            py = 0. ; py = lf[index,:]
            ind = np.where(py > 0)
            x = lhist[ind]
            y = np.log10(py[ind])
            ind = np.where(y < 0.)
            if (iz == 0):
                ax1.plot(x[ind],y[ind],color=cols[index-1],linestyle='-.')
            else:
                ax.plot(x[ind],y[ind],color=cols[index-1],linestyle='-.')

        # Legend
        if (iz == len(snap_list)-1):
            leg = plt.legend(loc=1, handlelength=0, handletextpad=0)
            renderer = fig.canvas.get_renderer()
            shift1 = max([t.get_window_extent(renderer).width for t in leg.get_texts()])
            shift2 = min([t.get_window_extent(renderer).width for t in leg.get_texts()])
            for item in leg.legendHandles:
                item.set_visible(False)
            for color,text in zip(cols,leg.get_texts()):
                text.set_color(color)
                text.set_ha('right') ; text.set_position((shift1-shift2,0))
                leg.draw_frame(False)           


#On splitting the legend
#line1, = plt.plot([1,2,3], label="Line 1", linestyle='--')
#line2, = plt.plot([3,2,1], label="Line 2", linewidth=4)
#
## Create a legend for the first line.
#first_legend = plt.legend(handles=[line1], loc=1)
#
## Add the legend manually to the current Axes.
#ax = plt.gca().add_artist(first_legend)
#
## Create another legend for the second line.
#plt.legend(handles=[line2], loc=4)
    
# Save figures
fig.subplots_adjust(hspace=0)
fig.savefig(plotfile)
print 'Output: ',plotfile
