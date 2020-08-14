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

#Survey cuts
zleg = []
bands = ['DEIMOS-R','MegaCam-i-atmos','eBOSS','DESI']
mcuts = [24.1,22.5]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 10.**-16.,8.*10.**-17.]
inleg = ['All','DEEP2 cuts','VVDS-Wide cuts','eBOSS cuts','DESI cuts']
############

#snap_list = [44, 42, 40, 37, 34] #MillGas
snap_list = [42, 40]
nvol = 3 #64

#path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
#model = 'MillGas/gp15newmg.anders/' #.anders 'MillGas/gp15newmg/'

path = '/gpfs/data/violeta/Galform_Out/v2.7.0/stable/'
#model = 'MillGas/gp17/'  #model = 'MillGas/gp17.spin/' 
model = 'MillGas/gp17.spin.ramt0.01/'

line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/desi_hod_o2/'

ntypes = len(inleg) 
cols = get_distinct(ntypes-1) ; cols.insert(0,'k')

# Initialize GSMF
mmin = 8.
mmax = 15.
dm = 0.1
mbins = np.arange(mmin,mmax,dm)
mhist = mbins + dm*0.5
# Initialize SSFR
smin = 3.
smax = 12.
ds = 0.1
sbins = np.arange(smin,smax,ds)
shist = sbins + ds*0.5

nlevel = 10
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
for iz,zsnap in enumerate(snap_list):

    gsmf = np.zeros(shape=(ntypes,len(mhist)))
    ssfrf = np.zeros(shape=(ntypes,len(shist)))
    smf = np.zeros(shape=(ntypes,len(shist),len(mhist)))

    volume = 0. ; firstpass = True
    for ivol in range(nvol):
        gfile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if(os.path.isfile(gfile)):
            # Read the relevant information from galaxies.hdf5
            f = h5py.File(gfile,'r') #; print gfile
            vol1 = f['Parameters/volume'].value ; volume = volume + vol1
            h0 =   f['Parameters/h0'].value 
            omega0 = f['Parameters/omega0'].value
            omegab = f['Parameters/omegab'].value
            lambda0 =f['Parameters/lambda0'].value
            zz   = f['Output001/redshift'].value
            tsf    = 3*f['Output001/tsfburst'].value 
            tburst = f['Output001/tburst'].value
            cen = f['Output001/type'].value
            f.close()

            if(firstpass):
                szz = "{:.2f}".format(zz) 
                inleg[0] = 'All, z='+szz
                firstpass = False

            set_cosmology(omega0=omega0,omegab=omegab, \
                              lambda0=lambda0,h0=h0, \
                              universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz)

            efile = path+model+'iz'+str(zsnap)+\
                '/ivol'+str(ivol)+'/elgs.hdf5'
            if(os.path.isfile(efile)):
                f = h5py.File(efile,'r') #; print efile
                group = f['Output001']
                sfr = group['mstardot'].value +\
                    group['mstardot_burst'].value 
                mass =group['mstars_tot'].value
                lum_ext = group['L_tot_'+line+'_ext'].value

                # GSMF
                ind = np.where(mass>mmin) 
                ll = np.log10(mass[ind]) 
                H, bins_edges = np.histogram(ll,bins=np.append(mbins,mmax))
                gsmf[0,:] = gsmf[0,:] + H

                for ii,ib in enumerate(bands):
                    index = ii + 1
                    fluxcut = fcuts[ii]
                    lcut = emission_line_luminosity(fluxcut,zz)

                    if ii<2:
                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value + tomag
                        icut = mcuts[ii] 
                        ind  = np.where((mag<icut) & \
                                            (lum_ext>lcut) &\
                                            (mass>mmin) & (sfr>smin))
                    else:
                        g = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                        r = f['Output001/mag_DES-r_o_tot_ext'].value + tomag
                        z = f['Output001/mag_DES-z_o_tot_ext'].value + tomag
                        rz = r-z ; gr = g-r

                        if ii==2: #eBOSS decam180 selection
                            ind = np.where((g>22.1) & (g<22.8) & \
                                               (gr>0.3) & (gr<0.7) & \
                                               (rz>0.25) & (rz<1.4) & \
                                               (rz>0.5*gr+0.4) & \
                                               (rz<0.5*gr+0.8) & \
                                               (lum_ext>lcut) &\
                                               (mass>mmin) & (sfr>smin))

                        elif ii==3: #DESI
                            ind  = np.where((r<23.4) & \
                                                (rz>0.3) & (gr>-0.3) & \
                                                (rz>0.9*gr+0.12) & \
                                                (rz<1.345-0.85*gr) & \
                                                (lum_ext>lcut) &\
                                                (mass>mmin) & (sfr>smin))
                    ll = np.log10(mass[ind]) 
                    H, bins_edges = np.histogram(ll,bins=np.append(mbins,mmax))
                    gsmf[index,:] = gsmf[index,:] + H


                # SFR
                ind = np.where((mass>mmin) & (sfr>smin)) 
                lmt = np.log10(mass[ind]) 
                llt = np.log10(sfr[ind]) #- lmt
                H, bins_edges = np.histogram(llt,bins=np.append(sbins,smax))
                ssfrf[0,:] = ssfrf[0,:] + H

                # SFR-GSMF ########
                ind = np.where((mass>mmin) & (sfr>smin)) 
                lst = np.log10(sfr[ind]) #- np.log10(mass[ind]) 
                lmt = np.log10(mass[ind]) 
                H, xedges, yedges = np.histogram2d(lst,lmt,\
                                                       bins=[np.append(sbins,smax),np.append(mbins,mmax)])
                smf[0,:,:] = smf[0,:,:] + H


                for ii,ib in enumerate(bands):
                    index = ii + 1
                    fluxcut = fcuts[ii]
                    lcut = emission_line_luminosity(fluxcut,zz)

                    if ii<2:
                        mag = f['Output001/mag_'+ib+'_o_tot_ext'].value + tomag
                        icut = mcuts[ii] 
                        ind  = np.where((mag<icut) & \
                                            (lum_ext>lcut) & \
                                            (mass>mmin) & (sfr>smin))
                    else:
                        g = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                        r = f['Output001/mag_DES-r_o_tot_ext'].value + tomag
                        z = f['Output001/mag_DES-z_o_tot_ext'].value + tomag
                        rz = r-z ; gr = g-r

                        if ii==2: #eBOSS decam180 selection
                            ind = np.where((g>22.1) & (g<22.8) & \
                                               (gr>0.3) & (gr<0.7) & \
                                               (rz>0.25) & (rz<1.4) & \
                                               (rz>0.5*gr+0.4) & \
                                               (rz<0.5*gr+0.8) & \
                                               (lum_ext>lcut) &\
                                               (mass>mmin) & (sfr>smin))

                        elif ii==3: #DESI
                            ind  = np.where((r<23.4) & \
                                                (rz>0.3) & (gr>-0.3) & \
                                                (rz>0.9*gr+0.12) & \
                                                (rz<1.345-0.85*gr) & \
                                                (lum_ext>lcut) & \
                                                (mass>mmin) & (sfr>smin))

                    lmt = np.log10(mass[ind]) 
                    llt = np.log10(sfr[ind]) #- lmt
                    H, bins_edges = \
                        np.histogram(llt,bins=np.append(sbins,smax))
                    ssfrf[index,:] = ssfrf[index,:] + H

                    lst = np.log10(sfr[ind]) #- np.log10(mass[ind]) 
                    lmt = np.log10(mass[ind]) 
                    H, xedges, yedges = np.histogram2d(lst,lmt,\
                                                       bins=[np.append(sbins,smax),np.append(mbins,mmax)])
                    smf[index,:,:] = smf[index,:,:] + H


                f.close()
            else:
                print 'NOT found:',efile
        else:
            print 'NOT found:',gfile
            
    if (volume<=0.):
        sys.exit('EXIT: No volume read.')

    gsmf = gsmf/volume/dm  # In Mpc^3/h^3
    ssfrf = ssfrf/volume/ds  
    smf = smf/volume/dm/ds  
    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Figure http://matplotlib.org/users/gridspec.html
    fig = plt.figure(figsize=(8.5,9.))
    gs = gridspec.GridSpec(3, 3)
    gs.update(wspace=0., hspace=0.)
    ax = plt.subplot(gs[1:,:-1])

    # SFRF vs M
    xtit="$log_{10}(\\rm M_{*}/M_{\odot}h^{-1})$"
    ytit="$log_{10}(\\rm SFR/M_{\odot}h^{-1}Gyr^{-1})$"
    xmin=8.5 ; xmax=11.9 ; ymin=3. ; ymax=11.9
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
    for ii in range(ntypes):
        zz = np.zeros(shape=(len(shist),len(mhist))) 
        pz = smf[ii,:,:] ; ind = np.where(pz>0.)
        zz = np.log10(pz) 
        #zz[ind] = convert_to_stdev(np.log10(pz[ind]))

        # Heat map
        #x,y = np.meshgrid(np.append(mbins,mmax),np.append(sbins,smax))
        #ax.pcolormesh(x, y, zz, cmap=plt.get_cmap(ptmap))
        # Contours
        xx,yy = np.meshgrid(mbins,sbins)
        #zz = 100.*zz ; al = [5.,16.,68.,95.] 
        al = [-4.5,-3.5,-2.5,-1.5] 
        # Smooth
        #zz = ndimage.gaussian_filter(zz, sigma=0.5)
        cs = ax.contour(xx, yy, zz, levels=al, \
                            linestyles='-',colors=cols[ii])
        cs.levels = [nf(val) for val in cs.levels]
        ax.clabel(cs, cs.levels, inline=1,inline_spacing=0,\
                      fontsize=10,fmt='%r')#fmt='%r %%')


    # GSMF
    axm = plt.subplot(gs[0, :-1],sharex=ax)
    ytit="$log_{10}(\Phi)$" ; axm.set_ylabel(ytit)
    axm.set_autoscale_on(False) ;  axm.minorticks_on()
    axm.set_ylim(-5.5,-1.) 
    plt.setp(axm.get_xticklabels(), visible=False)
    for ii in range(ntypes):
        py = gsmf[ii,:] ; ind = np.where(py>0.)
        x = mhist[ind] ; y = np.log10(py[ind])
        ind = np.where(y < 0.)
        axm.plot(x[ind],y[ind],color=cols[ii],\
                    label=inleg[ii])


    # SFRF
    axs = plt.subplot(gs[1:, 2],sharey=ax)
    xtit="$log_{10}(\Phi)$" ; axs.set_xlabel(xtit)
    axs.set_autoscale_on(False) ;  axs.minorticks_on()
    axs.set_xlim(-5.5,0.) 
    plt.setp(axs.get_yticklabels(), visible=False)
    for ii in range(ntypes):
        px = ssfrf[ii,:] ; ind = np.where(px>0.)
        y = shist[ind] ; x = np.log10(px[ind])
        ind = np.where(x < 0.)
        axs.plot(x[ind],y[ind],color=cols[ii])

    # Legend
    leg = axm.legend(bbox_to_anchor=(1., 1.),fontsize='small',\
                         handlelength=0,handletextpad=0)
    for item in leg.legendHandles:
        item.set_visible(False)
    for color,text in zip(cols,leg.get_texts()):
        text.set_color(color)
        leg.draw_frame(False)


    # Save figure
    plotfile = outdir+'plots/cuts/'+model+'sfr_f3/'+line+\
        '_sn'+str(zsnap)+'_all_novvdsdeep.pdf'
    fig.savefig(plotfile)
    print 'Output: ',plotfile
