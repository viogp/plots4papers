#! /usr/bin/env python

import numpy as np
import os.path, sys
import h5py
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import read_jc_obs as jc
from Cosmology import * 

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
model = str(sys.argv[1])

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/cuts/'+model+'att_compared_'
plotfile = outdir+line+'_deep2.pdf'
############################# Obs
obsh0 = 0.677
obs_dir = '/gpfs/data/violeta/lines/desi_hod_o2/lf_obs_data/lf_may16_comparat/individual_LF/'
#############################

snap_list = [61, 44, 42, 40, 37, 34] #MillGas #z=0, 0.6, 0.75, 0.9, 1.18
zleg = ['z=0.','z=0.6','z=0.75', 'z=0.9', 'z=1.18', 'z=1.5']
nvol = 64

obsnom = 'DEEP2' ; obands = 'R24.2' ; band = 'RK' ; icut = 24.1 ; fluxcut = 2.7*10.**-17. 

#obsnom = 'VVDSWIDE' ; obands = 'I22.5' ; band = 'm2' ; icut = 22.5 ; fluxcut = 3.5*10.**-17. 

#obsnom = 'VVDSDEEP' ; obands = 'I24' ; band = 'm2' ; icut = 24. ; fluxcut = 1.9*10.**-17.

inleg1 = ['Attenuated','Centrals','Quiescent','High sSFR','Big','Spheroid']
inleg2 = ['Intrinsic','Satellites','Bursty','Low sSFR','Small','Disk']
colors = ['k','r','g','b','c','m']  ; colors12 = ['k','k','r','r','g','g','b','b','c','c','m','m']

# Initialize histogram
lmin = 38.
lmax = 46.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15
fig = plt.figure(figsize=(14.,10.))

xtit = "${\\rm log}_{10}(L\\rm{"+lline+"}/h^{-2}erg\, s^{-1})$"
ytit = "${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dlog}_{10}L)$"

xmin = 40. ; xmax = 44.
ymin = -5.5 ; ymax = -1.

print 'Start'
# Loop over the redshifts of interest
jj = 230
for iz,zsnap in enumerate(snap_list):
    jj = jj + 1

    lf = np.zeros(shape=(len(inleg1),len(lhist)))
    lf_ext = np.zeros(shape=(len(inleg1),len(lhist)))

    volume = 0.
    ntot=0. ; ntotsf=0. ; no2=0. ; ncen=0. ; nsat=0. ; nli=0. 
    nqui=0. ; nbur=0. ; nsfh=0. ; nsfl=0. ; nsfhm=0. 
    nbig=0. ; nsmall=0. ; nspheroid=0. ; ndisk=0.
    mmin = 999. ; mmean =0. ; firstpass = True
    for ivol in range(nvol):
        gfile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (os.path.isfile(gfile)):
            #print gfile
            # Get some of the model constants
            f = h5py.File(gfile,'r')
            group = f['Parameters']
            vol1 = group['volume'].value ; volume = volume + vol1
            h0 = group['h0'].value 
            omega0 = group['omega0'].value
            omegab = group['omegab'].value
            lambda0 =group['lambda0'].value
            zz = f['Output001/redshift'].value
            tsf    = 3*f['Output001/tsfburst'].value 
            tburst = f['Output001/tburst'].value
            r50 = np.log10(f['Output001/rcomb'].value) +3. #kpc/h 
            sfrd = f['Output001/mstardot'].value
            sfrb = f['Output001/mstardot_burst'].value
            sfr = sfrd+sfrb
            mass = f['Output001/mstars_disk'].value +\
                f['Output001/mstars_bulge'].value
            mhhalo = f['Output001/mhhalo'].value
            lssfr = np.zeros(shape=(len(mass))) ; lssfr.fill(-999.)
            ind = np.where((sfr>0.) & (mass>0.))
            ntot = ntot + np.size(ind)
            lssfr[ind] = np.log10(sfr[ind]) - np.log10(mass[ind])

            set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                              h0=h0, universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz) 
            slim = 0.3/tHubble(zz) # Franx+08

            ind  = np.where((sfr>0.) & (mass>0.) & (10.**lssfr>slim))
            ntotsf = ntotsf + np.size(ind)

            f.close()

            efile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r')
                BoT = f['Output001/BoT'].value
                lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                lum = f['Output001/L_tot_'+line].value

                if(len(lum_ext) ==len(sfr)):
                    mag = f['Output001/mag'+band+'o_tot_ext'].value + tomag
                    lcut = emission_line_luminosity(fluxcut,zz)

                    index =0
                    # All
                    ind  = np.where((mag<icut) & (lum_ext>lcut))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H
                        no2 = no2 + np.shape(ind)[1]

                        m = np.log10(mhhalo[ind])
                        mmin1 = np.min(m)
                        mmean1= np.mean(m)
                        if firstpass:
                            mmin = mmin1
                            mmean= mmean1
                            firstpass = False
                        else:
                            mmin = np.min([mmin,mmin1])
                            mmean= np.mean([mmean,mmean1])

                    # Ext
                    indi = np.where((mag<icut) & (lum>lcut))
                    if (np.shape(indi)[1] > 0.):
                        ll = np.log10(lum[indi]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H
                        nli = nli + np.shape(indi)[1]

                    index = 1
                    # Centrals
                    cen = f['Output001/type'].value
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (cen<1))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H
                        ncen = ncen + np.shape(ind)[1]

                    # Satellites
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (cen>0))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H
                        nsat = nsat + np.shape(ind)[1]

                    index = 2
                    # Quiescent
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (tburst>tsf))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H
                        nqui = nqui + np.shape(ind)[1]
                        
                    # Bursty
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (tburst<tsf))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H
                        nbur = nbur + np.shape(ind)[1]
                                        
                    index = 3
                    # High sSFR
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (10.**lssfr>slim))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H
                        nsfh = nsfh + np.shape(ind)[1]

                    ind  = np.where((lum_ext>lcut) & (10.**lssfr>slim))
                    nsfhm = nsfhm + np.shape(ind)[1]

                    # Low sSFR
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (10.**lssfr<slim))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H
                        nsfl = nsfl + np.shape(ind)[1]
                        
                    index = 4
                    # Big.
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (r50>0.5))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H
                        nbig = nbig + np.shape(ind)[1]

                    # Small.
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (r50<=0.5))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H
                        nsmall = nsmall + np.shape(ind)[1]

                    index = 5
                    # Spheroid.
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (BoT>0.5))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf_ext[index,:] = lf_ext[index,:] + H
                        nspheroid = nspheroid + np.shape(ind)[1]

                    # Disk.
                    ind  = np.where((mag<icut) & (lum_ext>lcut) & (BoT<=0.5))
                    if (np.shape(ind)[1] > 0.):
                        ll = np.log10(lum_ext[ind]) + 40.
                        #ll = np.log10(lum[ind]) + 40.
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        lf[index,:] = lf[index,:] + H
                        ndisk = ndisk + np.shape(ind)[1]

                else:
                    print 'WARNING: arrays do not match in ivol=',ivol
                f.close()
            else:
                print 'NOT found',efile
        else:
            print 'NOT found',gfile

    lf = lf/dl/volume
    lf_ext = lf_ext/dl/volume
    print '---------------------------------------------'
    print 'z= ',zz,', Boundary sSFR = ',slim
    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)
    print ' '
    print '  mhalo_min= ',mmin,'  mhalo_mean= ',mmean
    print ' '
    print 'ntot=',ntot,' , ntotsf =',ntotsf
    print ' '
    print 'no2=',no2,' , nli =',nli
    print 'no2(%)=',no2,' (',no2*100./ntot,\
        '), nli(%)=',nli,' (',nli*100./ntot,')'
    print 'no2/volume=',no2/volume
    print ' '
    if (no2>0.):
        print 'ncen(%)=',ncen,' (',ncen*100./no2,\
            '), nsat(%)=',nsat,' (',nsat*100./no2,') , ncen+nsat=',ncen+nsat
        print ' '
        print 'nqui(%)=',nqui,' (',nqui*100./no2,\
            '), nbur(%)=',nbur,' (',nbur*100./no2,') , nqui+nbur=',nqui+nbur
        print ' '
        print 'nsfh(%)=',nsfh,' (',nsfh*100./no2,\
            '), nsfl(%)=',nsfl,' (',nsfl*100./no2,') , nsfh+nsfl=',nsfh+nsfl
        print 'nsfh(%total sf)=',nsfh*100./ntotsf, \
            'nsfhm(%total sf)=',nsfhm,' (',nsfhm*100./ntotsf,')'
        print ' '
        print 'nbig(%)=',nbig,' (',nbig*100./no2,\
            '), nsmall(%)=',nsmall,' (',nsmall*100./no2, \
            ') , nbig+nsmall=',nbig+nsmall
        print ' '
        print 'nspheroid(%)=',nspheroid,' (',nspheroid*100./no2,\
            '), ndisk(%)=',ndisk,' (',ndisk*100./no2, \
            ') , nspheroid+ndisk=',nspheroid+ndisk
    print '---------------------------------------------'

    # Plot
    ax = fig.add_subplot(jj)
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    ax.set_xlabel(xtit,fontsize = fs) ; ax.set_ylabel(ytit,fontsize = fs)
    ax.tick_params(labelsize=13)   
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 1.))
    ax.xaxis.set_ticks(np.arange(start, end, 0.2),minor=True)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    ax.text(40.1, -1.3, zleg[iz])

    # Plot the observations  O2_3728-VVDSWIDEI21.5-z1.234.txt
    ox, oy, el, eh = jc.read_jc_indlf(obs_dir,zz,h0=obsh0,line='O2_3728',survey=obsnom,band=obands)
    if hasattr(ox, "__len__"):
        ax.errorbar(ox,oy,yerr=[el,eh],fmt='o',ecolor='grey',color='grey',mec='grey')

    # Plot the model predictions
    for index in range(len(inleg1)):
        py1 = 0. ; py1 = lf_ext[index,:]
        ind = np.where(py1 > 0)
        x1 = lhist[ind] 
        y1 = np.log10(py1[ind])
        ind = np.where(y1 < 0.)
        ax.plot(x1[ind],y1[ind],color=colors[index],linestyle='-',label=inleg1[index])

        py = 0. ; py = lf[index,:]
        ind = np.where(py > 0)
        x = lhist[ind]
        y = np.log10(py[ind])
        ind = np.where(y < 0.)
        ax.plot(x[ind],y[ind],color=colors[index],linestyle='--',label=inleg2[index])

        #ax.plot(lhist,np.log10(py1+py),color=colors[index],linestyle='-.')
                
        # Legend
        leg = plt.legend(loc=1,prop={'size':10})
        for color,text in zip(colors12,leg.get_texts()):
            text.set_color(color)
            leg.draw_frame(False)

# Save figures
fig.tight_layout()
fig.savefig(plotfile)
print 'Output: ',plotfile
