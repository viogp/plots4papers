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

path = '/gpfs/data/violeta/Galform_Out/v2.7.0/stable/MillGas/'

snap_list = [44, 42, 41, 40, 39, 37, 34] #MillGas
nvol = 64

model = 'gp18' #'gp18.font','gp18.starvation','gp17'

surveys = ['eBOSS']#['DEEP2','VVDS-WIDE','VVDS-DEEP','eBOSS','DESI']

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/cosmicweb/'
plotfile = outdir+'plots/'+model+'/selections/gsmf_'+line+'_all.pdf'
#############################

ntypes = len(surveys)
cols = get_distinct(ntypes) 

# Initialize histogram
lmin = 8.
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15
fig = plt.figure(figsize=(9.,15.))

xtit = "${\\rm log}_{10}(M_{*}/M_{\odot}h^{-1})$"
ytit = "${\\rm log}_{10}(\Phi/h^3Mpc^{-3}dex^{-1})$"

xmin = 8.5 ; xmax = 12.5
ymin = -5.5 ; ymax = -0.5

# Loop over the redshifts of interest
jj = 420
for iz,zsnap in enumerate(snap_list):
    jj = jj + 1

    all_nm = np.zeros(shape=(len(lhist)))
    all_cen= np.zeros(shape=(len(lhist)))
    all_sat= np.zeros(shape=(len(lhist)))

    sf_nm = np.zeros(shape=(len(lhist)))
    sf_cen= np.zeros(shape=(len(lhist)))
    sf_sat= np.zeros(shape=(len(lhist)))

    nm = np.zeros(shape=(len(surveys),len(lhist)))
    cen= np.zeros(shape=(len(surveys),len(lhist)))
    sat= np.zeros(shape=(len(surveys),len(lhist)))

    o_all_nm = np.zeros(shape=(len(lhist))) ; o_all_nm.fill(-999.)
    o_all_cen= np.zeros(shape=(len(lhist))) ; o_all_cen.fill(-999.)
    o_all_sat= np.zeros(shape=(len(lhist))) ; o_all_sat.fill(-999.)

    o_sf_nm = np.zeros(shape=(len(lhist))) ; o_sf_nm.fill(-999.)
    o_sf_cen= np.zeros(shape=(len(lhist))) ; o_sf_cen.fill(-999.)
    o_sf_sat= np.zeros(shape=(len(lhist))) ; o_sf_sat.fill(-999.)

    o_nm = np.zeros(shape=(len(surveys),len(lhist))) ; o_nm.fill(-999.)
    o_cen= np.zeros(shape=(len(surveys),len(lhist))) ; o_cen.fill(-999.)
    o_sat= np.zeros(shape=(len(surveys),len(lhist))) ; o_sat.fill(-999.)

    volume = 0. ; firstpass = True
    for ivol in range(nvol):
        gfile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        #print gfile
        if (os.path.isfile(gfile)):
            # Get some of the model constants
            f = h5py.File(gfile,'r')
            group = f['Parameters']
            vol1 = group['volume'].value ; volume = volume + vol1
            h0 = group['h0'].value 
            omega0 =  group['omega0'].value
            omegab =  group['omegab'].value
            lambda0 = group['lambda0'].value
            zz   = f['Output001/redshift'].value
                                                                   
            set_cosmology(omega0=omega0,omegab=omegab, \
                              lambda0=lambda0, h0=h0, \
                              universe="Flat",include_radiation=False)

            slim = 1./tHubble(zz)

            tomag = band_corrected_distance_modulus(zz)
            if (firstpass):
                szz = "{:.2f}".format(zz)     

            efile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r')
                # Number of haloes per bin
                group = f['Output001']
                gtype  = group['type'].value
                lum_ext =group['L_tot_'+line+'_ext'].value
                BoT = group['BoT'].value
                ind = np.where(gtype == 0)
                mass = group['mstars_tot'].value
                sfr = group['mstardot'].value +\
                      group['mstardot_burst'].value # Msolar/h/Gyr

                ssfr = np.zeros(len(mass)) ; ssfr.fill(-999.)
                ind = np.where(mass>0.)
                ssfr[ind] = sfr[ind]/mass[ind]

                # All
                ind = np.where((mass>0.))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mass[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    all_nm = all_nm + H
                    
                # Satellites 
                ind = np.where((mass>0.) & (gtype>0))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mass[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    all_sat = all_sat + H
                    
                # Centrals 
                ind = np.where((mass>0.) & (gtype==0)) 
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mass[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    all_cen = all_cen + H

                # All SF
                ind = np.where((mass>0.) & (ssfr>=0.3*slim))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mass[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    sf_nm = sf_nm + H
                    
                # Satellites SF 
                ind = np.where((mass>0.) & (ssfr>=0.3*slim) & (gtype>0))
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mass[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    sf_sat = sf_sat + H
                    
                # Centrals SF
                ind = np.where((mass>0.) & (ssfr>=0.3*slim) & (gtype==0)) 
                if (np.shape(ind)[1] > 0):
                    ll = np.log10(mass[ind])
                    H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                    sf_cen = sf_cen + H

                for index,survey in enumerate(surveys):
                    if (survey == 'DEEP2'):
                        fluxcut = 2.7*10.**-17
                        mcut = 24.1
                        band = 'DEIMOS-R'
                        
                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        sel0 = (mag < mcut)
                        
                    elif (survey == 'VVDS-DEEP'):
                        fluxcut = 1.9*10.**-17.
                        mcut = 24.
                        band = 'MegaCam-i-atmos'
                        
                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        sel0 = (mag < mcut)
                        
                    elif (survey == 'VVDS-WIDE'):
                        fluxcut = 3.5*10.**-17.
                        mcut = 22.5
                        band = 'MegaCam-i-atmos'
                        
                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        sel0 = (mag < mcut)
                        
                    elif (survey == 'eBOSS'):
                        fluxcut = 10.**-16. #erg/s/cm^2
                        
                        g = group['mag_DES-g_o_tot_ext'].value + tomag 
                        r = group['mag_DES-r_o_tot_ext'].value + tomag 
                        z = group['mag_DES-z_o_tot_ext'].value + tomag 
                        rz = r-z ; gr = g-r
                        
                        sel0 = (g>21.825) & (g<22.825) & \
                            (gr>-0.068*rz + 0.457) & \
                            (gr< 0.112*rz + 0.773) & \
                            (rz> 0.218*gr + 0.571) & \
                            (rz<-0.555*gr + 1.901)
                        
                    elif (survey == 'DESI'):
                        fluxcut = 8.*10.**-17. #erg/s/cm^2
                        
                        g = group['mag_DES-g_o_tot_ext'].value + tomag 
                        r = group['mag_DES-r_o_tot_ext'].value + tomag 
                        z = group['mag_DES-z_o_tot_ext'].value + tomag 
                        rz = r-z ; gr = g-r
                        
                        sel0 = (r<23.4) & (rz>0.3) & (gr>-0.3) & \
                            (rz>0.9*gr+0.12) & (rz<1.345-0.85*gr)

                    lcut = emission_line_luminosity(fluxcut,zz)

                    # All
                    sel = sel0 & (lum_ext>lcut) & (mass>0.)
                    ind = np.where(sel)
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mass[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        nm[index,:] = nm[index,:] + H
                        
                    # Satellites 
                    sel = sel0 & (lum_ext>lcut)  & (mass>0.) & (gtype>0)
                    ind = np.where(sel)
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mass[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        sat[index,:] = sat[index,:] + H
                        
                    # Centrals 
                    sel = sel0 & (lum_ext>lcut)  & (mass>0.) & (gtype==0)
                    ind = np.where(sel) ; print np.shape(ind)[1],' Centrals'
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mass[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cen[index,:] = cen[index,:] + H

            if firstpass:
                firstpass = False
            f.close()
                
    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)
    if (volume>0.):
        ind = np.where(all_nm>0.)
        o_all_nm[ind] = np.log10(all_nm[ind]/dl/volume)

        ind = np.where(all_cen>0.)
        o_all_cen[ind] = np.log10(all_cen[ind]/dl/volume)

        ind = np.where(all_sat>0.)
        o_all_sat[ind] = np.log10(all_sat[ind]/dl/volume)

        ind = np.where(sf_nm>0.)
        o_sf_nm[ind] = np.log10(sf_nm[ind]/dl/volume)

        ind = np.where(sf_cen>0.)
        o_sf_cen[ind] = np.log10(sf_cen[ind]/dl/volume)

        ind = np.where(sf_sat>0.)
        o_sf_sat[ind] = np.log10(sf_sat[ind]/dl/volume)

        ind = np.where(nm>0.)
        o_nm[ind] = np.log10(nm[ind]/dl/volume)

        ind = np.where(cen>0.)
        o_cen[ind] = np.log10(cen[ind]/dl/volume)

        ind = np.where(sat>0.)
        o_sat[ind] = np.log10(sat[ind]/dl/volume)        

    # Plot
    ax = fig.add_subplot(jj)
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    ax.set_xlabel(xtit,fontsize = fs) ; ax.set_ylabel(ytit,fontsize = fs)
    ax.text(xmax-(xmax-xmin)*0.3, ymax-(ymax-ymin)*0.07, 'z='+szz)

    py = o_all_nm ; ind = np.where(py>-999.)
    y = py[ind] ; x = lhist[ind]
    ax.plot(x,y,color='k',linestyle='-',linewidth=2.5, label="All")

    py = o_all_sat ; ind = np.where(py>-999.)
    y = py[ind] ; x = lhist[ind]
    ax.plot(x,y,color='k',linestyle=':')

    py = o_all_cen ; ind = np.where(py>-999.)
    y = py[ind] ; x = lhist[ind]
    ax.plot(x,y,color='k',linestyle='--')

    # SF
    py = o_sf_nm ; ind = np.where(py>-999.)
    y = py[ind] ; x = lhist[ind]
    ax.plot(x,y,color='gray',linestyle='-',linewidth=2.5, label="SF")

    py = o_sf_sat ; ind = np.where(py>-999.)
    y = py[ind] ; x = lhist[ind]
    ax.plot(x,y,color='gray',linestyle=':')

    py = o_sf_cen ; ind = np.where(py>-999.)
    y = py[ind] ; x = lhist[ind]
    ax.plot(x,y,color='gray',linestyle='--')

    # Plot the model predictions
    for index,survey in enumerate(surveys):
        py = o_nm[index,:] ; ind = np.where(py>-999.)
        y = py[ind] ; x = lhist[ind]
        ax.plot(x,y,color=cols[index],linestyle='-',
                linewidth=2.5, label=survey)

        py = o_sat[index,:] ; ind = np.where(py>-999.)
        y = py[ind] ; x = lhist[ind]
        ax.plot(x,y,color=cols[index],linestyle=':')

        py = o_cen[index,:] ; ind = np.where(py>-999.)
        y = py[ind] ; x = lhist[ind]
        ax.plot(x,y,color=cols[index],linestyle='--')

        # Legend
        leg = plt.legend(loc=3,prop={'size':10})
        leg.draw_frame(False)
                
        # Output files
        tofile = zip(lhist, o_all_nm, o_all_sat, o_all_cen, 
                     o_sf_nm, o_sf_sat, o_sf_cen, 
                     o_nm[index,:],o_sat[index,:],o_cen[index,:]
                     )
        outfil = outdir+'selections/'+model+'/gsmf/gsmf_z'+str(szz)+'_'+survey+'.txt'

        with open(outfil, 'w') as outf:
            outf.write('# GSMF: '+model+', snap='+str(zsnap)+
                       ', '+line+', '+survey+' \n')
            extras = "# h0= %.2f, omega0= %.2f, omegab= %.2f, lambda0= %.2f \n" % (h0, omega0, omegab,lambda0)
            outf.write(extras)
            extras = "# Simulation box side (Mpc/h)= %.2f \n" % (pow(volume,1./3.))
            outf.write(extras)
            outf.write('# \n')
            outf.write('# log10(M*/Msunh-1, All total, All satellites, All centrals, SF total, SF satellites, SF centrals, Survey total, Survey satellites, Survey centrals) \n')
            np.savetxt(outf,tofile,fmt=('%.5f'), delimiter=',')
        outf.closed


# Save figures
fig.tight_layout()
fig.savefig(plotfile)
print 'Output: ',plotfile
print '        ',outfil
