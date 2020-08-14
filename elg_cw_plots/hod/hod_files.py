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

snap_list = [41]#[44, 42, 41, 40, 39, 37, 34] #MillGas
nvol = 2#64

model = 'gp18' #'gp18.font','gp18.starvation','gp17'

surveys = ['eBOSS']#['DEEP2','VVDS-WIDE','VVDS-DEEP','eBOSS','DESI']

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/gpfs/data/violeta/lines/cosmicweb/'
plotfile = outdir+'plots/'+model+'/hod/hod_'+line+'_all.pdf'
#############################

ntypes = len(surveys)
cols = get_distinct(ntypes) 

# Initialize histogram
lmin = 8.5
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

############################################
# Initialize the parameters for the figures
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15
fig = plt.figure(figsize=(9.,15.))

xtit = "${\\rm log}_{10}(M_{\\rm halo}/M_{\odot}h^{-1})$"
ytit = "$\\langle N_M\\rangle$"

xmin = 10. ; xmax = 15.
ymin = -3. ; ymax = 2.

# Loop over the redshifts of interest
jj = 420
for iz,zsnap in enumerate(snap_list):
    jj = jj + 1

    nm = np.zeros(shape=(len(surveys),len(lhist)))
    cen= np.zeros(shape=(len(surveys),len(lhist)))
    sat= np.zeros(shape=(len(surveys),len(lhist)))
    cs = np.zeros(shape=(len(surveys),len(lhist)))
    cd = np.zeros(shape=(len(surveys),len(lhist)))
    nh = np.zeros(shape=(len(lhist)))
    mhcd = np.zeros(shape=(len(surveys)))
    mhcs = np.zeros(shape=(len(surveys)))

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
            sfr = f['Output001/mstardot'].value +\
                f['Output001/mstardot_burst'].value
            mass = f['Output001/mstars_disk'].value +\
                f['Output001/mstars_bulge'].value
            lssfr = np.zeros(shape=(len(mass))) ; lssfr.fill(-999.)
            ind = np.where((sfr>0.) & (mass>0.))                   
            lssfr[ind] = np.log10(sfr[ind]) - np.log10(mass[ind])
                                                                   
            set_cosmology(omega0=omega0,omegab=omegab, \
                              lambda0=lambda0, h0=h0, \
                              universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz)
            slim = 0.3/tHubble(zz) # Franx+08
            if (firstpass):
                szz = "{:.2f}".format(zz)     

            efile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r')
                # Number of haloes per bin
                group = f['Output001']
                gtype  = group['type'].value
                mhhalo = group['mhhalo'].value
                lum_ext =group['L_tot_'+line+'_ext'].value
                BoT = group['BoT'].value
                ind = np.where(gtype == 0)
                ll = np.log10(mhhalo[ind])
                H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                nh = nh + H

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
                        sel0 = (mag <= mcut)
                        
                    elif (survey == 'VVDS-WIDE'):
                        fluxcut = 3.5*10.**-17.
                        mcut = 22.5
                        band = 'MegaCam-i-atmos'
                        
                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        sel0 = (mag <= mcut)
                        
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
                               (gr<1.1*rz-0.13) & (gr<1.6-1.18*rz)

                    lcut = emission_line_luminosity(fluxcut,zz)

                    # All
                    sel = sel0 & (lum_ext>lcut)
                    ind = np.where(sel)
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        nm[index,:] = nm[index,:] + H
                        
                    # Satellites 
                    sel = sel0 & (lum_ext>lcut) & (gtype>0)
                    ind = np.where(sel)
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        sat[index,:] = sat[index,:] + H
                        
                    # Centrals 
                    sel = sel0 & (lum_ext>lcut) & (gtype==0)
                    ind = np.where(sel) ; print np.shape(ind)[1],' Centrals'
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cen[index,:] = cen[index,:] + H

                    # Centrals Spheroids
                    sel = sel0 & (lum_ext>lcut) & (gtype==0) & (BoT>0.5)
                    ind = np.where(sel) ; print np.shape(ind)[1],' C-Sph'
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cs[index,:] = cs[index,:] + H

                        av = np.mean(mhhalo[ind])
                        if firstpass:
                            mhcs[index] = av
                        else:
                            mhcs[index] = (mhcs[index]+av)/2.

                    # Centrals Disk
                    sel = sel0 & (lum_ext>lcut) & (gtype==0) & (BoT<=0.5)
                    ind = np.where(sel) ; print np.shape(ind)[1],' C-Disk'
                    if (np.shape(ind)[1] > 0):
                        ll = np.log10(mhhalo[ind])
                        H, bins_edges = np.histogram(ll,bins=np.append(lbins,lmax))
                        cd[index,:] = cd[index,:] + H

                        av = np.mean(mhhalo[ind])
                        if firstpass:
                            mhcd[index] = av
                        else:
                            mhcd[index] = (mhcd[index]+av)/2.

            if firstpass:
                firstpass = False
            f.close()
                
    print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

    # Plot
    ax = fig.add_subplot(jj)
    ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
    ax.set_xlabel(xtit,fontsize = fs) ; ax.set_ylabel(ytit,fontsize = fs)
    ax.text(xmax-(xmax-xmin)*0.3, ymin+(ymax-ymin)*0.05, 'z='+szz)

    # Plot the model predictions
    for index,survey in enumerate(surveys):
        ig = survey
        indh = np.where(nh > 0) ; nhalos = sum(nh)
        x = lhist[indh]

        py = nm[index,:] ; nall = sum(py)
        yall = py[indh]/nh[indh] 
        y = np.log10(yall) 
        ax.plot(x,y,color=cols[index],linestyle='-',linewidth=2.5)

        py = sat[index,:] ; nsat = sum(py)
        ysat = py[indh]/nh[indh]
        y = np.log10(ysat) 
        ax.plot(x,y,color=cols[index],linestyle=':')

        py = cen[index,:] ; ncen = sum(py)
        ycen = py[indh]/nh[indh] 
        y = np.log10(ycen) 
        ax.plot(x,y,color=cols[index],linestyle='-')

        py = cs[index,:] ; ncs = sum(py)
        ycs = py[indh]/nh[indh]
        y = np.log10(ycs) 
        ax.plot(x,y,color=cols[index],linestyle='-.')

        py = cd[index,:]  ; ncd = sum(py)
        ycd = py[indh]/nh[indh]
        y = np.log10(ycd)
        if (mhcd[index]>0. and mhcs[index]>0.):
            extras = "<Md>=%.1f, <Ms>=%.1f" % (np.log10(mhcd[index]),np.log10(mhcs[index]))
            ax.plot(x,y,color=cols[index],\
                        linestyle='--',\
                        label=ig+extras)
        else:
            ax.plot(x,y,color=cols[index],\
                        linestyle='--', label=ig)

        # Legend
        leg = plt.legend(loc=2,prop={'size':10})
        for color,text in zip(cols,leg.get_texts()):
            text.set_color(color)
            leg.draw_frame(False)
                
        # Output files
        tofile = zip(x, yall,ysat, ycen, ycs, ycd)
        outfil = outdir+'hod/'+model+'/hod_z'+str(szz)+'_'+survey+'.txt'

        with open(outfil, 'w') as outf:
            outf.write('# Mean HOD: '+model+', snap='+str(zsnap)+', '+line+', '+ig+' \n')
            extras = "# h0= %.2f, omega0= %.2f, omegab= %.2f, lambda0= %.2f \n" % (h0, omega0, omegab,lambda0)
            outf.write(extras)
            extras = "# Simulation box side (Mpc/h)= %.2f \n" % (pow(volume,1./3.))
            outf.write(extras)
            outf.write('# \n')
            outf.write('# Total number of haloes='+str(nhalos)+' , Number of ELGs='+str(nall)+' \n')
            if (nall > 0.):
                extras = "# %.1f percent of ELGs are satellites and %.1f centrals \n" % (100.*nsat/nall, 100.*ncen/nall)
                outf.write(extras)
                extras = "# %.1f percent of central ELGs are spheroids and %.1f percent disks \n" % (100.*ncs/nall, 100.*ncd/nall)
                print ncs,ncd,ncs+ncd,nall,100.*ncs/nall, 100.*ncd/nall
                outf.write(extras)

            outf.write('# \n')
            outf.write('# log10(M*/Msunh-1), Total, Satellites, Centrals, Central spheroids, Central Disks \n')
            np.savetxt(outf,tofile,fmt=('%.5f'), delimiter=',')
        outf.closed


# Save figures
fig.tight_layout()
fig.savefig(plotfile)
print 'Output: ',plotfile
print '        ',outfil
