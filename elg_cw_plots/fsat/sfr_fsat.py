#! /usr/bin/env python
import numpy as np
import os.path, sys
import h5py
from Cosmology import *

# Get fraction of satellites as a function of a cut in L[OII]

model = 'gp19'
#model = 'gp19.font'
#model = 'gp19.starvation'

snap_list = [41, 39] #MillGas

nvol = 64

surveys = ['All','logL>39.5','eBOSS-SGC','VVDS-DEEP','DEEP2','DESI']

#################
# Set output file for the property.
prop = 'sfr'
info_prop = '# fsat as a function of stellar mass log(SFR/h^-1 Msun Gyr^-1) \n'

# Generate the x-axis with the corresponding property
pmin = 3.
pmax = 12.
dp = 0.1
pbins = np.arange(pmin,pmax,dp)
phist = pbins + dp*0.5

#################
line = 'OII3727' ; lline = '[OII]'
path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
#################
# Loop over redshifts, volumes and surveys
outpath = '/cosma5/data/durham/violeta/lines/cosmicweb/fsat/'+model+'/'

for sn in snap_list:
    outfile = outpath+prop+'_'+str(sn)+'.dat'
    with open(outfile, 'w' ) as outf:
        outf.write(info_prop)

    # Initialize the arrays to calculate the fraction of satellites
    fsat = np.zeros(shape=(len(pbins),len(surveys))) ; fsat.fill(-999.)
    ntot = np.zeros(shape=(len(pbins),len(surveys)))
    nsat = np.copy(ntot)

    volume = 0. ; firstpass = True 
    for ivol in range(nvol):
        gfile = path+model+'/iz'+str(sn)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if(os.path.isfile(gfile)):
            f = h5py.File(gfile,'r')
            group = f['Parameters']
            # Get this file's volume
            vol1 = group['volume'].value ; volume = volume + vol1
            # Set the cosmology
            h0 = group['h0'].value
            omega0 = group['omega0'].value
            omegab = group['omegab'].value
            lambda0 = group['lambda0'].value
            set_cosmology(omega0=omega0, omegab=omegab,
                          lambda0=lambda0, h0=h0,
                          universe="Flat", include_radiation=False)

            # Get galaxies' properties
            group = f['Output001']
            zz = group['redshift'].value
            gtype  = group['type'].value
            sfr = group['mstardot'].value + group['mstardot_burst'].value
            #mass = group['mstars_disk'].value + group['mstars_bulge'].value           
            f.close()

            # logSFR
            lsfr = np.zeros(shape=(len(sfr))) ; lsfr.fill(-999.)
            ind = np.where(sfr>0.)
            lsfr[ind] = np.log10(sfr[ind])

            # Get modulus to convert magnitudes to observed ones
            tomag = band_corrected_distance_modulus(zz) # To convert to obs. mags.

            ## sSFR
            #slim = 0.3/tHubble(zz) # sSFR boundary from Franx+08
            #lssfr = np.zeros(shape=(len(mass))) ; lssfr.fill(-999.)
            #ind = np.where((sfr>0.) & (mass>0.))
            #lssfr[ind] = np.log10(sfr[ind]) - np.log10(mass[ind])

            if (firstpass): #Store zz as a string for file names
                szz = "{:.2f}".format(zz)

            efile = path+model+'/iz'+str(sn)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if(os.path.isfile(efile)):
                # Get elgs' information
                f = h5py.File(efile,'r') ; group = f['Output001']
                lum_ext = group['L_tot_'+line+'_ext'].value
                # log(L/h^-2 erg s-1)
                logl = np.zeros(shape=len(lum_ext)) ; logl.fill(-999.)
                ind = np.where(lum_ext>0.)
                logl[ind] = np.log10(lum_ext[ind]) + 40.

                for isy,survey in enumerate(surveys):
                    if (survey == 'logL>39.5'):
                        sel0 = (logl > 39.5)

                    elif (survey == 'DEEP2'):
                        fluxcut = 2.7*10.**-17
                        mcut = 24.1
                        band = 'DEIMOS-R'
                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        lcut = emission_line_luminosity(fluxcut,zz)
                        sel0 = (mag < mcut) & (lum_ext > lcut)

                    elif (survey == 'VVDS-DEEP'):
                        fluxcut = 1.9*10.**-17.
                        mcut = 24.
                        band = 'MegaCam-i-atmos'
                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        lcut = emission_line_luminosity(fluxcut,zz)
                        sel0 = (mag < mcut) & (lum_ext > lcut)

                    elif (survey == 'VVDS-WIDE'):
                        fluxcut = 3.5*10.**-17.
                        mcut = 22.5
                        band = 'MegaCam-i-atmos'
                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        lcut = emission_line_luminosity(fluxcut,zz)
                        sel0 = (mag < mcut) & (lum_ext > lcut)

                    elif (survey == 'eBOSS-SGC'):
                        fluxcut = 10.**-16. #erg/s/cm^2

                        g = group['mag_DES-g_o_tot_ext'].value + tomag 
                        r = group['mag_DES-r_o_tot_ext'].value + tomag 
                        z = group['mag_DES-z_o_tot_ext'].value + tomag 
                        rz = r-z ; gr = g-r

                        lcut = emission_line_luminosity(fluxcut,zz)
                        sel0 = (g>21.825) & (g<22.825) & \
                            (gr>-0.068*rz + 0.457) & \
                            (gr< 0.112*rz + 0.773) & \
                            (rz> 0.218*gr + 0.571) & \
                            (rz<-0.555*gr + 1.901) & \
                            (lum_ext > lcut)

                    elif (survey == 'DESI'):
                        fluxcut = 8.*10.**-17. #erg/s/cm^2

                        g = group['mag_DES-g_o_tot_ext'].value + tomag 
                        r = group['mag_DES-r_o_tot_ext'].value + tomag 
                        z = group['mag_DES-z_o_tot_ext'].value + tomag 
                        rz = r-z ; gr = g-r

                        lcut = emission_line_luminosity(fluxcut,zz)
                        sel0 = (r<23.4) & (rz>0.3) & (gr>-0.3) & \
                            (gr<1.1*rz-0.13) & (gr<1.6-1.18*rz) & \
                            (lum_ext > lcut)

                    # Loop over the properties limits
                    for ip,pbin in enumerate(pbins):
                        if (survey == 'All'):
                            sel = (lsfr>pbin)
                        else:
                            sel = sel0 & (lsfr>pbin)

                        ind = np.where(sel) # All
                        if(np.shape(ind)[1] > 1):
                            ntot[ip,isy] = ntot[ip,isy] + np.shape(ind)[1]

                        ind = np.where(sel & (gtype >0)) # Satellites
                        if(np.shape(ind)[1] > 1):
                            nsat[ip,isy] = nsat[ip,isy] + np.shape(ind)[1]

            f.close()
            if (firstpass):
                firstpass = False

    # Get the percentage of satellites
    ind = np.where(ntot > 0.)
    if(np.shape(ind)[1] > 1):
        fsat[ind] = 100.*nsat[ind]/ntot[ind]

    # Write output
    with open(outfile, 'a' ) as outf:
        outf.write('# Volume = {} \n'.format(volume))
        outf.write('# redshift = {}, h0 = {}, omega0 = {}, omegab = {}, lambda0 = {} \n'.
                   format(szz,h0,omega0,omegab,lambda0))
        outf.write('# logSFR_min %satellites('+' '.join(surveys)+') \n') 

        tofile = np.column_stack([pbins,fsat])
        np.savetxt(outf,tofile, fmt='%.5f')

    print('Output: {}'.format(outfile))
