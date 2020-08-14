#! /usr/bin/env python

import os, sys
import h5py
import numpy as np
from Cosmology import *

path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp19'
#model = 'gp19.font'
#model = 'gp19.starvation'

info_file = model+'_inputs4cute_r.txt' 

info = open( info_file, 'w' )
#############################
line = 'OII3727' ; lline = '[OII]'
#############################

snap_list = [41, 39] #MillGas
nvol = 64

#inleg = ['DEEP2','VVDS-WIDE','VVDS-DEEP','eBOSS','DESI']

inleg = ['DEEP2','eBOSS-SGC','VVDS-DEEP','DESI']

############################################
# Loop over the redshifts of interest
for iz,zsnap in enumerate(snap_list):
    outpath = path+model+'/iz'+str(zsnap)+'/r/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # Delete previous files
    for inlegs in inleg:
        outfile = outpath+line+'_'+inlegs+'_4cute.dat'
        if(os.path.isfile(outfile)):             
            os.remove(outfile)

    volume = 0. ; firstpass = True
    ngals = np.zeros(shape=(len(inleg)), dtype=np.int)
    for ivol in range(nvol):
        gfile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (os.path.isfile(gfile)):
            # Get some of the model constants
            f = h5py.File(gfile,'r')
            group = f['Parameters']
            vol1 = group['volume'].value ; volume = volume + vol1
            h0 = group['h0'].value ; lambda0 =group['lambda0'].value
            omega0 = group['omega0'].value ; omegab = group['omegab'].value

            group = f['Output001']
            zz     = group['redshift'].value
            set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                          h0=h0, universe="Flat",include_radiation=False)
            tomag = band_corrected_distance_modulus(zz) 

            xgal   = group['xgal'].value   # Mpc/h
            ygal   = group['ygal'].value
            zgal   = group['zgal'].value
            vxgal  = group['vxgal'].value*(1.+zz)/H(zz)  # km/s
            vygal  = group['vygal'].value*(1.+zz)/H(zz)
            vzgal  = group['vzgal'].value*(1.+zz)/H(zz)

            ihhalo = group['ihhalo'].value
            tjm    = group['Trees/jm'].value
            tngals = group['Trees/ngals'].value
            jm = np.repeat(tjm,tngals)

            gtype  = group['type'].value

            sfr  = group['mstardot'].value + group['mstardot_burst'].value
            mass = group['mstars_disk'].value + group['mstars_bulge'].value
            lssfr = np.zeros(shape=(len(mass))) ; lssfr.fill(-999.)
            ind = np.where((sfr>0.) & (mass>0.))                   
            lssfr[ind] = np.log10(sfr[ind]) - np.log10(mass[ind])

            tsf    = 3*group['tsfburst'].value  
            tburst = group['tburst'].value
            burst = np.zeros(shape=(len(tsf)))
            ind = np.where(tburst<tsf)
            burst[ind] = 1.

            f.close()
    
            efile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (os.path.isfile(efile)):
                f = h5py.File(efile,'r')
                group = f['Output001']
                lum_ext = group['L_tot_'+line+'_ext'].value
                lum = group['L_tot_'+line].value

                for index,ileg in enumerate(inleg):
                    if (ileg == 'DEEP2'):
                        fluxcut = 2.7*10.**-17
                        mcut = 24.1
                        band = 'DEIMOS-R'

                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        sel0 = (mag < mcut)

                    elif (ileg == 'VVDS-DEEP'):
                        fluxcut = 1.9*10.**-17.
                        mcut = 24.
                        band = 'MegaCam-i-atmos'

                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        sel0 = (mag < mcut)

                    elif (ileg == 'VVDS-WIDE'):
                        fluxcut = 3.5*10.**-17.
                        mcut = 22.5
                        band = 'MegaCam-i-atmos'

                        mag = group['mag_'+band+'_o_tot_ext'].value + tomag
                        sel0 = (mag < mcut)

                    elif (ileg == 'eBOSS-SGC'):
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

                    elif (ileg == 'DESI'):
                        fluxcut = 8.*10.**-17. #erg/s/cm^2

                        g = group['mag_DES-g_o_tot_ext'].value + tomag 
                        r = group['mag_DES-r_o_tot_ext'].value + tomag 
                        z = group['mag_DES-z_o_tot_ext'].value + tomag 
                        rz = r-z ; gr = g-r

                        sel0 = (r<23.4) & (rz>0.3) & (gr>-0.3) & \
                            (gr<1.1*rz-0.13) & (gr<1.6-1.18*rz)

                    lcut = emission_line_luminosity(fluxcut,zz)

                    sel = sel0 & (lum_ext>lcut)
                    #ind = np.where(sel)
                    ind = np.where(lum_ext>lcut)

                    sel = sel0 & (lum>lcut)
                    indi = np.where(sel)

                    outfile = outpath+line+'_'+inleg[index]+'_4cute.dat'
                    if firstpass:
                        slim = 0.3/tHubble(zz) # Franx+08

                    if(np.shape(ind)[1] > 0): #####
                        ngals[index] = ngals[index] + np.size(ind)
                        vols = np.zeros(shape=(np.size(ind))) ; vols.fill(ivol)
                        tofile = np.column_stack((xgal[ind],ygal[ind],zgal[ind],\
                                                  vxgal[ind],vygal[ind],vzgal[ind],\
                                                  lssfr[ind],\
                                                  jm[ind],ihhalo[ind],vols,\
                                                  gtype[ind],burst[ind]))

                        with open(outfile,'a') as f_handle:
                            if(np.shape(tofile)[1] == 12):
                                np.savetxt(f_handle, tofile, \
                                           fmt=('%.5f %.5f %.5f %.5f %.5f %.5f %.5f %d %d %d %d %d'))

                f.close()

            # Continue loop over subvolumes
            firstpass = False

    lbox = pow(volume,1./3.)
    for index,ileg in enumerate(inleg):
        outfile = outpath+line+'_'+ileg+'_4cute.dat'
        info.write(str(outfile+ '\n') )
info.close()
print info_file    

