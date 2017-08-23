#! /usr/bin/env python

import os.path, sys
import h5py
import numpy as np
from Cosmology import *

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
#model = 'MillGas/gp14/'
model = 'MillGas/gp15newmg/'
info_file = 'inputs4cute_gp15newmg.txt' ; info = open( info_file, 'w' )

#############################
line = 'OII3727' ; lline = '[OII]'
#############################

snap_list = [44, 42, 40, 37, 34] #MillGas
nvol = 64

bands = ['RK','m2','m2','eBOSS','DESI']
mcuts = [24.1,22.5,24]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17.,10.**-16.,8.*10.**-17.]

inleg = ['DEEP2','VVDSWIDE','VVDSDEEP','eBOSS','DESI']

############################################
# Loop over the redshifts of interest
for iz,zsnap in enumerate(snap_list):
    outpath = path+model+'iz'+str(zsnap)+'/'

    volume = 0. ; firstpass = True
    ngals = np.zeros(shape=(len(bands)), dtype=np.int)
    for ivol in range(nvol):
        gfile = path+model+'iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
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
                lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                lum = f['Output001/L_tot_'+line].value

                for index,ib in enumerate(bands):
                    if index<3:
                        mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                        icut = mcuts[index] ; fluxcut = fcuts[index]
                        lcut = emission_line_luminosity(fluxcut,zz)
                        ind  = np.where((mag<icut) & (lum_ext>lcut))
                        indi = np.where((mag<icut) & (lum>lcut))
                    else:
                        g = f['Output001/magDgo_tot_ext'].value + tomag
                        r = f['Output001/magDro_tot_ext'].value + tomag
                        z = f['Output001/magDzo_tot_ext'].value + tomag
                        rz = r-z ; gr = g-r

                        if index==3: #eBOSS decam180 selection
                            fluxcut = fcuts[index]
                            lcut = emission_line_luminosity(fluxcut,zz)

                            ind = np.where((g>22.1) & (g<22.8) & \
                                               (gr>0.3) & (gr<0.7) & \
                                               (rz>0.25) & (rz<1.4) & \
                                               (rz>0.5*gr+0.4) & \
                                               (rz<0.5*gr+0.8) & (lum_ext>lcut))
                            indi = np.where((g>22.1) & (g<22.8) & \
                                       (gr>0.3) & (gr<0.7) & \
                                       (rz>0.25) & (rz<1.4) & \
                                       (rz>0.5*gr+0.4) & \
                                       (rz<0.5*gr+0.8) & (lum>lcut))

                        elif index==4: #DESI                            
                            fluxcut = fcuts[index]
                            lcut = emission_line_luminosity(fluxcut,zz)

                            ind = np.where((r<23.4) & \
                                               (rz>0.3) & (gr>-0.3) & \
                                               (rz>0.9*gr+0.12) & \
                                               (rz<1.345-0.85*gr) & \
                                               (lum_ext>lcut))
                            indi = np.where((r<23.4) & \
                                                (rz>0.3) & (gr>-0.3) & \
                                                (rz>0.9*gr+0.12) & \
                                                (rz<1.345-0.85*gr) & \
                                                (lum>lcut))
                        

                    outfile = outpath+line+'_'+inleg[index]+'_4cute.dat'
                    if firstpass:
                        slim = 0.3/tHubble(zz) # Franx+08


                    ngals[index] = ngals[index] + np.size(ind)
                    vols = np.zeros(shape=(np.size(ind))) ; vols.fill(ivol)
                    tofile = np.column_stack((xgal[ind],ygal[ind],zgal[ind],\
                                              vxgal[ind],vygal[ind],vzgal[ind],\
                                              lssfr[ind],\
                                              jm[ind],ihhalo[ind],vols,\
                                              gtype[ind],burst[ind]))

                    with open(outfile,'a') as f_handle:
                        np.savetxt(f_handle, tofile, \
                                   fmt=('%.5f %.5f %.5f %.5f %.5f %.5f %.5f %d %d %d %d %d'))

                f.close()

            # Continue loop over subvolumes
            firstpass = False

    lbox = pow(volume,1./3.)
    for index,ib in enumerate(bands):
        outfile = outpath+line+'_'+inleg[index]+'_4cute.dat'
        info.write(str(outfile+ '\n') )
info.close()
print info_file    

