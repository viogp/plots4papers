import os.path, sys    #! /usr/bin/env python
import h5py
import numpy as np
from Cosmology import *

path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
#model = 'gp19'
#model = 'gp19.font'
model = 'gp19.starvation'

info_file = model+'_inputs4cute_gflux_z.txt' 

info = open( info_file, 'w' )
#############################
line = 'OII3727' ; lline = '[OII]'
#############################

snap_list = [41, 39] #MillGas
nvol = 64

fcuts = [10.**-18., 10.**-17., 5.*10.**-17., 10.**-16., 5*10.**-16.]
gcuts = [22.8, 22.2]
inleg = ['fm18','fm17','f5m17','fm16','f5m16','g22.8','g22.2']

ntypes = len(inleg)
############################################
# Loop over the redshifts of interest
for iz,zsnap in enumerate(snap_list):
    outpath = path+model+'/iz'+str(zsnap)+'/z/'

    volume = 0. ; firstpass = True
    ngals = np.zeros(shape=(ntypes), dtype=np.int)
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
                lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                lum = f['Output001/L_tot_'+line].value

                for index in range(ntypes):
                    if (index<ntypes-2):
                        fluxcut = fcuts[index]
                        lcut = emission_line_luminosity(fluxcut,zz)
                        ind  = np.where(lum_ext>lcut)
                        inlegs = inleg[index]

                        if (firstpass):
                            print inlegs, fluxcut, lcut
                    else:
                        gcut = gcuts[ntypes-1-index]
                        g = f['Output001/mag_DES-g_o_tot_ext'].value + tomag
                        ind = np.where((g>20.) & (g<gcut))
                        inlegs = 'g'+str(gcut)

                        if (firstpass):
                            print inlegs, gcut

                    outfile = outpath+line+'_'+inlegs+'gflux_4cute_z.dat'

                    xzs = xgal[ind] + vxgal[ind]

                    ngals[index] = ngals[index] + np.size(ind)
                    vols = np.zeros(shape=(np.size(ind))) ; vols.fill(ivol)
                    tofile = np.column_stack((xzs,ygal[ind],zgal[ind],\
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
    for index  in range(ntypes):
        outfile = outpath+line+'_'+inleg[index]+'gflux_4cute_z.dat'
        print outfile
        info.write(str(outfile+ '\n') )
info.close()
print(' - Info file:',info_file)
