# Output the x,y,z coordinates of galaxies selected in OII flux into 
# an ascii file, either in redshift or real-space

import os.path, sys    #! /usr/bin/env python
import h5py
import numpy as np
from Cosmology import *

Testing = False

propname = 'elg'

nvol = 64

sn_list = ['39','41']
surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']

if Testing:
    nvol = 2
    sn_list = ['39']

#############################
path = '/cosma6/data/dp004/dc-gonz3/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp19/'

line = 'OII3727' ; lline = '[OII]'
############################################
# Path to mass cuts
ndpath = '/cosma6/data/dp004/dc-gonz3/lines/cosmicweb/selections/'

# Generate output files with a header
for iz, sn in enumerate(sn_list):
    surveys = [surveys1[iz],surveys2[iz]]

    for survey in surveys:
        outm = ndpath+model+'ascii_files/'+propname+\
               'cut_'+survey+'_sn'+sn+'.dat'
        print('Output: {}'.format(outm))
        outf = open(outm, 'w')
        outf.write('# xgal,ygal,zgal (Mpc/h), vxgal,vygal,vzgal (Km/s), log10(massh),log10(mass/Msun/h), log10(sfr/Msun/h/Gyr), log10(lum/h^-2 erg/s) ,log10(lum_ext/h^-2 erg/s), type (0= Centrals; 1,2= Satellites) \n')
        outf.close()

for iz, sn in enumerate(sn_list):
    surveys = [surveys1[iz],surveys2[iz]]

    volume = 0.
    for ivol in range(nvol):
        gfile = path+model+'iz'+sn+'/ivol'+str(ivol)+'/galaxies.hdf5'
        if (os.path.isfile(gfile)):        
            # Get some of the model constants
            f = h5py.File(gfile,'r')
            group = f['Parameters']
            vol1 = group['volume'][()] ; volume = volume + vol1
            h0 = group['h0'][()] ; lambda0 =group['lambda0'][()]
            omega0 = group['omega0'][()] ; omegab = group['omegab'][()]
    
            group = f['Output001']
            zz     = group['redshift'][()]
            set_cosmology(omega0=omega0,omegab=omegab,\
                              lambda0=lambda0,h0=h0,\
                              universe="Flat",\
                              include_radiation=False)
            tomag = band_corrected_distance_modulus(zz)
    
            xgal   = group['xgal'][:]   # Mpc/h
            ygal   = group['ygal'][:]
            zgal   = group['zgal'][:]
            vxgal  = group['vxgal'][:]*(1.+zz)/H(zz)  # km/s
            vygal  = group['vygal'][:]*(1.+zz)/H(zz)
            vzgal  = group['vzgal'][:]*(1.+zz)/H(zz)
    
            mhhalo = group['mhhalo'][:]   # Msun/h
            gtype  = group['type'][:] # 0= Centrals; 1,2= Satellites
    
            mdisk = group['mstars_disk'][:] # Msun/h
            mbulge = group['mstars_bulge'][:]
            mass1 = mdisk + mbulge
    
            sdisk = group['mstardot'][:] # Msolar/h/Gyr
            sbulge = group['mstardot_burst'][:]
            sfr1 = sdisk + sbulge
    
            f.close()
    
            gfile = path+model+'iz'+sn+'/ivol'+str(ivol)+'/elgs.hdf5'
            if (not os.path.isfile(gfile)):        
                print('STOP {} not found'.format(gfile)) ; sys.exit()
            f = h5py.File(gfile,'r')
            group = f['Output001']
    
            # Get the log10 of the luminosities (10^40 h^-2 erg/s)
            lum = group['L_tot_'+line][:] 
            llum = np.zeros(shape=(len(lum))) ; llum.fill(-999.)
            ind = np.where(lum>0.)
            llum[ind] = np.log10(lum[ind]) +40.
    
            lum_ext = group['L_tot_'+line+'_ext'][:] 
            llum_ext = np.zeros(shape=(len(lum_ext))) ; llum_ext.fill(-999.)
            ind = np.where(lum_ext>0.)
            llum_ext[ind] = np.log10(lum_ext[ind]) +40.
    
            for survey in surveys:
                if (survey == 'DEEP2'):
                    fluxcut = 2.7*10.**-17
                    mcut = 24.1
                    band = 'DEIMOS-R'
                    
                    mag = group['mag_'+band+'_o_tot_ext'][:] + tomag
                    sel0 = (mag < mcut)
                        
                elif (survey == 'VVDS-DEEP'):
                    fluxcut = 1.9*10.**-17.
                    mcut = 24.
                    band = 'MegaCam-i-atmos'
                    
                    mag = group['mag_'+band+'_o_tot_ext'][:] + tomag
                    sel0 = (mag <= mcut)
                    
                elif (survey == 'VVDS-WIDE'):
                    fluxcut = 3.5*10.**-17.
                    mcut = 22.5
                    band = 'MegaCam-i-atmos'
                        
                    mag = group['mag_'+band+'_o_tot_ext'][:] + tomag
                    sel0 = (mag <= mcut)
                        
                elif (survey == 'eBOSS-SGC'):
                    fluxcut = 10.**-16. #erg/s/cm^2
                    
                    g = group['mag_DES-g_o_tot_ext'][:] + tomag 
                    r = group['mag_DES-r_o_tot_ext'][:] + tomag 
                    z = group['mag_DES-z_o_tot_ext'][:] + tomag 
                    rz = r-z ; gr = g-r
                    
                    sel0 = (g>21.825) & (g<22.825) & \
                           (gr>-0.068*rz + 0.457) & \
                           (gr< 0.112*rz + 0.773) & \
                           (rz> 0.218*gr + 0.571) & \
                           (rz<-0.555*gr + 1.901)
                        
                elif (survey == 'DESI'):
                    fluxcut = 8.*10.**-17. #erg/s/cm^2
                        
                    g = group['mag_DES-g_o_tot_ext'][:] + tomag 
                    r = group['mag_DES-r_o_tot_ext'][:] + tomag 
                    z = group['mag_DES-z_o_tot_ext'][:] + tomag 
                    rz = r-z ; gr = g-r
                    
                    sel0 = (r<23.4) & (rz>0.3) & (gr>-0.3) & \
                           (gr<1.1*rz-0.13) & (gr<1.6-1.18*rz)
    
                cut = emission_line_luminosity(fluxcut,zz)
                if cut>0.: 
                    lcut = np.log10(cut) + 40.
                else:
                    print('STOP: L_cut={} <0 ???'.format(cut))
                    sys.exit()
    
                ind = np.where((mhhalo>0.) & (mass1>0.) & (sfr1>0.) &
                               sel0 & (llum_ext>lcut))
                if (np.shape(ind)[1]<1): continue
    
                massh = np.log10(mhhalo[ind])
                mass = np.log10(mass1[ind])
                sfr = np.log10(sfr1[ind])
    
                tofile = np.column_stack((xgal[ind],\
                                          ygal[ind],\
                                          zgal[ind],\
                                          vxgal[ind],\
                                          vygal[ind],\
                                          vzgal[ind],\
                                          massh,mass,sfr,\
                                          llum[ind],llum_ext[ind], gtype[ind]))
    
                outm = ndpath+model+'ascii_files/'+propname+\
                       'cut_'+survey+'_sn'+sn+'.dat'
                with open(outm,'a') as outf:
                    np.savetxt(outf, tofile, fmt ='%.5f')
    
            f.close()

lbox = pow(volume,1./3.)
print('{}, Box side (Mpc/h) ={}'.format(zz,lbox))
