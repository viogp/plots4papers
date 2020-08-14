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

#################
# Set output file for the property.
prop = 'g_fixedlim'
info_prop = '# fsat as a function of g-band with fixed limits \n'

# Fixed limits from Comparat+15
gmin = [21.,20.,20.,20.5,22.1]
gmax = [22.5,23.,22.5,22.8,22.8]

llim = np.arange(39.5,42.,0.5)

surveys = ['All']
for ll in llim:
    surveys.append('logL>'+str(ll))

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
    fsat = np.zeros(shape=(len(gmin),len(surveys))) ; fsat.fill(-999.)
    ntot = np.zeros(shape=(len(gmin),len(surveys)))
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
            f.close()

            # Get modulus to convert magnitudes to observed ones
            tomag = band_corrected_distance_modulus(zz) # To convert to obs. mags.

            if (firstpass): #Store zz as a string for file names
                szz = "{:.2f}".format(zz)

            efile = path+model+'/iz'+str(sn)+'/ivol'+str(ivol)+'/elgs.hdf5'
            if(os.path.isfile(efile)):
                # Get elgs' information
                f = h5py.File(efile,'r') ; group = f['Output001']
                g = group['mag_DES-g_o_tot_ext'].value + tomag 
                lum_ext = group['L_tot_'+line+'_ext'].value
                # log(L/h^-2 erg s-1)
                logl = np.zeros(shape=len(lum_ext)) ; logl.fill(-999.)
                ind = np.where(lum_ext>0.)
                logl[ind] = np.log10(lum_ext[ind]) + 40.

                for isy,survey in enumerate(surveys):
                    if (isy > 0):
                        sel0 = (logl > llim[isy-1])

                    # Loop over the properties limits
                    for ig in range(len(gmin)):
                        if (survey == 'All'):
                            sel = (g > gmin[ig]) & (g < gmax[ig])
                        else:
                            sel = sel0 & (g > gmin[ig]) & (g < gmax[ig])
                        
                        ind = np.where(sel) # All
                        if(np.shape(ind)[1] > 1):
                            ntot[ig,isy] = ntot[ig,isy] + np.shape(ind)[1]
                        
                        ind = np.where(sel & (gtype >0)) # Satellites
                        if(np.shape(ind)[1] > 1):
                            nsat[ig,isy] = nsat[ig,isy] + np.shape(ind)[1]

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
        outf.write('# gmin gmax %satellites('+' '.join(surveys)+') \n') 

        tofile = np.column_stack([gmin,gmax,fsat])
        np.savetxt(outf,tofile, fmt='%.5f')

    print('Output: {}'.format(outfile))
