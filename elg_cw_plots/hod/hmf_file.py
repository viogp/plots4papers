import numpy as np
import os.path, sys
import h5py
from Cosmology import * 

path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'

snap_list = [39, 41]
nvol = 64

model = 'gp19'

#############################
line = 'OII3727' ; lline = '[OII]'
outdir = '/cosma5/data/durham/violeta/lines/cosmicweb/'
#############################

# Initialize histogram
lmin = 8.5
lmax = 16.
dl = 0.1
lbins = np.arange(lmin,lmax,dl)
lhist = lbins + dl*0.5

allbins = np.append(lbins,lmax)
allb_low = allbins[:-1] 
allb_high = allbins[1:]

############################################

# Loop over the redshifts of interest
for iz,zsnap in enumerate(snap_list):
    # Initialize output file
    outfil = outdir+'hod/'+model+'/hmf_sn'+str(zsnap)+'.txt'
    #outfil = outdir+'hod/'+model+'/hmf_sn'+str(zsnap)+'_2vols.txt'
    with open(outfil, 'w') as outf:
        outf.write('# HMF: '+model+', snap='+str(zsnap)+' \n')

    # Initialize array with number of haloes
    nh = np.zeros(shape=(len(lhist)))

    volume = 0. 
    for ivol in range(nvol):
        gfile = path+model+'/iz'+str(zsnap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
        # Check the existance of the file
        if (not os.path.isfile(gfile)):
            print('WARNING: File not found {}'.format(gfile))
            continue

        # Get some of the model constants
        f = h5py.File(gfile,'r')
        group = f['Parameters']
        vol1 = group['volume'].value ; volume = volume + vol1
        h0 = group['h0'].value 
        omega0 =  group['omega0'].value
        omegab =  group['omegab'].value
        lambda0 = group['lambda0'].value
        
        zz   = f['Output001/redshift'].value
        mhhalo = f['Output001/mhhalo'].value
        gtype  = f['Output001/type'].value

        ind = np.where(gtype == 0)
        ll = np.log10(mhhalo[ind])
        H, bins_edges = np.histogram(ll,bins=allbins)
        nh = nh + H
                
    # Output 
    indh = np.where(nh > 0) 
    tofile = zip(lhist[indh], allb_low[indh], allb_high[indh], nh[indh])
    with open(outfil, 'a') as outf:
        extras = "# h0= %.2f, omega0= %.2f, omegab= %.2f, lambda0= %.2f \n" % (h0, omega0, omegab,lambda0)
        outf.write(extras)
        extras = "# Simulation box side (Mpc/h)= %.2f \n" % (pow(volume,1./3.))
        outf.write(extras)
        outf.write('# \n')
        extras = "# dl %.2f \n" % (dl)
        outf.write(extras)
        outf.write('# log(Mh/Msusn/h)_midpoint, log(Mh)_low, log(Mh)_high, Number of haloes \n')
        
        np.savetxt(outf,tofile,fmt=('%.5f'))
    outf.closed

print('Output: {}'.format(outfil))
