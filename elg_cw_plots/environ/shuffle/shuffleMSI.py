# Output the x,y,z original and shuffled coordinates
# for the Millenium-II simulation

import os.path, sys    
import h5py
import numpy as np
from Cosmology import *

# Functions

def add2sat(xgal,ygal,zgal,gtype,add=True):
    '''
    Add or remove the position of the central to the satellites
    '''
    xlast = xgal[0]
    ylast = ygal[0]
    zlast = zgal[0]

    for i in range(len(xgal)):
        if (gtype[i] > 0):
            if (add):
                xgal[i] = xgal[i] + xlast 
                ygal[i] = ygal[i] + ylast
                zgal[i] = zgal[i] + zlast
            else:
                xgal[i] = xgal[i] - xlast
                ygal[i] = ygal[i] - ylast
                zgal[i] = zgal[i] - zlast
        else:
            xlast = xgal[i]
            ylast = ygal[i]
            zlast = zgal[i]
        

def correct_periodic(xgal,ygal,zgal,gtype,lbox,halfbox=False):
    '''
    Correct satellites for periodic boundary conditions
    '''
    lbox2 = lbox/2.
    
    for arr in [xgal,ygal,zgal]:
        if halfbox:
            ind = np.where((arr > lbox2) & (gtype > 0))
            arr[ind] = arr[ind] - lbox

            ind = np.where((arr < -lbox2) & (gtype > 0))
            arr[ind] = arr[ind] + lbox
        else:
            ind = np.where((arr >= lbox) & (gtype > 0))
            arr[ind] = arr[ind] - lbox

            ind = np.where((arr < 0) & (gtype > 0))
            arr[ind] = arr[ind] + lbox

def mhhalo2sat(mh,gtype):
    '''
    Assign the new halo mass to the satellite galaxies 
    '''
    mhlast = mh[0]

    for i in range(len(mh)):
        if (gtype[i] > 0):
            mh[i] = mhlast
        else:
            mhlast = mh[i]

############################################

Testing = True

nvol = 64
sn_list = ['39','41','61']

if Testing:
    nvol = 2
    sn_list = ['39']
    lbox_test = 500.

############################################

path = '/cosma6/data/dp004/dc-gonz3/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp19/'

############################################

nbin = 70
mmin, mmax = 9., 16.

for iz, sn in enumerate(sn_list):
    volume = 0. ; iifil = -1
    for ivol in range(nvol):
        gfile = path+model+'iz'+sn+'/ivol'+str(ivol)+'/galaxies.hdf5'

        if (not os.path.isfile(gfile)):
            continue
        iifil += 1
        if Testing: 
            if iifil>2: break

        # Read the halo information
        f = h5py.File(gfile,'r')
        group = f['Output001']
        vol1 = f['Parameters/volume'][()] ; volume = volume + vol1
        tjm = group['Trees/jm'][:]
        tngals = group['Trees/ngals'][:] 
        mtot = group['mstars_disk'][:] + group['mstars_bulge'][:]
        imgal = np.zeros(len(mtot)) ; imgal.fill(-999.) 
        ind = np.where(mtot>0.)  
        imgal[ind] = np.log10(mtot[ind])

        ixgal   = group['xgal'][:]   # Mpc/h
        iygal   = group['ygal'][:]
        izgal   = group['zgal'][:]
        imhhalo = np.log10(group['mhhalo'][:])   # log10(M/Msun/h)
        igtype  = group['type'][:] # 0= Centrals; 1,2= Satellites
        ijm     = np.repeat(tjm,tngals) # tree index 
        iihhalo = group['ihhalo'][:]    # dhalo index within a tree        

        f.close()

        # Array with subvolumes
        ivols = np.zeros(shape=len(ixgal),dtype=int) ; ivols.fill(ivol)

        # Sort by jm, then by ihhalo, then by gtype
        indsort = np.lexsort((igtype,iihhalo,ijm))

        if (iifil==0):
            xgal   = ixgal[indsort]
            ygal   = iygal[indsort]
            zgal   = izgal[indsort]
            mhhalo = imhhalo[indsort]
            gtype  = igtype[indsort]
            jm     = ijm[indsort]
            ihhalo = iihhalo[indsort]
            mgal   = imgal[indsort]
            vols   = ivols
        else:
            xgal   = np.append(xgal,ixgal[indsort])
            ygal   = np.append(ygal,iygal[indsort])
            zgal   = np.append(zgal,izgal[indsort])
            mhhalo = np.append(mhhalo,imhhalo[indsort])
            gtype  = np.append(gtype,igtype[indsort])
            jm     = np.append(jm,ijm[indsort])
            ihhalo = np.append(ihhalo,iihhalo[indsort])
            mgal   = np.append(mgal,imgal[indsort])
            vols   = np.append(vols,ivols)

    lbox = pow(volume,1./3.)
    if Testing: lbox=lbox_test
    print('sn={}, Box side (Mpc/h) ={}'.format(sn,lbox))

    # Initialise arrays to contain shuffled properties
    s_xgal = np.copy(xgal)
    s_ygal = np.copy(ygal)
    s_zgal = np.copy(zgal)
    s_mhhalo = np.copy(mhhalo)
    s_mgal  = np.copy(mgal)

    # Change the position of the satellites to be relative to the halo
    add2sat(s_xgal,s_ygal,s_zgal,gtype,add=False)

    # Correct for periocid boundary conditions (for satellites)
    correct_periodic(s_xgal,s_ygal,s_zgal,gtype,lbox,halfbox=True) 

    # Find the bin each halo mass belongs to
    imass = ((np.copy(s_mhhalo)-mmin)/(mmax-mmin)*nbin).astype(int)

    # Shuffle the centrals in each halo mass bin
    for i in range(nbin):
        # Select all central galaxies in the halo mass bin
        ind = np.where((gtype == 0) & (i == imass))[0]
        if (len(ind) < 1): continue

        # Get shuffled indexes
        ind_shuffle = np.copy(ind)
        np.random.shuffle(ind_shuffle)

        # Reassigned properties
        s_xgal[ind]   = s_xgal[ind_shuffle]
        s_ygal[ind]   = s_ygal[ind_shuffle]
        s_zgal[ind]   = s_zgal[ind_shuffle]
        s_mhhalo[ind] = s_mhhalo[ind_shuffle]

    # Add the central position
    add2sat(s_xgal,s_ygal,s_zgal,gtype,add=True)

    # Correct satellites again for periodic boundaries
    correct_periodic(s_xgal,s_ygal,s_zgal,gtype,lbox,halfbox=False) 
    print('Length xgal - length s_gal = {}'.format(len(xgal)-len(s_xgal)))

    # Change the mass of the host halo for the satellites
    mhhalo2sat(s_mhhalo,gtype) 

    # Save the shuffled information into a file
    outff = path+model+'iz'+sn+'/shuffled.hdf5'

    hf = h5py.File(outff,'w')
    hf.create_dataset('vols',data=vols)
    hf.create_dataset('jm',data=jm)
    hf.create_dataset('ihhalo',data=ihhalo)
    hf.create_dataset('gtype',data=gtype)
    hf.create_dataset('xgal',data=xgal)
    hf.create_dataset('ygal',data=ygal)
    hf.create_dataset('zgal',data=zgal)
    hf.create_dataset('mhhalo',data=mhhalo)
    hf.create_dataset('mgal',data=mgal)
    hf.create_dataset('s_xgal',data=s_xgal)
    hf.create_dataset('s_ygal',data=s_ygal)
    hf.create_dataset('s_zgal',data=s_zgal)
    hf.create_dataset('s_mhhalo',data=s_mhhalo)
    hf.close()

    print('Output: {}'.format(outff))
