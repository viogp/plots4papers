#!/bin/env python

import os.path, sys
from numpy import *
import h5py
from Cosmology import *
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from distinct_colours import get_distinct

#### John's scripts###########
def density_image(x, y, min_count, max_count, shape, log=True, range=None):    
    """
    Generate a floating point image of a 2D particle distribution

    x, y                 - particle coordinates
    min_count, max_count - range of counts that map to 0 and 1 in
                           the output image
    shape                - shape of the output array, (nx,ny)
    log                  - whether to log scale the counts
    range                - range of coordinates to plot

    Returns:

    2D image array with values in range 0-1
    """
    if range is None:
        range = ((amin(x), amax(x)), (amin(y), amax(y)))

    image,xedge,yedge = histogram2d(x, y, range=range, bins=shape)
    image = asarray(image, dtype=float)
    if log:
        # Log scaling - negative values get set to zero in output
        ind = image > 0
        image[logical_not(ind)] = 0.0        
        lgmin = log10(min_count)
        lgmax = log10(max_count)
        image[ind] = (log10(image[ind]) - lgmin) / (lgmax-lgmin)
    else:
        # Linear scaling - clamp to range 0-1
        image = (image-min_count)-(max_count-min_count)
        image[image<0.0] = 0.0
        image[image>1.0] = 1.0
        
    return image.transpose()

def extract_particles(isnap, zlow, zup, lmin, lmax, shift, cm):

    basename = "/gpfs/data/Millgas/data/dm/500/snapdir_%03d/500_dm_%03d" % (isnap, isnap)

    # Read in particle data in a 10Mpc/h slice
    # First loop over snapshot files.
    x = [] ; y = []
    for ifile in range(512):

        # Read coordinates from one file
        fname = "%s.%d.hdf5" % (basename,ifile)
        print "Reading: %s" % fname
        f = h5py.File(fname, "r")
        pos_x = f["PartType1/Coordinates"].value[:,0]
        pos_y = f["PartType1/Coordinates"].value[:,1]
        pos_z = f["PartType1/Coordinates"].value[:,2]
        f.close()

        # Discard particles not in 0<z<10 slice and outside lmin,lmac
        ind = where((pos_z>zlow) & (pos_z<zup)  & \
                        (pos_x>lmin) & (pos_x<lmax) & \
                        (pos_y>lmin) & (pos_y<lmax) )
        x = append(x,pos_x[ind]) 
        y = append(y,pos_y[ind])
           
    # Make image as 2D float array
    img = density_image(x-shift,y-shift,min_count=0.1, max_count=10000, 
                        shape=(1024,1024), log=True)    

    # Display
    return plt.imshow(img, extent=(lmin-shift, lmax-shift, \
                                       lmin-shift, lmax-shift), \
                          cmap=cm, vmin=0.0, vmax=1.0, \
                          interpolation="nearest", origin="lower")

##############################
plotdir =  '/gpfs/data/violeta/lines/desi_hod_o2/plots/pretty_images/'
sn =['44','42','40','37']
zz = ['0.6','0.75','0.9','1.2']

snap = 39 # 0.988
zlow = 10.
zup  = 20.

path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/'
model = 'MillGas/gp15newmg/' #'MillGas/gp14/'
line = 'OII3727' ; lline = '[OII]'


# Plot only a section
#lmin = 400. ; lmax= 450. ; shift=400.
#lmin = 200. ; lmax= 250. ; shift=200.
#lmin = 300. ; lmax= 350. ; shift=300.
lmin = 0. ; lmax= 50. ; shift=0.

##############################

# DM
extract_particles(snap, zlow, zup, lmin, lmax, shift, plt.cm.Greys)
plt.xlabel("x(Mpc $h^{-1})$") ; plt.ylabel("y(Mpc $h^{-1})$")
plt.xlim((lmin-shift,lmax-shift)) ;plt.ylim((lmin-shift,lmax-shift))

##############################

# Look for OII emitters
bands = ['RK','m2','m2']
mcuts = [24.1,22.5,24]
fcuts = [2.7*10.**-17., 3.5*10.**-17., 1.9*10.**-17.]

inleg = ['Halo cut','DEEP2','VVDS-WIDE','VVDS-DEEP']
ntypes = len(inleg)
cols = get_distinct(5)

# Loop over volumes
firstpass = True 
xh =[] ; yh =[] ; lh =[]
xd2=[] ; yd2=[] ; ld2=[]
xvw=[] ; yvw=[] ; lvw=[]
xvd=[] ; yvd=[] ; lvd=[]
for ivol in range(64): #(64):
    # OII emitters 
    pfile = path+model+'iz'+str(snap)+'/ivol'+str(ivol)+'/galaxies.hdf5'
    if (os.path.isfile(pfile)):
        #print pfile
        # Get some of the model constants
        f = h5py.File(pfile,'r')
        zz   = f['Output001/redshift'].value
        xgal = f['Output001/xgal'].value
        ygal = f['Output001/ygal'].value
        zgal = f['Output001/zgal'].value
        group = f['Parameters']
        h0 = group['h0'].value 
        omega0 = group['omega0'].value
        omegab = group['omegab'].value
        lambda0 =group['lambda0'].value
        mhhalo = np.log10(f['Output001/mhhalo'].value)
        cen = f['Output001/type'].value 
        f.close()

        set_cosmology(omega0=omega0,omegab=omegab,lambda0=lambda0, \
                          h0=h0, universe="Flat",include_radiation=False)
        tomag = band_corrected_distance_modulus(zz)
        
        efile = path+model+'/iz'+str(snap)+'/ivol'+str(ivol)+'/elgs.hdf5'
        if (os.path.isfile(efile)):
            f = h5py.File(efile,'r')            

            # Haloes
            ind = where((cen<1) & (mhhalo>11.75) & \
                            (xgal>lmin) & (xgal<lmax) & \
                            (ygal>lmin) & (ygal<lmax) & \
                            (zgal>zlow) & (zgal<zup))

            xh = append(xh,xgal[ind])
            yh = append(yh,ygal[ind])
            lh = append(lh,mhhalo[ind])

            for index,ib in enumerate(bands):
                mag = f['Output001/mag'+ib+'o_tot_ext'].value + tomag
                lum_ext = f['Output001/L_tot_'+line+'_ext'].value
                icut = mcuts[index] ; fluxcut = fcuts[index]
                lcut = emission_line_luminosity(fluxcut,zz)
                if(len(mag) == len(zgal)): 
                    ind = where((mag<icut) & (lum_ext>lcut)  & \
                                    (xgal>lmin) & (xgal<lmax) &\
                                    (ygal>lmin) & (ygal<lmax) &\
                                    (zgal>zlow) & (zgal<zup))

                    if (index ==0): 
                        xd2 = append(xd2,xgal[ind])
                        yd2 = append(yd2,ygal[ind])
                        ll = np.log10(lum_ext[ind])
                        ld2 = append(ld2,ll)
                    elif (index ==1): 
                        xvw = append(xvw,xgal[ind])
                        yvw = append(yvw,ygal[ind])
                        ll = np.log10(lum_ext[ind])
                        lvw = append(lvw,ll)
                    else: 
                        xvd = append(xvd,xgal[ind])
                        yvd = append(yvd,ygal[ind])
                        ll = np.log10(lum_ext[ind])
                        lvd = append(lvd,ll)
                else:
                    if(index==0):
                        print 'Different lengths in elgs.hdf5 and galaxies.hdf5'
                        print len(mag),len(zgal)
                        print efile

            f.close()
    else:
        print pfile,' not found'

volume = 500.*500.*(zup-zlow)


# DEEP2
print 'DEEP2 density=',len(xd2)/volume
for i in range(len(xd2)):
    x = xd2[i]-shift
    y = yd2[i]-shift
    plt.plot(x,y,"o",c=cols[0],markeredgecolor=cols[0],\
                 markersize=ld2[i]*4.,alpha=0.8)

# VVDSWIDE
print 'VVDSWIDE density=',len(xvw)/volume
for i in range(len(xvw)):
    x = xvw[i]-shift
    y = yvw[i]-shift
    #plt.plot(x,y,"^",c=cols[1],markeredgecolor=cols[1],\
    #             markersize=lvw[i]*5.,alpha=0.8)

# Haloes
print 'Haloes density=',len(xh)/volume, min(lh),max(lh)
for i in range(len(xh)):
    x = xh[i]-shift
    y = yh[i]-shift
    plt.plot(x,y,"o",c='none',markeredgecolor='r',\
                 markersize=lh[i]/1.5,linewidth=3)


# Save figures
plotfile = plotdir+model+'deep2z0.9.pdf'
plt.savefig(plotfile)
print 'Output: ',plotfile
