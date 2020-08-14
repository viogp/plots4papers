# Output the x,y,z original and shuffled coordinates
# for the Millenium-II simulation

import os.path, sys    
import h5py
import numpy as np
from Corrfunc.theory.xi import xi # For test plot
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec 
import mpl_style 
plt.style.use(mpl_style.style1)

Testing = False

sn_list = ['39','41']
nd = 0.01   # Number density cut
lbox = 500.

if Testing:
    sn_list = ['39']

############################################

path = '/cosma6/data/dp004/dc-gonz3/Galform_Out/v2.7.0/stable/MillGas/'
model = 'gp19/'

############################################

nbin = 70
mmin, mmax = 9., 16.

for iz, sn in enumerate(sn_list):
    inff = path+model+'iz'+sn+'/shuffled.hdf5'

    hf = h5py.File(inff,'r')
    xgal = hf['xgal'][:]
    ygal = hf['ygal'][:]
    zgal = hf['zgal'][:]
    mgal = hf['mgal'][:]
    s_xgal = hf['s_xgal'][:]
    s_ygal = hf['s_ygal'][:]
    s_zgal = hf['s_zgal'][:]
    hf.close()

    # Plot to test the shuffling
    mass2cut = np.sort(mgal)[::-1] # Sort in reverse order
    mcut = mass2cut[int(nd*lbox**3)]
    print('Mcut for plot ={} ({} galaxies)'.format(mcut,int(nd*lbox**3)))

    mask = (mgal>mcut)
    x1 = xgal[mask]
    y1 = ygal[mask]
    z1 = zgal[mask]

    x2 = s_xgal[mask]
    y2 = s_ygal[mask]
    z2 = s_zgal[mask]

    # Plot the clustering as a test
    fig, ax0 = plt.subplots(nrows=1, ncols=1, figsize=(10,10))
    gs = gridspec.GridSpec(2, 1)
    gs.update(wspace=0., hspace=0.)
    ax = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0],sharex=ax)

    rbins = np.logspace(-2.,2.,40)
    rbins_mid = rbins[:-1] + (rbins[1]-rbins[0])/2.

    nthreads = 4
    xi1_cf = xi(lbox, nthreads, rbins, x1, y1, z1)
    xi2_cf = xi(lbox, nthreads, rbins, x2, y2, z2)
    ratios = xi1_cf['xi']/xi2_cf['xi']
    print('Clustering ratios = {}'.format(ratios))

    ax.set_xlim([-1.,1.5]) ; ax.set_xlim([-2.,3.5])    
    ax.set_ylabel(r'$\rm log_{10}\xi(r)$')
    ax.plot(np.log10(rbins_mid), np.log10(xi1_cf['xi']),label='nd='+str(nd))
    ax.plot(np.log10(rbins_mid), np.log10(xi2_cf['xi']),label='Shuffled')
    
    ax2.set_ylim([0.5,1.5])
    ax2.set_ylabel(r'$\rm log(r/h^{-1}Mpc$')
    ax2.set_ylabel(r'$\rm Ratios$')
    ax2.plot([-3,3], [1,1], color='k',ls=':',linewidth=1)
    ax2.plot(np.log10(rbins_mid), ratios)

    ax.legend(loc='best',frameon=False)
    plotff = path+model+'iz'+sn+'/shuffled.pdf'
    plt.savefig(plotff,bbox_inches='tight')
    print('Plot at {}'.format(plotff))
