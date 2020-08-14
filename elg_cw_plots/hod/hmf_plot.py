import numpy as np
import os.path, sys
import matplotlib ; matplotlib.use('Agg') 
from matplotlib import pyplot as plt

path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
snap_list = [39, 41]
model = 'gp19'

plotfile = path+'plots/'+model+'/hmf.pdf'

# Initialize the parameters for the figure ------------------ 
plt.rcParams['legend.numpoints'] = 1 
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 20 
plt.rcParams['lines.linewidth'] = 2 
fig = plt.figure(figsize=(8.5,9.))  
xtit = "${\\rm log}_{10}(\\rm{M/M_{\odot}}h^{-1})$" 
ytit = "${\\rm log}_{10}(\Phi/ Mpc^{-3}h^3 {\\rm dlog}_{10}M)$"  
xmin = 10. ; xmax = 16. 
ymin = -6.5 ; ymax = 0.  
jj = 111 
ax = fig.add_subplot(jj) 
#ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit,fontsize = fs) ; ax.set_ylabel(ytit,fontsize = fs) #-------------------------------------------------------------

# Loop over the redshifts of interest
for iz,zsnap in enumerate(snap_list):
    infil = path+'hod/'+model+'/hmf_sn'+str(zsnap)+'.txt'
    #infil = path+'hod/'+model+'/hmf_sn'+str(zsnap)+'_2vols.txt'
    print(infil)
    # Read the volume and dl of the simulation
    with open(infil,'r') as ff:
        for il,line in enumerate(ff):
            if il == 2:
                if(line.split()[1] != 'Simulation'):
                    sys.exit('STOP unexpected header \n'+infil+' \n')
                volume = float(line.split()[-1])**3.
            elif il == 4:
                if(line.split()[1] != 'dl'):
                    sys.exit('STOP unexpected header \n'+infil+' \n')
                dl = float(line.split()[-1])
                break

    # Read the midpoint Mh and the number of haloes
    mh, nrh = np.loadtxt(infil,usecols=(0,3),unpack=True)
    y = nrh/volume/dl
    ind = np.where(y>0.)
    ax.plot(mh[ind],np.log10(y),label='sn='+str(zsnap))

ax.legend(loc=1,frameon= False) 
fig.savefig(plotfile)
print('Output: {}'.format(plotfile))
