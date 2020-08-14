#! /usr/bin/env python

import numpy as np
import os.path, sys
from scipy import interpolate
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)
##################################

Testing = False

space = 'r' #'r' or 'z'

model = 'gp19/'

# Min. number of galaxy pairs to be considered
npairs = 1.

# Paths 
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
plotpath = inpath+'plots/'+model+'selections/xi/'
filepath = inpath+'selections/'

snapnum = ['39','41']
zzs = ['0.99','0.83']
surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']

cw_list = ['Vweb','Pweb']  

if Testing:
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI'] 
    snapnum = ['39'] ; zzs = ['0.99']
    cw_list = ['Vweb'] 
##################################

# Initialize the parameters for the figures
xmin = -2. ; xmax = 2.
ymin = 0. ; ymax = 8.

if (space == 'r'):
    xtit = "${\\rm log}_{10}(r/h^{-1}{\\rm Mpc})$"
    ytit = "DD" 
else:
    xtit = "${\\rm log}_{10}(s/h^{-1}{\\rm Mpc})$"
    ytit = "DD" 

cut = 'elg'

elabels = ['Total','Voids','Sheets','Filaments','Knots'] 
cols = ['k','greenyellow','limegreen','forestgreen','darkolivegreen']
lstyle = ['-','-','-','-','-'] 
lwidth = [4,2,2,2,2] 

##################################
# Loop over the different files
for cw in cw_list:
    for iiz,sn in enumerate(snapnum):
        zz = zzs[iiz] 

        if (space == 'r'):
            epath = filepath+model+'iz'+sn+'/'+space+'/'
    
        surveys = [surveys1[iiz],surveys2[iiz]] 
        for survey in surveys:
            # Initialize the parameters for the figures
            fig = plt.figure(figsize=(6.5,7.))
            jj = 111 ; ax = fig.add_subplot(jj)
            ax.set_autoscale_on(False) ;  ax.minorticks_on()
            ax.set_xlabel(xtit) ; ax.set_xlim(xmin,xmax) 
            ax.set_ylabel(ytit) ; ax.set_ylim(ymin,ymax)
            ztext = survey+', z = '+zz
            ax.text(xmax-0.4*(xmax-xmin),ymin+0.07*(ymax-ymin), ztext)
        
            # Loop over environments
            for iie,elabel  in enumerate(elabels):
                if (elabel == 'Total'):
                    efile = epath+cut+'cut_'+survey+\
                        '_sn'+sn+'_CUTExi_'+space+'.dat'
                else:
                    efile = epath+elabel+'_'+cw+'_'+cut+'cut_'+survey+\
                        '_sn'+sn+'_CUTExi_'+space+'.dat'

                if (not os.path.isfile(efile)):
                    continue
    
                # r   xi(r)   error(r)   DD(r)
                in_r,in_dd = np.loadtxt(efile, usecols=(0,3), unpack=True)

                ind = np.where((in_dd >npairs) & (in_r>=10**xmin))
                if (np.shape(ind)[1] > 0):
                    xg = np.log10(in_r[ind]) ; yg = np.log10(in_dd[ind])    
                    ax.plot(xg,yg,color=cols[iie],
                            linewidth=lwidth[iie],linestyle=lstyle[iie],
                            label=elabel)
    
            # Legends
            leg = plt.legend(loc=2) 
            for color,text in zip(cols,leg.get_texts()): 
                text.set_color(color) 
            leg.draw_frame(False)
    
            # Save figure
            plotfile = plotpath+'dd_'+space+'_env'+cw+\
                       '_'+survey+'_sn'+sn+'.pdf'
            fig.savefig(plotfile)
            print('Output: {}'.format(plotfile))
