import sys, os.path
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mpl_style
plt.style.use(mpl_style.style1)

model = 'gp19/'

sn_list = ['41','39'] ; zz_list = ['0.83','0.99']
surveys1 = ['VVDS-DEEP','DEEP2']  
surveys2 = ['eBOSS-SGC','DESI']
nds = ['-2.0','-3.0','-4.2']
cuts = ['m','sfr','lo2']

verbose = False

##############
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
hodpath = inpath+'hod/'

###############
# Initialize the parameters for the figures
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15
fig = plt.figure(figsize=(9.,16.))

xtit = "${\\rm log}_{10}(M_{\\rm halo}/M_{\odot}h^{-1})$"
ytit = "${\\rm log}_{10}(\\langle N_M\\rangle)$"

xmin = 10. ; xmax = 15.
ymin = -3. ; ymax = 2.

cols = ['navy','royalblue','palegreen']
cutlabel = ['Mass cut','SFR cut','L[OII] cut']

for iiz, sn in enumerate(sn_list):
    szz = zz_list[iiz]

    surveys = ['All',surveys1[iiz],surveys2[iiz]]

    ## Loop over the files
    for survey in surveys:
        plotfile = inpath+'plots/'+model+'selections/hod/hods_'+survey+'_sn'+sn+'.pdf'
        fig, ax = plt.subplots(len(nds), #row,colum
                               #sharex='all', sharey='all',
                               figsize=(9.,15.)) 

        for iind,nd in enumerate(nds):
            ifile = 0
            ax[iind].set_xlim(xmin,xmax) ; ax[iind].set_ylim(ymin,ymax)
            ax[iind].set_xlabel(xtit,fontsize = fs) ; ax[iind].set_ylabel(ytit,fontsize = fs)
            ax[iind].text(11.6,1.5,survey+', z='+szz+', nd='+nd, fontsize='small')

            for ic, cut in enumerate(cuts):
                hfile = hodpath+model+cut+'cut_'+survey+'_nd'+nd+'_sn'+sn+'.dat'
                # Check if the file exists
                if (not os.path.isfile(hfile)):
                    if verbose:
                        print('WARNING: {} not found'.format(infile))
                    continue
                ifile += 1
                mh,nall,ncen,nsat = np.loadtxt(hfile,unpack='True')

                # All
                ind = np.where(nall>-999.)
                ax[iind].plot(mh[ind],nall[ind],color=cols[ic],linestyle='-',
                        linewidth=2.5,label=cutlabel[ic])
                # Centrals
                ind = np.where(ncen>-999.)
                ax[iind].plot(mh[ind],ncen[ind],color=cols[ic],linestyle='--')
                # Satellites
                ind = np.where(nsat>-999.)
                ax[iind].plot(mh[ind],nsat[ind],color=cols[ic],linestyle=':')

            if (ifile > 0):
                leg = ax[iind].legend(loc=2, fontsize='small',
                                handlelength=0, handletextpad=0)
                for item in leg.legendHandles:
                    item.set_visible(False)
                for color,text in zip(cols,leg.get_texts()):
                    text.set_color(color)
                    leg.draw_frame(False)
            
        # Save figures
        fig.savefig(plotfile)
        print('Output: {}'.format(plotfile))

