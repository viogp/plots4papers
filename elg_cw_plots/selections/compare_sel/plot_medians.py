import sys, os.path
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt

model = 'gp19/'

sn_list = ['41','39']
surveys = ['All','VVDS-DEEP','DEEP2','eBOSS-SGC','DESI'] 
nds = ['-2.0','-3.0','-4.2']
cuts = ['m','sfr','lo2']

verbose = False

# Paths
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
ndpath = inpath+'selections/'

# Input file
infile = ndpath+model+'medians.dat'
dinfo = np.loadtxt(infile,usecols=(0,),unpack=True,dtype='str')
massh,mass,sfr,lum_ext = np.loadtxt(infile,usecols=(1,2,3,4),unpack=True)

# Prepare the plot
oplot = ndpath+model+'medians.pdf'
symbols=['*','^','s','p','o']
colours=['navy','royalblue','lightsteelblue']
pattern=['none','full','left']
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col') 
ax3.set_xlabel('M*') ; ax4.set_xlabel('SFR')
ax1.set_ylabel('Mhalo') ; ax3.set_ylabel('Mhalo')
ax2.set_ylabel('L[OII]') ; ax4.set_ylabel('L[OII]')

# Find the variables
itop = 0 ; ibottom = 0
for ii,info in enumerate(dinfo):
    dd = info.split('cut')[0]
    iic = cuts.index(dd)

    dd = info.split('_')[1]
    iis = surveys.index(dd)

    ddd = info.split('nd')[1]
    dd = ddd.split('_')[0]
    iin = nds.index(dd)

    sn = info.split('sn')[1]

    # Plot the medians
    if (sn == '41'):
        if (itop == 0):
            ax1.plot(massh[ii],mass[ii],color=colours[iin],
                     marker=symbols[iis],fillstyle=pattern[iic],
                     label='sn='+sn)
            ax2.plot(sfr[ii],lum_ext[ii],color=colours[iin],
                     marker=symbols[iis],fillstyle=pattern[iic],
                     label='sn='+sn)
            itop += 1
        else:
            ax1.plot(massh[ii],mass[ii],color=colours[iin],
                     marker=symbols[iis],fillstyle=pattern[iic])
            ax2.plot(sfr[ii],lum_ext[ii],color=colours[iin],
                     marker=symbols[iis],fillstyle=pattern[iic])
    else:
        if (ibottom == 0):
            ax3.plot(massh[ii],mass[ii],color=colours[iin],
                     marker=symbols[iis],fillstyle=pattern[iic],
                     label='sn='+sn)
            ax4.plot(sfr[ii],lum_ext[ii],color=colours[iin],
                     marker=symbols[iis],fillstyle=pattern[iic],
                     label='sn='+sn)
            ibottom += 1
        else:
            ax3.plot(massh[ii],mass[ii],color=colours[iin],
                     marker=symbols[iis],fillstyle=pattern[iic])
            ax4.plot(sfr[ii],lum_ext[ii],color=colours[iin],
                     marker=symbols[iis],fillstyle=pattern[iic])
        
ax1.legend(loc=2,frameon= False)
ax2.legend(loc=2,frameon= False)
ax3.legend(loc=2,frameon= False)
ax4.legend(loc=2,frameon= False)
fig.savefig(oplot)
print('Plot with medians: {}'.format(oplot))
