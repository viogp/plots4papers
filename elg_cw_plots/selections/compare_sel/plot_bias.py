import sys, os.path
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

model = 'gp19/'
space = 'r'

sn_list = ['41','39']
surveys = ['All','DEEP2','VVDS-DEEP','eBOSS-SGC','DESI']
nds = ['-2.0','-3.0','-4.2']
cuts = ['m','sfr']
lst = ['-','--']

# Path
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
npath = inpath+'selections/'+model

# Prepare the plot
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 12.0 ; fs = 15  

# Find the variables
itop = 0 ; ibottom = 0
for iz in sn_list:
    # Input file
    infile = npath+'iz'+iz+'/'+space+'/'+'bias8.dat'
    bs, ns = np.loadtxt(infile, usecols=(0,2), unpack=True)
    cs, ss = np.loadtxt(infile, usecols=(3,4), unpack=True,dtype='str')

    # Plot
    plotfile = inpath+'plots/'+model+'xi/bias_sn'+iz+'.pdf'
    fig = plt.figure(figsize=(9.,9.))   
    xtit = 'nd' ; plt.xlabel(xtit)  
    ytit = 'bias' ; plt.ylabel(ytit)  
 
    for iis,survey in enumerate(surveys):
        for iic,cut in enumerate(cuts): 
            ind = np.where((ss == survey) & (cs == cut))
            if(np.shape(ind)[1]>0):
                x = ns[ind]
                y = bs[ind]
                inleg = survey+', '+cut+' cut'
                plt.plot(x,y,linestyle=lst[iic],label=inleg)

    # Save the figure with legend
    plt.legend(loc=1,frameon= False)
    fig.savefig(plotfile)
    print('Bias plot: {}'.format(plotfile))
    plt.close()
