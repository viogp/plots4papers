import sys
import numpy as np
from Cosmology import *
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

prop = 'lo2' 
xtit = "${\\rm log}_{10}(L\\rm{[OII]}/h^{-2}erg\, s^{-1})$" 
xmin = 38.
xmax = 43.
#################

model = 'gp19'
#model = 'gp19.font'
#model = 'gp19.starvation'

snap_list = [41, 39] #MillGas

surveys = ['All','eBOSS-SGC','VVDS-DEEP','DEEP2','DESI']

#################

# Prep plot

path = '/cosma5/data/durham/violeta/lines/cosmicweb/fsat/'+model+'/'
outplot = path+prop+'.pdf'

ytit = ('% satellites (>x)')

fig = plt.figure(figsize=(8.,9.))
ax = plt.subplot()
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
ax.set_xlim(xmin,xmax) ; ax.set_ylim(0.,30.)

# Skip second colour to match the other plots
cols1 = get_distinct(len(surveys)+1)
cols = get_distinct(len(surveys))
for i in range(2,len(cols1)):
    cols[i-1] = cols1[i]

#################

# Loop over the snapshots, reading the percenage of satellites

lstyle = ['-','--']
for iz,sn in enumerate(snap_list):
    #infile = path+'test.dat'
    infile = path+prop+'_'+str(sn)+'.dat'

    # Read cosmology from header line
    with open(infile,'r') as ff:
        l1 = ff.readline() ; l2 = ff.readline()
        l3 = ff.readline()
        zz = float(l3.split('redshift = ')[1].split(',')[0].strip())
        omega0 = float(l3.split('omega0 = ')[1].split(',')[0].strip())
        omegab = float(l3.split('omegab = ')[1].split(',')[0].strip())
        lambda0 = float(l3.split('lambda0 = ')[1].split(',')[0].strip())
        h0 = float(l3.split('h0 = ')[1].split(',')[0].strip())

    set_cosmology(omega0=omega0, omegab=omegab, 
                  lambda0=lambda0, h0 =h0,
                  universe="Flat", include_radiation=False)

    # Read data
    data = np.loadtxt(infile)
    xdata = data[:,0]
    
    # Print out logL and corresponding fluxes
    for x in np.arange(xmin,xmax,0.5):
        print('logL= {} , Flux(z={})= {}'.format(x,zz,logL2flux(x,zz))) 

    for isy,survey in enumerate(surveys):
        y = data[:,isy+1]
        ind = np.where(y >= 0.)
        if (np.shape(ind)[1]>1):
            yy = y[ind]
            xx = xdata[ind]

            if (iz==0):
                ax.plot(xx,yy,color=cols[isy],linestyle=lstyle[iz],label=survey)
            else:
                ax.plot(xx,yy,color=cols[isy],linestyle=lstyle[iz])

leg = plt.legend(loc=1)
for color,text in zip(cols,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)

# Save figure
plt.savefig(outplot)
print('Output: {}'.format(outplot))


