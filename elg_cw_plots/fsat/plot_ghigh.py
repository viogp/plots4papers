import sys
import numpy as np
from Cosmology import *
import matplotlib ; matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

prop = 'g_highlim' 
xtit = "$g$ (fixed high limit)"
xmin = 19.
xmax = 23.
#################

model = 'gp19'
#model = 'gp19.font'
#model = 'gp19.starvation'

snap_list = [41, 39] #MillGas

llim = np.arange(39.5,42.,0.5)
surveys = ['All']
for ll in llim:
    surveys.append('logL>'+str(ll))

#################

# Prep plot

path = '/cosma5/data/durham/violeta/lines/cosmicweb/fsat/'+model+'/'
outplot = path+prop+'.pdf'

ytit = ('% satellites (>x)')

fig = plt.figure(figsize=(8.,9.))
ax = plt.subplot()
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
ax.set_xlim(xmin,xmax) ; ax.set_ylim(0.,30.)

cols = get_distinct(len(surveys))

#################

# Loop over the snapshots, reading the percenage of satellites

lstyle = ['-','--']
for iz,sn in enumerate(snap_list):
    #infile = path+'test.dat'
    infile = path+prop+'_'+str(sn)+'.dat'

    # Read data
    data = np.loadtxt(infile)
    xdata = data[:,0]
    
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


