import sys,os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

epath = '/gpfs/data/violeta/lines/cosmicweb/env_files/'

inleg = ['Mass cut','SFR cut','eBOSS (mass)','eBOSS (SFR)']
nfiles = len(inleg)

endings = ['_m8.5_sn41.dat','_s3.0_sn41.dat']
survey = 'eBOSS'
root = ['mcut_','scut_']

# Plot histogram: 0 Void; 1 sheet; 2 filament; 3 Knots.
fig = plt.figure(figsize=(8.5,9.))
jj = 111 ; ax = fig.add_subplot(jj)

x = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']
plt.xticks(x,elabels) 
ax.tick_params(axis='x',which='minor',bottom='off')

ytit = "$log_{10}(\Phi/ \\rm{h^{3} Mpc^{-3}})$"
ax.set_ylabel(ytit)
ymin = -7. ; ymax = -1. ; ax.set_ylim(ymin,ymax)

cols = get_distinct(nfiles)

# Bins
emin = -0.5 ; emax = 3.5 ; dm = 1. 
frac = 0.85 ; lbar = dm*frac/nfiles
ebins = np.arange(emin,emax, dm) 

vol = 500.**3. # Simulation volume in (Mpc/h)^3
for ii,ending in enumerate(endings):
    infile = epath+root[ii]+survey+ending
    if (not os.path.isfile(infile)):
        print('STOP: {} not found'.format(infile)) ; sys.exit()

    xgal,ygal,zgal,env = np.loadtxt(infile,unpack=True)

    hist, bins_edges = np.histogram(env, bins=np.append(ebins,emax))

    ind = np.where(hist > 0)
    logden = np.log10(hist[ind]/vol)

    print infile, logden
    xenv = ebins + 0.5*dm*(1.-frac) + 0.5*lbar + ii*lbar
    ax.bar(xenv, abs(ymin-logden), lbar, bottom=ymin, \
               color=cols[ii], label=inleg[ii])

for ending in endings:
    ii += 1

    fil = epath+survey+'_centrals'+ending #; print ii,fil
    if (not os.path.isfile(fil)):
        print('STOP: {} not found'.format(fil)) ; sys.exit()
    cenv = np.loadtxt(fil,usecols=[3],unpack=True)

    fil = epath+survey+'_satellites'+ending
    if (not os.path.isfile(fil)):
        print('STOP: {} not found'.format(fil)) ; sys.exit()
    senv = np.loadtxt(fil,usecols=[3],unpack=True)

    env = np.append(cenv,senv)

    hist, bins_edges = np.histogram(env, bins=np.append(ebins,emax))
    ind = np.where(hist > 0)
    logden = np.log10(hist[ind]/vol)

    xenv = ebins + 0.5*dm*(1.-frac) + 0.5*lbar + ii*lbar

    ax.bar(xenv, abs(ymin-logden), lbar, bottom=ymin, \
               color=cols[ii], label=inleg[ii])
#    ax.bar(xenv, abs(ymin-slogden), lbar, bottom=clogden, \
#               color=cols[ii],hatch='//')

# Legend
leg = plt.legend(loc=1)
#for color,text in zip(cols,leg.get_texts()):
#    text.set_color(cols)
leg.draw_frame(False)

# Save figure
plotfile = epath + 'stack_satcen.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile
