import sys,os
from glob import glob
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from distinct_colours import get_distinct
import mpl_style
plt.style.use(mpl_style.style1)

model = 'gp19'

#cw = 'Vweb'
cw = 'Pweb'

#snap = '41' ; survey = 'eBOSS'
#snap = '41' ; survey = 'VVDS-DEEP'
#snap = '39' ; survey = 'DEEP2'
snap = '39' ; survey = 'DESI'

endings = ['_m8.5_sn'+snap+'.dat','_s3.0_sn'+snap+'.dat']
root = ['mcut_','scut_']
##########################################

path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
epath = path+'env_files/'+cw+'/'+model+'/'
plotfile = path+'plots/'+model+'/environ/'+cw+'_'+survey+'_persatcen.pdf'

##########################################
inleg = ['Mass cut','SFR cut',survey+' (mass)',survey+' (SFR)']
nfiles = len(inleg)

if(snap == '41'):
    ztext = 'z = 0.83'
elif(snap == '39'):
    ztext = 'z = 0.99'

# Plot histogram: 0 Void; 1 sheet; 2 filament; 3 Knots.
fig = plt.figure(figsize=(8.5,9.))
jj = 111 ; ax = fig.add_subplot(jj)

x = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']
plt.xticks(x,elabels) 
ax.tick_params(axis='x',which='minor',bottom='off')

ytit = "Fraction"
ax.set_ylabel(ytit)
ymin = 0. ; ymax = 1. ; ax.set_ylim(ymin,ymax)

ax.text(-0.5, 0.95, ztext)

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
    per = hist/float(len(env))

    print(infile, per)
    xenv = ebins + 0.5*dm*(1.-frac) + 0.5*lbar + ii*lbar
    ax.bar(xenv, per, lbar, \
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
    per = hist/float(len(env))

    xenv = ebins + 0.5*dm*(1.-frac) + 0.5*lbar + ii*lbar

    ax.bar(xenv, per, lbar, \
               color=cols[ii], label=inleg[ii])
#    ax.bar(xenv, abs(ymin-slogden), lbar, bottom=clogden, \
#               color=cols[ii],hatch='//')

# Legend
leg = plt.legend(loc=1)
#for color,text in zip(cols,leg.get_texts()):
#    text.set_color(cols)
leg.draw_frame(False)

# Save figure
fig.savefig(plotfile)
print 'Output: ',plotfile
