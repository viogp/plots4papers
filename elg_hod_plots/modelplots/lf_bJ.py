#! /usr/bin/env python

import numpy as np
import os.path,sys
import h5py
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from distinct_colours import get_distinct 
import mpl_style
plt.style.use(mpl_style.style1)

path6 = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/'
nvol = 64
# 
#plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/'
#models = ['gp14.r595','gp15newmg']
#inleg = ['GP14','This work']

plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/anders.'
models = ['gp15newmg','gp15newmg.anders'] ; inleg=models

# Bands for LF
iband = 'BJ'

# Limits for LFs
mmin = -40.
mmax = -3.
dm   = 0.25
mbins = np.arange(mmin,mmax,dm)
xlf = mbins + dm*0.5

plf  = np.zeros(shape = (len(models),len(mbins)))
plfe = np.zeros(shape = (len(models),len(mbins)))

vols = np.zeros(shape = (len(models)))

for index,model in enumerate(models):
    volh = 0.
    for ivol in range(nvol):
        infile = path6+model+'/iz61/ivol'+str(ivol)+'/galaxies.hdf5'

        if (os.path.isfile(infile)):
            f = h5py.File(infile,'r')
            volh = volh + f['Parameters/volume'].value
            f.close()
        else:
            print 'Not found: ',infile

        infile = path6+model+'/iz61/ivol'+str(ivol)+'/tosedfit.hdf5'

        if (os.path.isfile(infile)):
            f = h5py.File(infile,'r')
            group = f['Output001']

            mag  = group['mag'+iband+'r_tot'].value
            H, bins_edges = np.histogram(mag,bins=np.append(mbins,mmax))
            plf[index,:] = plf[index,:] + H

            mage = group['mag'+iband+'r_tot_ext'].value
            H, bins_edges = np.histogram(mage,bins=np.append(mbins,mmax))
            plfe[index,:] = plfe[index,:] + H

            f.close()
        else:
            print 'Not found: ',infile

    # Take logs and normalize
    print 'Total volume considered = (', np.power(volh,1./3.), ' Mpc/h)^3'

    ind = np.where(plf[index,:] > 0)
    plf[index,ind] = np.log10(plf[index,ind]/dm/volh)

    ind = np.where(plfe[index,:] > 0)
    plfe[index,ind] = np.log10(plfe[index,ind]/dm/volh)

######################
# Make plots
xmin = -14.
xmax = -23.
ymin = -5.5
ymax = -1.

xtit = "${{\\rm M_{AB}(b_J)}\, -\, 5log{\\rm h}}$"
ytit = "$log(\Phi/ \\rm{h^{3} Mpc^{-3} mag^{-1}})$"
fig = plt.figure(figsize=(8.5,9.))
jj = 111
ax = fig.add_subplot(jj) ; ax.set_autoscale_on(False)
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
#ax.tick_params(labelsize=fs-2)   
start, end = ax.get_xlim()

cols = get_distinct(len(models)) 
colors = cols ; g = ['grey'] ; colors.extend(g)

# Observational data
dobs = '/cosma/home/violeta/Galform2/galform-2.6.0/Obs_Data2/'
fobs ='lfBJ_norberg02.data' #'lfk_z0_driver12.data'

# Plot observations
file = dobs+fobs
# Norberg 2002
den, mag, err = np.loadtxt(file,usecols=[0,1,2],unpack=True)
ind = np.where(den > 0.)
ox = mag[ind]-0.089 #Chenkde to AB
oy = np.log10(den[ind])
oerr = err[ind]
eh = np.log10(den[ind]+oerr) - np.log10(den[ind])
el = np.log10(den[ind]) - np.log10(den[ind]-oerr) 
ind = np.where(np.isinf(el) | np.isnan(el))
el[ind] = 999.
ax.errorbar(ox,oy,yerr=[el,eh],fmt='o', ecolor='grey',color='grey', mec='grey', label='2dF, Norberg+2002')

# For Chi2
ind = np.where((den > 0.) & (mag<-14.) & (mag>-22.))
ox = mag[ind]-0.089 
oy = den[ind] ; oerr = err[ind] 

# Plot predictions
for im in range(len(models)):
    py = plfe[im,:]
    ind = np.where(py < 0.)
    x = xlf[ind] ; y = 10**py[ind] 
    ax.plot(x,py[ind],color=cols[im],label=inleg[im])

    # Chi2
    my = np.interp(ox,x,y)
    chi = (my-oy)**2/oerr**2
    rchi2 = np.sum(chi)/len(chi)
    print models[im],'  Chi2',np.sum(chi),rchi2
    
                
# Legend
leg = plt.legend(loc=3)
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)


# Save figures
#fig.tight_layout()
plotfile = plotdir + 'lf_bJ.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile

#Output information
#Total volume considered = ( 500.0  Mpc/h)^3
#Total volume considered = ( 500.0  Mpc/h)^3
#gp14.r595   Chi2 6049.95690098 208.619203482
#gp15newmg   Chi2 4308.89452314 148.582569763
#Output:  /gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/lf_bJ.pdf
