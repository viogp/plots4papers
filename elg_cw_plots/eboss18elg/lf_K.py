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

path = '/cosma5/data/durham/violeta/Galform_Out/v2.7.0/stable/MillGas/'
nvol = 64

#plotdir = '/cosma5/data/durham/violeta/lines/desi_hod_o2/plots/modelplots/gp18.bh.'
plotdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/modelplots/thiswork.' 
#plotdir = '/cosma5/data/durham/violeta/lines/cosmicweb/plots/modelplots/nominfrac.' 
models = ['gp19','gp19.font','gp19.starvation','gp17'] 
#models= ['gp18','gp18.e0p01nominfrac','gp18.e0.nominfrac','gp18.e0p1.nominfrac']
inleg = ['This work','10% stripping ','Starvation','GP18']
#inleg = models

# Bands for LF
iband = 'UKIRT-K' ; iband6 = 'K'

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
        infile = path+model+'/iz61/ivol'+str(ivol)+'/galaxies.hdf5'

        if (os.path.isfile(infile)):
            f = h5py.File(infile,'r')
            volh = volh + f['Parameters/volume'].value
            f.close()
        else:
            print 'Not found: ',infile

        infile = path+model+'/iz61/ivol'+str(ivol)+'/tosedfit.hdf5'

        if (os.path.isfile(infile)):
            f = h5py.File(infile,'r')
            group = f['Output001']

            try:
                mag  = group['mag_'+iband+'_r_tot'].value
            except:
                mag  = group['mag'+iband6+'r_tot'].value
            H, bins_edges = np.histogram(mag,bins=np.append(mbins,mmax))
            plf[index,:] = plf[index,:] + H

            try:
                mage = group['mag_'+iband+'_r_tot_ext'].value
            except:
                mage = group['mag'+iband6+'r_tot_ext'].value
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
xmin = -15.
xmax = -25.
ymin = -5.5
ymax = -1.

xtit = "${{\\rm M_{AB}("+iband+")}\, -\, 5log_{10}{\\rm h}}$"
ytit = "$log_{10}(\Phi/ \\rm{h^{3} Mpc^{-3} mag^{-1}})$"
fig = plt.figure(figsize=(8.5,9.))
jj = 111
ax = fig.add_subplot(jj) ; ax.set_autoscale_on(False)
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)
#ax.tick_params(labelsize=fs-2)   
start, end = ax.get_xlim()

cols = get_distinct(len(models)+1) 
colors = cols #; g = ['grey'] ; colors.extend(g)
colors[len(colors)-1] = 'k'        

# Observational data
dobs = '/cosma/home/dphlss/violeta/Galform2/galform/Obs_Data2/'
fobs = 'lfk_z0_driver12.data'

# Plot observations
file = dobs+fobs
# Driver observations
mag, den, err, num = np.loadtxt(file,unpack=True)
ind = np.where(den > 0.)
x = mag[ind]
y = np.log10(den[ind]/0.5) # Observations are per 0.5 mag
eh = np.log10(den[ind]+err[ind]) - np.log10(den[ind])
el = np.log10(den[ind]) - np.log10(den[ind]-err[ind]) 
ind = np.where(np.isinf(el) | np.isnan(el))
el[ind] = 999.
ax.errorbar(x,y,yerr=[el,eh],fmt='o', s=3, 
            ecolor='k',color='k', mec='k', label='GAMA, Driver+2012')

# For Chi2
ind = np.where((den > 0.) & (mag<-15.) & (mag>-23.))
ox = mag[ind]-0.089 
oy = den[ind] ; oerr = err[ind] 

# Plot predictions
lsty = ['-','--',':','-']   
lwdt = [3.,1.5,1.5,1.5]
for im in range(len(models)):
    py = plfe[im,:]
    ind = np.where(py < 0.)
    x = xlf[ind] ; y = 10**py[ind] 
    ax.plot(x,py[ind],cols[im],label=inleg[im],\
                linestyle=lsty[im],linewidth=lwdt[im])
    
    # Chi2
    my = np.interp(ox,x,y)
    chi = (my-oy)**2/oerr**2
    rchi2 = np.sum(chi)/len(chi)
    print inleg[im],'  Chi2',np.sum(chi),rchi2
         
# Legend
leg = plt.legend(loc=3)
for color,text in zip(colors,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)


# Save figures
#fig.tight_layout()
plotfile = plotdir + 'lf_K.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile

#gp14.r595   Chi2 4886.93721148 305.433575718
#gp15newmg   Chi2 12965.566347 810.34789669
#This work   Chi2 14285.2265612 892.826660072
#10% stripping    Chi2 13067.0347076 816.689669225
#Starvation   Chi2 12523.1485008 782.696781302
#GP18   Chi2 12967.3617295 810.460108094
