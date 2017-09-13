#! /usr/bin/env python
#~/lines/desi_hod_o2/modelplots/ssfr_m.py 
import numpy as np
import h5py, os.path, sys
from Cosmology import *
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from distinct_colours import get_distinct,pault_cmap
import matplotlib.gridspec as gridspec
from scipy import ndimage
from stats import *
import mpl_style
plt.style.use(mpl_style.style1) ; ptmap=pault_cmap(1)

path6 = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/'
path = '/gpfs/data/violeta/Galform_Out/v2.6.0/aquarius_trees/MillGas/'
nvol = 64

plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/'
models = ['gp14','gp15newmg']
inleg = ['GP14','This work']

#plotdir = '/gpfs/data/violeta/lines/desi_hod_o2/plots/modelplots/anders.'
#models = ['gp15newmg','gp15newmg.anders'] ; inleg=models

# Initialize GSMF
mmin = 8.5
mmax = 15.
dm = 0.1
mbins = np.arange(mmin,mmax,dm)
mhist = mbins + dm*0.5
# Initialize SSFR
smin = -5.
smax = 3.
ds = 0.1
sbins = np.arange(smin,smax,ds)
shist = sbins + ds*0.5

nlevel = 10

gsmf  = np.zeros(shape=(len(models),len(mbins)))
ssfrf = np.zeros(shape=(len(models),len(sbins)))
smf   = np.zeros(shape=(len(models),len(sbins),len(mbins)))

###########################################
# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()


# Loop over redshifts
for index,model in enumerate(models):
    print model
    volume = 0. ; first = True
    for ivol in range(nvol):
        gfile = path+model+'/iz61/ivol'+str(ivol)+'/galaxies.hdf5'
        if(os.path.isfile(gfile)):
            # Read the relevant information from galaxies.hdf5
            f = h5py.File(gfile,'r') #; print gfile
            vol1 = f['Parameters/volume'].value ; volume = volume + vol1
            h0 =   f['Parameters/h0'].value 
            omega0 = f['Parameters/omega0'].value
            omegab = f['Parameters/omegab'].value
            lambda0 =f['Parameters/lambda0'].value


            group = f['Output001']
            zz   = group['redshift'].value
            set_cosmology(omega0=omega0,omegab=omegab, \
                              lambda0=lambda0,h0=h0, \
                              universe="Flat",include_radiation=False)
            slim = 1./tHubble(zz)            

            mdisk = group['mstars_disk'].value
            mbulge = group['mstars_bulge'].value
            mass1 = mdisk + mbulge
            if (model=='gp14'):
                mass1 = mass1/0.81

            sdisk = group['mstardot'].value # Msolar/h/Gyr
            sbulge = group['mstardot_burst'].value
            sfr1 = sdisk + sbulge
            #if (model=='gp14'):
            #    sfr1 = sfr1/1.26

            ind = np.where((mass1>10**mmin) & (sfr1>10**smin))
            mass = np.log10(mass1[ind])
            sfr = np.log10(sfr1[ind])
            ssfr = sfr - mass

            # GSMF
            H, bins_edges = np.histogram(mass,bins=np.append(mbins,mmax))
            gsmf[index,:] = gsmf[index,:] + H

            # sSFR
            H, bins_edges = np.histogram(ssfr,bins=np.append(sbins,smax))
            ssfrf[index,:] = ssfrf[index,:] + H

            # sSFR-GSMF
            H, xedges, yedges = np.histogram2d(ssfr,mass,\
                                                       bins=[np.append(sbins,smax),np.append(mbins,mmax)])
            smf[index,:,:] = smf[index,:,:] + H

            f.close()
        else:
            print 'NOT found:',gfile
            
    if (volume>0.):
        gsmf[index,:] = gsmf[index,:]/volume/dm  # In Mpc^3/h^3
        ssfrf[index,:] = ssfrf[index,:]/volume/ds  
        smf[index,:,:] = smf[index,:,:]/volume/dm/ds  
        print 'Side of the explored box (Mpc/h) = ',pow(volume,1./3.)

# Figure http://matplotlib.org/users/gridspec.html
fig = plt.figure(figsize=(8.5,9.))
gs = gridspec.GridSpec(3, 3)
gs.update(wspace=0., hspace=0.)
ax = plt.subplot(gs[1:,:-1])
cols = get_distinct(len(models)) 
colors = cols ; g = ['grey'] ; colors.extend(g)

# SFRF vs M
xtit="$log_{10}(\\rm M_{*}/M_{\odot}h^{-1})$"
ytit="$log_{10}(\\rm sSFR/Gyr^{-1})$"
xmin=mmin ; xmax=11.9 ; ymin=smin ; ymax=0.9
ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)


# Plot Franx+08 cut
y = [np.log10(slim),np.log10(slim)] ; x= [xmin,xmax]
ax.plot(x,y,'k--')
y = [np.log10(0.3*slim),np.log10(0.3*slim)] ; x= [xmin,xmax]
ax.plot(x,y,'k:')

for ii in range(len(models)):
    zz = np.zeros(shape=(len(shist),len(mhist))) 
    pz = smf[ii,:,:] ; ind = np.where(pz>0.)
    zz = np.log10(pz) 

    # Heat map
    #x,y = np.meshgrid(np.append(mbins,mmax),np.append(sbins,smax))
    #ax.pcolormesh(x, y, zz, cmap=plt.get_cmap(ptmap))
    # Contours
    xx,yy = np.meshgrid(mbins,sbins)
    #zz = 100.*zz ; al = [5.,16.,68.,95.] 
    al = [-3.5,-2.5,-1.5] 
    # Smooth
    #zz = ndimage.gaussian_filter(zz, sigma=0.5)
    cs = ax.contour(xx, yy, zz, levels=al, \
                        linestyles='-',colors=cols[ii])
    cs.levels = [nf(val) for val in cs.levels]
    ax.clabel(cs, cs.levels, inline=1,inline_spacing=0,\
                  fontsize=10,fmt='%r')#fmt='%r %%')

# GSMF ###################################
axm = plt.subplot(gs[0, :-1],sharex=ax)
ytit="$log_{10}(\Phi)$" ; axm.set_ylabel(ytit)
axm.set_autoscale_on(False) ;  axm.minorticks_on()
axm.set_ylim(-5.5,-1.) 
plt.setp(axm.get_xticklabels(), visible=False)

# Plot observations from Baldry+2012
dobs = '/cosma/home/violeta/Galform2/galform-2.6.0/Obs_Data2/'
file = dobs+'mf/mf_baldry_2012.txt'
oh = 0.7 
lm,p3,dp3 = np.loadtxt(file,usecols=[0,1,2],unpack=True)            
xobs = lm + np.log10(oh)
yobs = xobs*0. - 999.
indx = np.where( p3 > 0)
yobs[indx] = np.log10(p3[indx]) -3. -3.*np.log10(oh)
lerr = yobs*0. - 999. 
indx = np.where( (p3-dp3) > 0)
lerr[indx]  = np.log10(p3[indx] - dp3[indx])-3. -3.*np.log10(oh)
herr = yobs*0. + 999.
indx = np.where( (p3+dp3) > 0)
herr[indx]  = np.log10(p3[indx] + dp3[indx])-3. -3.*np.log10(oh)
plt.errorbar(xobs, yobs, yerr=[yobs-lerr,herr-yobs], fmt='o', ecolor='grey',color='grey', mec='grey',label="Baldry+2012, z<0.06")

# Models
for ii in range(len(models)):
    py = gsmf[ii,:] ; ind = np.where(py>0.)
    x = mhist[ind] ; y = np.log10(py[ind])
    ind = np.where(y < 0.)
    axm.plot(x[ind],y[ind],color=cols[ii],label=inleg[ii])


# SFRF
axs = plt.subplot(gs[1:, 2],sharey=ax)
xtit="$log_{10}(\Phi)$" ; axs.set_xlabel(xtit)
axs.set_autoscale_on(False) ;  axs.minorticks_on()
axs.set_xlim(-4.4,0.0)
start, end = axs.get_xlim()
axs.xaxis.set_ticks(np.arange(-4., end, 1.))
plt.setp(axs.get_yticklabels(), visible=False)

## Gruppioni+15 data, z=0.15
#ofil = '/cosma/home/violeta/compare_sams/nifty/tuning/Obs_Data/sfrf_z0.15.dat'
##log(SFR/(Msun/yr))_low  log(SFR/(Msun/yr))_high  log(Phi/Mpc^3 dex^-1) error (Gruppioni+ 2015, z=[0.0,0.3])
#oh = 0.71 
#ls,hs,lphi,error = np.loadtxt(file,unpack=True)            

for ii in range(len(models)):
    px = ssfrf[ii,:] ; ind = np.where(px>0.)
    y = shist[ind] ; x = np.log10(px[ind])
    ind = np.where(x < 0.)
    axs.plot(x[ind],y[ind],color=cols[ii],label=inleg[ii])

# Legend
leg = axs.legend(bbox_to_anchor=(1., 1.4),fontsize='small')
for color,text in zip(cols,leg.get_texts()):
    text.set_color(color)
    leg.draw_frame(False)


# Save figures
plotfile = plotdir + 'ssfr_m.pdf'
fig.savefig(plotfile)
print 'Output: ',plotfile

