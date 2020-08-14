import sys,os.path
import subprocess
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import ndimage
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

propname = 'sfr_m'

##########################################
model = 'gp19/'
path = '/cosma6/data/dp004/dc-gonz3/lines/cosmicweb/'
proppath = path+'selections/'+model+'ascii_files/'
hodpath = path+'hod/'
volume = np.power(500.,3.)
##########################################

cut = 'elg'

al = np.array([-4.2,-3.,-2.]) 

# Initialize the GSMF
mmin = 8.5 ; mmax = 15. ; dm = 0.1
mbins = np.arange(mmin,mmax,dm) 
mhist = mbins + dm*0.5
# Initialize the sSFR
smin = 3. ; smax = 13. ; ds = 0.1
sbins = np.arange(smin,smax,ds)
shist = sbins + ds*0.5  

ylabels = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']

#cols = ['k','greenyellow','limegreen','forestgreen','darkolivegreen']
cols = ['k','#e7d4e8','#7fbf7b','#af8dc3','#1b7837']
lwidth = [4,2.5,2.5,2.5,2.5] 

surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']
snaps = ['39','41']
zzs = ['0.99','0.83']
cw_list = ['Vweb','Pweb']

if Testing:
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI']
    snaps = ['39']
    cw_list = ['Vweb'] 

##########################################

# Loop over the different files
for cw in cw_list:
    epath = path+'env_files/'+model+cw+'/'
    plotroot = path+'plots/'+model+'environ/props/'+cw+'/elgs/'
    print('\n Plots: {}* \n'.format(plotroot))

    for iiz, sn in enumerate(snaps):
        surveys = [surveys1[iiz],surveys2[iiz]]

        for survey in surveys:
            # Initialize the parameters for the figures
            fig = plt.figure(figsize=(7.,7.))
            gs = gridspec.GridSpec(3, 3)
            gs.update(wspace=0., hspace=0.) 
            ax = plt.subplot(gs[1:,:-1]) 

            # Fig. SFRF vs M   
            xtit="${\\rm log}_{10}(M_{*}/h^{-1}{\\rm M}_{\odot})$" 
            ytit="${\\rm log}_{10}(SFR/h^{-1}{\\rm M}_{\odot}{\\rm Gyr}^{-1})$"  
            if (survey == 'DESI'):
                xmin=mmin ; xmax=11.9 ; ymin=8.5 ; ymax=11. 
            else:
                xmin=mmin ; xmax=11.9 ; ymin=6. ; ymax=smax 
            ax.set_xlim(xmin,xmax) ; ax.set_ylim(ymin,ymax) 
            ax.set_xlabel(xtit) ; ax.set_ylabel(ytit)

            # GSMF   
            axm = plt.subplot(gs[0, :-1],sharex=ax) 
            ytit="$log_{10}(\Phi)$" ; axm.set_ylabel(ytit) 
            axm.set_autoscale_on(False) ;  axm.minorticks_on()
            if (survey == 'DESI'):
                axm.set_ylim(-5.5,-2.)
            else:
                axm.set_ylim(-5.5,-1.)
            plt.setp(axm.get_xticklabels(), visible=False)  

            # SFRF
            axs = plt.subplot(gs[1:, 2],sharey=ax)  
            xtit="$log_{10}(\Phi)$" ; axs.set_xlabel(xtit) 
            axs.set_autoscale_on(False) ;  axs.minorticks_on()
            if (survey == 'DESI'):
                axs.set_xlim(-4.4,-2.)
            else:
                axs.set_xlim(-4.4,0.0)
            start, end = axs.get_xlim() 
            axs.xaxis.set_ticks(np.arange(-4., end, 1.))
            plt.setp(axs.get_yticklabels(), visible=False) 

            # Files to read
            end = cut+'cut_'+survey+'_sn'+sn+'.dat'
            efile = epath+end
            pfile = proppath+end

            # Check if files exist and has more than one line
            if (not os.path.isfile(efile) or not os.path.isfile(pfile)):
                if Testing: print('Jumping {} or {}'.format(efile,pfile))
                continue
            wcl_line = subprocess.check_output(["wc", "-l",efile])
            wcl = int(wcl_line.split()[0])
            pcl_line = subprocess.check_output(["wc", "-l",pfile])
            pcl = int(pcl_line.split()[0])
            if (wcl <= 1 or pcl <= 1):
                if Testing: print('Jumping {} or {}'.format(efile,pfile))
                continue

            # Read property file
            # xgal 0, ygal 1, zgal 2 (Mpc/h), vxgal 3,vygal 4,vzgal 5 (Km/s),
            # log10(massh) 6, log10(mass/Msun/h) 7, log10(sfr/Msun/h/Gyr) 8, 
            # log10(lum/h^-2 erg/s) 9, log10(lum_ext/h^-2 erg/s) 10,
            # type 11 (0= Centrals; 1,2= Satellites) 
            px, py, pz, mass, sfr = np.loadtxt(pfile, usecols=(0,1,2,7,8), unpack=True)


            # Total values: GSMF
            H, bins_edges = np.histogram(mass,bins=np.append(mbins,mmax))
            gsmf = H/volume/dm  # In Mpc^3/h^3 

            # Total values: sSFR
            H, bins_edges = np.histogram(sfr,bins=np.append(sbins,smax))
            sfrf = H/volume/ds  # In Mpc^3/h^3 

            # Total: sSFR-GSMF
            H, xedges, yedges = np.histogram2d(sfr,mass,                     
                                               bins=[np.append(sbins,smax),    
                                                     np.append(mbins,mmax)])
            smf = H/volume/dm/ds

            # Plot SMF vs SFR                                                  
            matplotlib.rcParams['contour.negative_linestyle'] = '-'      
            zz = np.zeros(shape=(len(shist),len(mhist))) ; zz.fill(-999.)      
            ind = np.where(smf>0.)                                             
            zz[ind] = np.log10(smf[ind])                                       
                                                                                   
            ind = np.where(zz>-999.)                                           
            if (np.shape(ind)[1]>1):                                           
                xx,yy = np.meshgrid(mbins,sbins) # Contours
                cs = ax.contour(xx, yy, zz, levels=al,
                                colors=cols[0],
                                linewidths=lwidth[0])                          
                                                                                   
            # GSMF                                                             
            ind = np.where(gsmf>0.) ; y = np.log10(gsmf[ind])
            x = mhist[ind] 
            ind = np.where(y < 0.)                                             
            axm.plot(x[ind],y[ind],color=cols[0],
                     linestyle='-',linewidth=lwidth[0])

            # SFRF                                                             
            ind = np.where(sfrf>0.) ; x = np.log10(sfrf[ind])
            y = shist[ind] 
            ind = np.where(x < 0.)
            inleg = survey+', z='+zzs[iiz]
            axs.plot(x[ind],y[ind],color=cols[0],
                     linestyle='-',linewidth=lwidth[0],
                     label=inleg) 

            # Read the environmental file
            xx, yy, zz, fenv = np.loadtxt(efile,unpack=True)
            env = fenv.astype(int)

            # Chech that the coordinates have the same size
            if ((len(xx) != len(px)) or 
                (len(yy) != len(py)) or
                (len(zz) != len(pz))):
                print('[PROBLEM] Different lengths coordinates: {}\n {}\n'.
                      format(efile,pfile)) ; continue
            # Chech that the coordinates are ordered in the same way
            if ((not np.allclose(xx,px,atol=1e-08,equal_nan=True)) or 
                (not np.allclose(yy,py,atol=1e-08,equal_nan=True)) or
                (not np.allclose(zz,pz,atol=1e-08,equal_nan=True))):
                print('[PROBLEM] Files with different coordinates: {}\n {}\n'.
                      format(efile,pfile)) ; continue

            # Loop over type of environment
            for iie, ienv in enumerate(np.unique(env)):
                ind = np.where(env == ienv)
                emass = mass[ind]
                esfr = sfr[ind]

                # Total values: GSMF
                H, bins_edges = np.histogram(emass,bins=np.append(mbins,mmax))
                gsmf = H/volume/dm  # In Mpc^3/h^3 

                # Total values: sSFR
                H, bins_edges = np.histogram(esfr,bins=np.append(sbins,smax))
                sfrf = H/volume/ds  # In Mpc^3/h^3 

                # Total: sSFR-GSMF
                H, xedges, yedges = np.histogram2d(esfr,emass,                     
                                                   bins=[np.append(sbins,smax),    
                                                         np.append(mbins,mmax)])
                smf = H/volume/dm/ds

                # Plot SMF vs SFR                                              
                zz = np.zeros(shape=(len(shist),len(mhist))) ; zz.fill(-999.)      
                ind = np.where(smf>0.)                                             
                zz[ind] = np.log10(smf[ind])                                       
                                                                                   
                ind = np.where(zz>-999.)                                           
                if (np.shape(ind)[1]>1):                                           
                    xx,yy = np.meshgrid(mbins,sbins) # Contours
                    cs = ax.contour(xx, yy, zz, levels=al,
                                    colors=cols[iie+1],
                                    linewidths=lwidth[iie+1])
                                                                                   
                # GSMF                                                             
                ind = np.where(gsmf>0.)                                  
                x = mhist[ind] ; y = np.log10(gsmf[ind])                             
                ind = np.where(y < 0.)                                             
                axm.plot(x[ind],y[ind],color=cols[iie+1],
                         linestyle='-',linewidth=lwidth[iie+1])

                # SFRF                                                             
                ind = np.where(sfrf>0.)                                  
                y = shist[ind] ; x = np.log10(sfrf[ind])                             
                ind = np.where(x < 0.)                                             
                axs.plot(x[ind],y[ind],color=cols[iie+1],
                         linestyle='-',linewidth=lwidth[iie+1],
                         label=elabels[iie]) 

            # Legend
            leg = axs.legend(bbox_to_anchor=(1., 1.4),fontsize='small',
                             handlelength=0, handletextpad=0)                      
            for item in leg.legendHandles:                                         
                item.set_visible(False)                                            
            for color,text in zip(cols,leg.get_texts()):                        
                text.set_color(color)                                              
                leg.draw_frame(False) 

            # Save figure
            plotfile = plotroot+'sfr_m_'+survey+'_sn'+sn+'.pdf'
            if (Testing): print(plotfile)
            fig.savefig(plotfile)
            plt.close()
