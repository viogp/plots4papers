import sys,os.path
import subprocess
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

propname = 'lo2'
xtit = '${\\rm log}_{10}(L{\\rm [OII]}/h^{-2}{\\rm erg\, s}^{-1})$'

# Initialize histogram  
lmin = 38. ; lmax = 46. ; dl = 0.1

##########################################
model = 'gp19/'
path = '/cosma6/data/dp004/dc-gonz3/lines/cosmicweb/'
proppath = path+'selections/'+model+'ascii_files/'
volume = np.power(500.,3.)
##########################################

cut = 'shuffled_elg'
# Initialize the parameters for the figures
xmin = 40.2 ; xmax = 43.7
ymin = -5.9 ; ymax = -1.

lbins = np.arange(lmin,lmax,dl)  
xhist = lbins + dl*0.5 

ytit = "${\\rm log}_{10}(\Phi/h^3{\\rm Mpc}^{-3}{\\rm dex}^{-1})$" 

ylabels = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']

#cols = ['k','greenyellow','limegreen','forestgreen','darkolivegreen'] 
cols = ['k','grey','#e7d4e8','r','#7fbf7b','b','#af8dc3','c','#1b7837','m']
#lstyle = ['-','-','--','-.',':']
lstyle = ['-','-','-','-','-']
slstyle = [':']
lwidth = [4,2.5,2.5,2.5,2.5] 

surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']
snaps = ['39','41']
zzs = ['0.99','0.83']
cw_list = ['Vweb'] #,'Pweb']

if Testing:
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI']
    snaps = ['39']
    cw_list = ['Vweb'] 

##########################################

# Loop over the different files
for cw in cw_list:
    epath = path+'env_files/'+model+cw+'/'
    plotroot = path+'plots/'+model+'environ/props/'+cw+'/elgs/shuffled/'
    print('\n Plots: {}* \n'.format(plotroot))

    for iiz, sn in enumerate(snaps):
        surveys = [surveys1[iiz],surveys2[iiz]]

        for survey in surveys:
            if (survey=='DESI'):
                xmin = 41. ; xmax = 43.
                ymin = -5.9 ; ymax = -2.

            # Initialize the parameters for the figures
            fig = plt.figure(figsize=(6.5,7.))
            jj = 111 ; ax = fig.add_subplot(jj)
            ax.set_autoscale_on(False) ;  ax.minorticks_on() 
            ax.set_xlabel(xtit) ; ax.set_xlim(xmin,xmax)
            ax.set_ylabel(ytit) ; ax.set_ylim(ymin,ymax)
            #ax.text(xmin+0.02*(xmax-xmin),ymin+0.02*(ymax-ymin), ztext) 

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

            # Original files
            oend = 'elgcut_'+survey+'_sn'+sn+'.dat'
            oefile = epath+oend
            opfile = proppath+oend

            if (not os.path.isfile(oefile) or not os.path.isfile(opfile)):
                if Testing: print('Jumping {} or {}'.format(oefile,opfile))
                continue
            wcl_line = subprocess.check_output(["wc", "-l",oefile])
            wcl = int(wcl_line.split()[0])
            pcl_line = subprocess.check_output(["wc", "-l",opfile])
            pcl = int(pcl_line.split()[0])
            if (wcl <= 1 or pcl <= 1):
                if Testing: print('Jumping {} or {}'.format(oefile,opfile))
                continue

            # Read property file
            # xgal 0, ygal 1, zgal 2 (Mpc/h), vxgal 3,vygal 4,vzgal 5 (Km/s),
            # log10(massh) 6, log10(mass/Msun/h) 7, log10(sfr/Msun/h/Gyr) 8, 
            # log10(lum/h^-2 erg/s) 9, log10(lum_ext/h^-2 erg/s) 10,
            # type 11 (0= Centrals; 1,2= Satellites) 
            px, py, pz, lum_ext = np.loadtxt(pfile, usecols=(0,1,2,10), unpack=True)
            opx, opy, opz, olum_ext = np.loadtxt(opfile, usecols=(0,1,2,10), unpack=True)

            # Plot total 
            H, bins_edges = np.histogram(olum_ext,bins=np.append(lbins,lmax)) 
            yhist = H/dl/volume
            ind = np.where(yhist>0.)                
            if (np.shape(ind)[1]>1):
                ax.plot(xhist[ind],np.log10(yhist[ind]),
                        color=cols[0],linestyle=lstyle[0],
                        linewidth=lwidth[0],label='All '+survey) #+', z = '+zzs[iiz])

            H, bins_edges = np.histogram(lum_ext,bins=np.append(lbins,lmax)) 
            yhist = H/dl/volume
            ind = np.where(yhist>0.)                
            if (np.shape(ind)[1]>1):
                ax.plot(xhist[ind],np.log10(yhist[ind]),
                        color=cols[1],linestyle=slstyle[0],
                        linewidth=lwidth[-1],label='Shuffled')

            # Read the environmental file
            xx, yy, zz, fenv = np.loadtxt(efile,unpack=True)
            env = fenv.astype(int)

            # Check that the coordinates have the same size
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

            oxx, oyy, ozz, fenv = np.loadtxt(oefile,unpack=True)
            oenv = fenv.astype(int)

            # Check that the coordinates have the same size
            if ((len(oxx) != len(opx)) or 
                (len(oyy) != len(opy)) or
                (len(ozz) != len(opz))):
                print('[PROBLEM] Different lengths coordinates: {}\n {}\n'.
                      format(oefile,opfile)) ; continue
            # Chech that the coordinates are ordered in the same way
            if ((not np.allclose(oxx,opx,atol=1e-08,equal_nan=True)) or 
                (not np.allclose(oyy,opy,atol=1e-08,equal_nan=True)) or
                (not np.allclose(ozz,opz,atol=1e-08,equal_nan=True))):
                print('[PROBLEM] Files with different coordinates: {}\n {}\n'.
                      format(oefile,opfile)) ; continue

            # Loop over type of environment
            icol=1
            for iie, ienv in enumerate(np.unique(env)):
                # Plot original hist 
                icol += 1
                ind = np.where(oenv == ienv)
                prop = olum_ext[ind]
                H, bins_edges = np.histogram(prop,bins=np.append(lbins,lmax)) 
                yhist = H/dl/volume
                ind = np.where(yhist>0.)                
                if (np.shape(ind)[1]>1):
                    ax.plot(xhist[ind],np.log10(yhist[ind]),
                            color=cols[icol],linestyle=lstyle[iie+1],
                            linewidth=lwidth[iie+1],label=elabels[iie])
                
                # Plot original hist 
                icol += 1
                ind = np.where(env == ienv)
                prop = lum_ext[ind] #print(prop) ; sys.exit()
                H, bins_edges = np.histogram(prop,bins=np.append(lbins,lmax)) 
                yhist = H/dl/volume
                ind = np.where(yhist>0.)                
                if (np.shape(ind)[1]>1):
                    ax.plot(xhist[ind],np.log10(yhist[ind]),
                            color=cols[icol],linestyle=slstyle[0],
                            linewidth=lwidth[-1],label='S_'+elabels[iie])

            # Legend
            leg = plt.legend(loc=1)
            for color,text in zip(cols,leg.get_texts()):
                text.set_color(color) 
            leg.draw_frame(False)

            # Save figure
            plotfile = plotroot+'lfo2_'+survey+'_sn'+sn+'.pdf'
            if (Testing): print(plotfile)
            fig.savefig(plotfile)
            plt.close()

