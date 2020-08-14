import sys,os.path
import subprocess
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

propname = 'hod'

##########################################
model = 'gp19/'
path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
proppath = path+'selections/'+model+'ascii_files/'
hodpath = path+'hod/'
volume = np.power(500.,3.)
##########################################

cut = 'elg'
# Initialize the parameters for the figures
xmin = 11. ; xmax = 15.
ymin = -3. ; ymax = 2.

xtit = '${\\rm log}_{10}(M_{\\rm halo}/h^{-1}{\\rm M}_{\odot})$'
ytit = '${\\rm log}_{10}\\langle N_{\\rm M}\\rangle$'

ylabels = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']

#cols = ['k','greenyellow','limegreen','forestgreen','darkolivegreen']
cols = ['k','#e7d4e8','#7fbf7b','#af8dc3','#1b7837'] 
lstyle = ['-',':','--']
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
        # Read the number of haloes 
        hmffile = hodpath+model+'hmf_sn'+sn+'.txt' # Full volume
        if Testing: hmffile = hodpath+model+'hmf_sn'+sn+'_2vols.txt'

        mhist,mlow,mhigh,nhs = np.loadtxt(hmffile,unpack='True')  

        surveys = [surveys1[iiz],surveys2[iiz]]
        for survey in surveys:
            if (survey == 'DESI'):
                ymax = 0.2
            ztext = 'All '+survey+', z = '+zzs[iiz]
            # Initialize the parameters for the figures
            fig = plt.figure(figsize=(6.5,7.))
            jj = 111 ; ax = fig.add_subplot(jj)
            ax.set_autoscale_on(False) ;  ax.minorticks_on() 
            ax.set_xlabel(xtit) ; ax.set_xlim(xmin,xmax)
            ax.set_ylabel(ytit) ; ax.set_ylim(ymin,ymax)
            #ax.text(xmax-0.5*(xmax-xmin),ymax-0.05*(ymax-ymin), ztext) 

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
            px, py, pz, massh, mass, x = np.loadtxt(pfile, usecols=(0,1,2,6,7,11), unpack=True)
            gtype = x.astype(int)

            # Initialize array for hod
            hod = np.zeros(shape=(3,len(nhs))) ; hod.fill(-999.) 

            # Calculate HODs
            for ii,nh in enumerate(nhs):
                # All
                ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii])) 
                numa = np.shape(ind)[1] 
                if (numa >= 1 and nh > 0):
                    hod[0,ii] = float(numa)/nh

                # Centrals
                ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]) & (gtype == 0)) 
                numc = np.shape(ind)[1] 
                if (numc >= 1 and nh > 0):
                    hod[1,ii] = float(numc)/nh

                # Satellites
                ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]) & (gtype > 0)) 
                nums = np.shape(ind)[1] 
                if (nums >= 1 and nh > 0):
                    hod[2,ii] = float(nums)/nh

            #Plot the 3 HODs:
            for ii in range(3):
                if (ii == 1): continue #skip centrals
                yy = hod[ii,:]
                ind = np.where(yy>0.)
                x = mhist[ind]
                y = np.log10(yy[ind])  

                if (ii == 0):
                    ax.plot(x,y,color=cols[0],linestyle=lstyle[ii],
                            linewidth=lwidth[0],label=ztext)
                else:
                    ax.plot(x,y,color=cols[0],linestyle=lstyle[ii],
                            linewidth=lwidth[0])

            # Read the environmental file
            xx, yy, zz, fenv = np.loadtxt(efile,unpack=True)
            env = fenv.astype(int)

            # ChecK that the coordinates have the same size
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
                # Initialize array for hod
                hod = np.zeros(shape=(3,len(nhs))) ; hod.fill(-999.) 

                # Calculate HODs
                for ii,nh in enumerate(nhs):
                    # All
                    ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]) &
                                   (env == ienv)) 
                    numa = np.shape(ind)[1] 
                    if (numa >= 1 and nh > 0):
                        hod[0,ii] = float(numa)/nh

                    # Centrals
                    ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]) & 
                                   (gtype == 0)  & (env == ienv)) 
                    numc = np.shape(ind)[1] 
                    if (numc >= 1 and nh > 0):
                        hod[1,ii] = float(numc)/nh

                    # Satellites
                    ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]) & 
                                   (gtype > 0)  & (env == ienv)) 
                    nums = np.shape(ind)[1] 
                    if (nums >= 1 and nh > 0):
                        hod[2,ii] = float(nums)/nh

                #Plot the 3 HODs:
                for ii in range(3):
                    if (ii == 1): continue #skip centrals
                    yy = hod[ii,:]
                    ind = np.where(yy>0.)
                    x = mhist[ind]
                    y = np.log10(yy[ind])  

                    if (ii == 0):
                        ax.plot(x,y,color=cols[iie+1],linestyle=lstyle[ii],
                            linewidth=lwidth[iie+1],label=elabels[iie])
                    else:
                        ax.plot(x,y,color=cols[iie+1],linestyle=lstyle[ii],
                            linewidth=lwidth[iie+1])

            # Legend
            leg = plt.legend(loc=2)
            for item in leg.legendHandles:
                item.set_visible(False) 
            for color,text in zip(cols,leg.get_texts()):
                text.set_color(color) 
            leg.draw_frame(False)

            # Save figure
            plotfile = plotroot+'hod_'+survey+'_sn'+sn+'.pdf'
            if (Testing): print(plotfile)
            fig.savefig(plotfile)
            plt.close()
