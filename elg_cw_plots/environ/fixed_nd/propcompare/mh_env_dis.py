import sys,os.path
import subprocess
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

propname = 'lmh'
xtit = "${\\rm log}_{10}(M_{\\rm halo}/h^{-1}{\\rm M}_{\odot})$"
xmin = 10.5 ; xmax = 15.

##########################################

model = 'gp19/'
path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
proppath = path+'selections/'+model+'ascii_files/'
##########################################

# Separation of environments in the plot
emin = -0.5 ; emax = 3.5 ; dm = 1.
sep = 0.85 
ebins = np.arange(emin,emax, dm)  

# Initialize the parameters for the figures
ylabels = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']

cuts = ['m','sfr','lo2']

cols = ['darkred','dodgerblue','palegreen']
lstyle = [':','--','-']

surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']
snaps = ['39','41']
cw_list = ['Vweb','Pweb']
nd_list = ['-2.0','-3.0','-4.2']

if Testing:
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI']
    snaps = ['39'] #; cuts = ['lo2'] 
    cw_list = ['Vweb'] ; nd_list=['-3.0']

##########################################

# Loop over the different files
for cw in cw_list:
    epath = path+'env_files/'+model+cw+'/'
    plotroot = path+'plots/'+model+'environ/props/'+cw+'/'+propname+'_'
    print('\n Plots: {}* \n'.format(plotroot))

    for iis,survey in enumerate(surveys1):
        inleg = [surveys2[iis],survey,'All']
        numinleg = len(inleg)*len(cuts)
        lbar = dm*sep/numinleg

        iz = snaps[iis] 
        for nd in nd_list:
            # Initialize the parameters for the figures
            fig = plt.figure(figsize=(7.,7.))
            jj = 111 ; ax = fig.add_subplot(jj)
            plt.yticks(ylabels,elabels)
            #plt.yticks(ylabels,[' ',' ',' ',' '])
            ax.tick_params(axis='y',which='minor',right=False,left=False)
            ax.set_xlabel(xtit) ; ax.set_xlim(xmin,xmax)

            if(iz == '41'):
                ztext = 'z = 0.83; $10^{'+nd+'}h^3{\\rm Mpc}^{-3}$'
            elif(iz == '39'):
                ztext = 'z = 0.99; $10^{'+nd+'}h^3{\\rm Mpc}^{-3}$'
            ax.text(xmax-0.5*(xmax-xmin),emax-0.02*(emax-emin), ztext) 

            ii = -1
            for ic, cut in enumerate(cuts):
                end = cut+'cut_All_nd'+nd+'_sn'+iz+'.dat'
                allfile = epath+end
                allprop = proppath+end

                end = cut+'cut_'+survey+'_nd'+nd+'_sn'+iz+'.dat'
                elgfile1 = epath+end
                elg1prop = proppath+end

                end = cut+'cut_'+surveys2[iis]+'_nd'+nd+'_sn'+iz+'.dat'
                elgfile2 = epath+end
                elg2prop = proppath+end

                files = [elgfile2,elgfile1,allfile]
                fprop = [elg2prop,elg1prop,allprop]

                for iif,efile in enumerate(files):
                    ii += 1
                    pfile = fprop[iif] 

                    # Check if files exist and has more than one line
                    if (not os.path.isfile(efile) or not os.path.isfile(pfile)):
                        if Testing: print('Jumping {}'.efile)
                        continue
                    wcl_line = subprocess.check_output(["wc", "-l",efile])
                    wcl = int(wcl_line.split()[0])
                    pcl_line = subprocess.check_output(["wc", "-l",pfile])
                    pcl = int(pcl_line.split()[0])
                    if (wcl <= 1 or pcl <= 1):
                        if Testing: print('Jumping {}'.efile)
                        continue

                    # Read the file
                    xx, yy, zz, fenv = np.loadtxt(efile,unpack=True)
                    env = fenv.astype(int)
                    
                    # Read property file
                    # xgal 0, ygal 1, zgal 2 (Mpc/h), vxgal 3,vygal 4,vzgal 5 (Km/s),
                    # log10(massh) 6, log10(mass/Msun/h) 7, log10(sfr/Msun/h/Gyr) 8, 
                    # lum 9,lum_ext 10 (10^40 h^-2 erg/s),
                    # type 11 (0= Centrals; 1,2= Satellites) 
                    px, py, pz, lmh = np.loadtxt(pfile, usecols=(0,1,2,6), unpack=True)

                    #print(min(lmass)) ; sys.exit()

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
                    for ienv in np.unique(env):
                        ind = np.where(env == ienv)
                        if (np.shape(ind)[1] <= 1):
                            continue
                        prop = lmh[ind]

                        # Plot
                        yenv = ebins[ienv] + 0.5*dm*(1.-sep) + 0.5*lbar + ii*lbar
                        parts = ax.violinplot(prop, [yenv],vert=False, widths=0.2,
                                              showmeans=False, showmedians=False,showextrema=False)

                        # Change colour of the body of the plot
                        for pc in parts['bodies']:
                            pc.set_color(cols[ic])
                        # Change the colour of the lines
                        quartile1, median, quartile3 = np.percentile(prop, [10, 50, 90])
                        # Median
                        ax.scatter(median,yenv,marker='|', color=cols[ic], s=100)
                        # Quartiles
                        ax.hlines(yenv, quartile1, quartile3, color=cols[ic], 
                                  linestyle='-', lw=5)
                        # Extremes
                        if ((ic == 0) and (ienv == 0)):
                            ax.hlines(yenv, np.min(prop), np.max(prop), 
                                      color=cols[ic], linestyle=lstyle[iif], lw=2,
                                      label=inleg[iif])
                        else:
                            ax.hlines(yenv, np.min(prop), np.max(prop), 
                                      color=cols[ic], linestyle=lstyle[iif], lw=2)

            newnd = nd.replace('.','p')
            plotfile = plotroot+survey+'_nd'+newnd+'_sn'+iz+'_env2.pdf'

            # Legend
            leg = plt.legend(loc=4)
            for ll in leg.legendHandles:
                ll.set_color('k')
            leg.draw_frame(False)

            # Save figure
            fig.savefig(plotfile)
            plt.close()
            if Testing: 
                print('Plot: {}'.format(plotfile))

