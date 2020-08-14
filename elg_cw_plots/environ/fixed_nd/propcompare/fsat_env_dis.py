import sys,os.path
import subprocess
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = False

propname = 'fsat'
ytit = "Percentage of satellites"
ymin = 0. ; ymax = 50.

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
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15

xlabels = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']

cw_list = ['Vweb','Pweb'] 
nd_list = ['-2.0','-3.0','-4.2']
cuts = ['m','sfr','lo2']

cols = ['darkred','dodgerblue','palegreen']
symbols = ['s','x','+'] 

surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']
snaps = ['39','41']

if Testing:
    cw_list = ['Vweb']  ; nd_list = ['-2.0']
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI']
    snaps = ['39']

##########################################
# Output fraction summary
envsumfile = path+'env_files/'+model+'fsat_env_fractions.txt'
sumfile = open(envsumfile,'w')
sumfile.write('# Percentage of satellites, LSE, Selection \n')
sumfile.close()
sumfile = open(envsumfile,'a')

# Loop over the different files
for cw in cw_list:
    epath = path+'env_files/'+model+cw+'/'
    plotroot = path+'plots/'+model+'environ/props/'+cw+'/'+propname+'_'
    print('\n Plots: {}* \n'.format(plotroot))

    for iis,survey in enumerate(surveys1):
        inleg = ['All',survey,surveys2[iis]]
        numinleg = len(inleg)*len(cuts)
        lbar = dm*sep/numinleg

        iz = snaps[iis] 
        for nd in nd_list:
            # Initialize the parameters for the figures
            fig = plt.figure(figsize=(8.5,9.))
            jj = 111 ; ax = fig.add_subplot(jj)
            plt.xticks(xlabels,elabels)
            ax.tick_params(axis='x',which='minor',right=False,left=False)
            ax.set_ylabel(ytit,fontsize=fs) ; ax.set_ylim(ymin,ymax)

            if(iz == '41'):
                ztext = 'z = 0.83; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
            elif(iz == '39'):
                ztext = 'z = 0.99; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
            ax.text(emax-0.43*(emax-emin),ymax-0.04*(ymax-ymin), ztext) 

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

                inleg = ['All',survey,surveys2[iis]]
                files = [allfile,elgfile1,elgfile2]
                fprop = [allprop,elg1prop,elg2prop]

                for iif,efile in enumerate(files):
                    ii += 1
                    pfile = fprop[iif] 
                    
                    # Check if files exist and has more than one line
                    if (not os.path.isfile(efile) or not os.path.isfile(pfile)):
                        if Testing:print('Jumping {}'.format(efile))
                        continue
                    wcl_line = subprocess.check_output(["wc", "-l",efile])
                    wcl = int(wcl_line.split()[0])
                    pcl_line = subprocess.check_output(["wc", "-l",pfile])
                    pcl = int(pcl_line.split()[0])
                    if (wcl <= 1 or pcl <= 1):
                        if Testing:print('Jumping {}'.format(efile))
                        continue

                    # Read the file
                    xx, yy, zz, fenv = np.loadtxt(efile,unpack=True)
                    env = fenv.astype(int)
                    
                    # Read property file
                    # xgal 0, ygal 1, zgal 2 (Mpc/h), vxgal 3,vygal 4,vzgal 5 (Km/s),
                    # log10(massh) 6, log10(mass/Msun/h) 7, log10(sfr/Msun/h/Gyr) 8, 
                    # lum 9,lum_ext 10 (10^40 h^-2 erg/s),
                    # type 11 (0= Centrals; 1,2= Satellites) 
                    px, py, pz, gtype = np.loadtxt(pfile, usecols=(0,1,2,11), unpack=True)

                    #print(max(gtype)) ; sys.exit()

                    # Check that the coordinates have the same size
                    if ((len(xx) != len(px)) or 
                        (len(yy) != len(py)) or
                        (len(zz) != len(pz))):
                        print('[PROBLEM] Different lengths coordinates: {}\n {}\n'.
                              format(efile,pfile)) ; continue
                    # Check that the coordinates are ordered in the same way
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
                        gtotal = np.shape(ind)[1]

                        ind = np.where((env == ienv) & (gtype>0))                        
                        fsat = 100.*float(np.shape(ind)[1])/float(gtotal)

                        # Plot
                        xenv = ebins[ienv] + 0.5*dm*(1.-sep) + 0.5*lbar + ii*lbar

                        if (ienv == 0 and ic == 0):
                            parts = ax.scatter(xenv, fsat, c=cols[ic], 
                                               marker=symbols[iif],label=inleg[iif])
                        else:
                            parts = ax.scatter(xenv, fsat, c=cols[ic], 
                                               marker=symbols[iif])

                        # Write to summary file
                        root = efile.split('/')[-1]
                        sumfile.write('{:.2f} {} {} {}\n'.format(fsat,elabels[ienv], root.split('.dat')[0],cw))

            newnd = nd.replace('.','p')
            plotfile = plotroot+survey+'_nd'+newnd+'_sn'+iz+'_env2.pdf'

            # Legend
            leg = plt.legend(loc=2)
            for ll in leg.legendHandles:
                ll.set_color('k')
            leg.draw_frame(False)

            # Save figure
            fig.savefig(plotfile)
            plt.close()
            if Testing:
                print('Plot: {}'.format(plotfile))

sumfile.close()
print('Enviroment fractions: {}'.format(envsumfile))
