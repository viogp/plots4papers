import sys,os.path
import subprocess
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

Testing = True

model = 'gp19/'
path = '/cosma6/data/dp004/dc-gonz3/lines/cosmicweb/'
##########################################

# Bins and separation of histograms in the plot
emin = -0.5 ; emax = 3.5 ; dm = 1.
sep = 0.85 
ebins = np.arange(emin,emax, dm)

# Initialize the parameters for the figures
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['axes.labelsize'] = 10.0 ; fs = 15

xlabels = np.array([0,1,2,3])
elabels = ['Voids','Sheets','Filaments','Knots']

ytit = 'Fraction'
ymin = 0. ; ymax = 1.

cuts = ['m','sfr','lo2']
cols = ['darkred','dodgerblue','palegreen']
hatching = [' ','/','o'] 

surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']
snap_list = ['39','41']

cw_list = ['Vweb','Pweb']
nd_list = ['-2.0','-3.0','-4.2']

if Testing:
    surveys1 = ['VVDS-DEEP'] ; surveys2 = ['eBOSS-SGC'] ;snap_list = ['41']
    #surveys1 = ['DEEP2'] ; surveys2 = ['DESI'] ; snap_list = ['39']
    cw_list = ['Pweb'] 
    #nd_list = ['-4.2']
    
##########################################
# Output fraction summary
envsumfile = path+'env_files/'+model+'env_fractions.txt'
sumfile = open(envsumfile,'w')
sumfile.write('# Fraction in Voids,Sheets,Filaments,Knots \n')
sumfile.close()
sumfile = open(envsumfile,'a')

# Loop over the different files
for cw in cw_list:
    epath = path+'env_files/'+model+cw+'/'

    for iis,survey in enumerate(surveys1):
        inleg = ['Mass cut, All','Mass cut, '+survey,'Mass cut, '+surveys2[iis],
                 'SFR cut, All','SFR cut, '+survey,'SFR cut, '+surveys2[iis],
                 'L[OII] cut, All','L[OII] cut, '+survey,'L[OII] cut, '+surveys2[iis]]
        numinleg = len(inleg)
        lbar = dm*sep/numinleg

        iz = snap_list[iis]

        for nd in nd_list:
            # Initialize the parameters for the figures
            fig = plt.figure(figsize=(8.5,9.))
            jj = 111 ; ax = fig.add_subplot(jj)
            plt.xticks(xlabels,elabels)
            ax.tick_params(axis='x',which='minor',bottom=False)
            ax.set_ylabel(ytit,fontsize=fs) ; ax.set_ylim(ymin,ymax)

            if(iz == '41'):
                ztext = 'z = 0.83; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
            elif(iz == '39'):
                ztext = 'z = 0.99; $10^{'+nd+'}{\\rm Mpc}^{-3}{\\rm h}^{3}$'
            ax.text(1.7, 0.95, ztext)

            ii = -1
            for ic, cut in enumerate(cuts):
                allfile = epath+cut+'cut_All_nd'+nd+'_sn'+iz+'.dat'
                elgfile1 = epath+cut+'cut_'+survey+\
                         '_nd'+nd+'_sn'+iz+'.dat'
                elgfile2 = epath+cut+'cut_'+surveys2[iis]+\
                         '_nd'+nd+'_sn'+iz+'.dat'
                files = [allfile,elgfile1,elgfile2]

                for iie,efile in enumerate(files):
                    ii += 1

                    # Check if exists and has more than one line
                    if (not os.path.isfile(efile)):
                        continue
                    wcl_line = subprocess.check_output(["wc", "-l",efile])
                    wcl = int(wcl_line.split()[0])
                    if (wcl <= 1):
                        continue

                    # Read the file
                    xx,yy,zz,fenv = np.loadtxt(efile,unpack=True)
                    env = fenv.astype(int)

                    # Histograms with the fraction of galaxies
                    hist, bins_edges = np.histogram(env, bins=np.append(ebins,emax))
                    frac = hist/float(len(env))
                    if (ii < numinleg and 'All' in efile): 
                        sumfile.write('{}_{} : {} \n'.format(cw,efile.split('/')[-1],frac))
                    elif ('All' not in efile):
                        sumfile.write('{}_{} : {} \n'.format(cw,efile.split('/')[-1],frac))

                    # Plot
                    xenv = ebins + 0.5*dm*(1.-sep) + 0.5*lbar + ii*lbar
                    if ('All' in efile):
                        ax.bar(xenv, frac, lbar, \
                               color=cols[ic], label=inleg[ii])
                    else:
                        ax.bar(xenv, frac, lbar, \
                               color=cols[ic], label=inleg[ii],\
                               hatch=hatching[iie],fill=False,edgecolor=cols[ic])

            newnd = nd.replace('.','p')
            plotfile = path+'plots/'+model+'environ/'+cw+\
                       '_'+survey+'_nd'+newnd+'_sn'+iz+'_env2.pdf'

            # Legend
            leg = plt.legend(loc=2)
            leg.draw_frame(False)

            # Save figure
            fig.savefig(plotfile)
            print('Plot: {}'.format(plotfile))
            plt.close()

sumfile.close()
print('Enviroment fractions: {}'.format(envsumfile))
