import sys,os.path
import subprocess
import numpy as np
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

model = 'gp19/'
path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
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

cols = ['darkred','dodgerblue']
hatching = ['/','//'] 
##########################################
# Output fraction summary
envsumfile = path+'env_files/'+model+'env_fractions.txt'
sumfile = open(envsumfile,'w')
sumfile.write('# File : Fraction in Voids,Sheets,Filaments,Knots \n')
sumfile.close()
sumfile = open(envsumfile,'a')

# Loop over the different files
for cw in ['Vweb','Pweb']:
    epath = path+'env_files/'+model+cw+'/'
    
    for survey in ['eBOSS-SGC']:#,'DESI','eBOSS-SGC','VVDS-DEEP']:
        inleg = ['Mass cut',survey+' (mass)','SFR cut',survey+' (SFR)']
        numinleg = len(inleg)
        lbar = dm*sep/numinleg

        for iz in ['41']:#,'39']:
            for nd in ['-2.0','-3.0','-4.2']:
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
                ax.text(1.8, 0.95, ztext)

                ii = -1
                for ic, cut in enumerate(['m','sfr']):
                    allfile = epath+cut+'cut_All_nd'+nd+'_sn'+iz+'.dat'
                    elgfile = epath+cut+'cut_'+survey+\
                             '_nd'+nd+'_sn'+iz+'.dat'
                    files = [allfile,elgfile]

                    for efile in files:
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
                            sumfile.write('{} : {} \n'.format(efile.split('/')[-1],frac))
                        elif ('All' not in efile):
                            sumfile.write('{} : {} \n'.format(efile.split('/')[-1],frac))

                        # Plot
                        xenv = ebins + 0.5*dm*(1.-sep) + 0.5*lbar + ii*lbar
                        if ('All' in efile):
                            ax.bar(xenv, frac, lbar, \
                                   color=cols[ic], label=inleg[ii])
                        else:
                            ax.bar(xenv, frac, lbar, \
                                   color=cols[ic], label=inleg[ii],\
                                   hatch=hatching[ic],fill=False,edgecolor=cols[ic])

                newnd = nd.replace('.','p')
                plotfile = path+'plots/'+model+'environ/'+cw+\
                           '_'+survey+'_nd'+newnd+'_sn'+iz+'_env.pdf'

                # Legend
                leg = plt.legend(loc=2)
                #for color,text in zip(zip(cols,cols),leg.get_texts()):
                #  text.set_color(cols)
                leg.draw_frame(False)

                # Save figure
                fig.savefig(plotfile)
                print 'Plot: ',plotfile

sumfile.close()
print 'Enviroment fractions: ',envsumfile
