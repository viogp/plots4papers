import sys, os.path
import subprocess
import numpy as np
from Cosmology import *

space = 'r' # 'r' or 'z'

model = 'gp19'
print('WARNING : check that the cosmology corresponds to the model')
# Read the redshift and cosmology
set_cosmology(omega0=0.27,omegab=0.05,lambda0=0.73, \
              h0=0.70, universe="Flat",include_radiation=False)

sn_list = ['41','39']
zz_list = [0.83,0.99]
surveys1 = ['VVDS-DEEP','DEEP2']
surveys2 = ['eBOSS-SGC','DESI']
cuts = ['shuffled_elg']

verbose = True

# Testing
testing = False
if testing:
    verbose = True
    sn_list = ['39'] ; zz_list = [0.99]
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI']
    nds = ['-2.0']
#------------------

# Paths
inpath = '/cosma6/data/dp004/dc-gonz3/lines/cosmicweb/selections/'+\
         model

# File with all the file names to input CUTE
info_file = model+'_inputs4cute_shuffled_elg_env_'+space+'.txt'
info = open(info_file,'w')

# Loop over all the files
for iis,sn in enumerate(sn_list):
    zz = zz_list[iis]
    surveys = [surveys1[iis],surveys2[iis]]

    for cut in cuts:
        for survey in surveys:
            infile = inpath+'/ascii_files/'+\
                     cut+'cut_'+survey+'_sn'+sn+'.dat'

            # Check if the file exists
            if (not os.path.isfile(infile)):
                if verbose:
                    print('WARNING: {} not found'.format(infile))
                continue
            # Jump files with only the header (1 line)
            wcl_line = subprocess.check_output(["wc","-l",infile])
            wcl = int(wcl_line.split()[0])
            if (wcl <= 1):
                if verbose:
                    print('WARNING: {} has too few lines'.format(infile))
                continue

            print('Reading {}'.format(infile))

            # Read the ascii files with the number density selections
            if (space == 'r'): # r-space
                x,y,z = np.loadtxt(infile,usecols=(0,1,2),unpack=True)
            else: # z-space
                x1,y,z,vx = np.loadtxt(infile,
                                       usecols=(0,1,2,3),unpack=True)
                x  = x1 + vx*(1.+zz)/H(zz)

            # Write the input file for CUTE
            outfile = inpath+'/iz'+sn+'/'+space+'/'+\
                     cut+'cut_'+survey+'_sn'+sn+\
                     '_4cute_'+space+'.dat'

            tofile = np.column_stack((x,y,z))
            with open(outfile,'w') as outf:
                np.savetxt(outf, tofile, fmt='%.5f')

            # Write file name in info_file
            info.write(outfile+' \n')
info.close()
print(' - Info file: {}'.format(info_file))
