import sys, os.path
import subprocess
import numpy as np
from stats import percentiles

model = 'gp19/'

sn_list = ['41','39']
surveys1 = ['VVDS-DEEP','DEEP2']  
surveys2 = ['eBOSS-SGC','DESI'] 
nds = ['-2.0','-3.0','-4.2']
cuts = ['m','sfr','lo2'] 

verbose = False

# Paths
inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/'
ndpath = inpath+'selections/'

# Output file
ofile = ndpath+model+'medians.dat'
with open(ofile,'w') as of:
    of.write('#sample median(log10(massh),log10(mass/Msun/h), log10(sfr/Msun/h/Gyr), log10(lum_ext/h^-2 erg/s)) \n')

# Read the information from the different files
for iiz, sn in enumerate(sn_list):
    surveys = ['All',surveys1[iiz],surveys2[iiz]] 

    for iin,nd in enumerate(nds):
        for iic,cut in enumerate(cuts):
            for iis,survey in enumerate(surveys):
                mtext = cut+'cut_'+survey+'_nd'+nd+'_sn'+sn
                infile = ndpath+model+'ascii_files/'+mtext+'.dat'

                # Check the existance of the files
                if (not os.path.isfile(infile)):
                    if verbose:
                        print('WARNING: {} not found'.format(infile))
                    continue
                # Jump files with only a one line header
                wcl_line = subprocess.check_output(["wc", "-l",infile])
                wcl = int(wcl_line.split()[0])
                if (wcl <= 1):
                    if verbose:
                        print('WARNING: {} has too few lines'.format(infile))
                    continue

                #massh,mass,sfr,lum_ext
                data = np.loadtxt(infile,usecols=(6,7,8,10),unpack=True)
                val0 = np.shape(data)[0] ; val1 = np.shape(data)[1] 

                if (cut != 'lo2'):
                    # Convert luminosity to log10(luminosity)+40.
                    lum_ext = data[val0-1,:] 
                    dd = np.zeros(shape=val1) ; dd.fill(-999.)
                    ind = np.where(lum_ext > 0.) 
                    dd[ind] = np.log10(lum_ext[ind]) + 40.
                    data[val0-1,:] = dd

                meds = np.zeros(shape=(val0)) ; meds.fill(-999.)

                for ii in range(val0):
                    dd = data[ii,:] 
                    ind = np.where(dd>-999.) 
                    median = percentiles(0.5,dd[ind]) 
                    meds[ii] = median
                    mtext = mtext+' '+str(median)

                # Write output file
                with open(ofile, 'a') as of:
                    of.write(mtext+' \n')

print('File with medians: {}'.format(ofile))
