import sys, os.path
import subprocess
import numpy as np
from Cosmology import *

Testing = False

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
cuts = ['elg']

elabels = ['Voids','Sheets','Filaments','Knots'] 
cw_list = ['Vweb','Pweb']

# Testing
verbose = False
if Testing:
    verbose = True
    sn_list = ['39'] ; zz_list = [0.99]
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI']
    cw_list = ['Vweb'] 

#------------------

# Paths
path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
inpath = path+'selections/'+model

# File with all the file names to input CUTE
info_file = model+'_inputs4cute_elg_env_'+space+'.txt'
info = open(info_file,'w')

# Loop over all the files
for cw in cw_list:
    for iis,sn in enumerate(sn_list):
        zz = zz_list[iis]
        surveys = [surveys1[iis],surveys2[iis]]
    
        for cut in cuts:
            for survey in surveys:
                end = cut+'cut_'+survey+'_sn'+sn+'.dat' 
                infile = inpath+'/ascii_files/'+end
                efile = path+'env_files/'+model+'/'+cw+'/'+end

                # Check if files exist
                if (not os.path.isfile(infile) or not os.path.isfile(efile)):
                    if verbose:
                        print('WARNING: file not found {} or {}'.format(infile,efile))
                    continue
                # Jump files with only the header (1 line)
                wcl_line = subprocess.check_output(["wc","-l",infile])
                wcl = int(wcl_line.split()[0])
                ecl_line = subprocess.check_output(["wc","-l",efile])
                ecl = int(wcl_line.split()[0])
                if (wcl <= 1 or ecl <= 1):
                    if verbose:
                        print('WARNING: Too few lines in {} or {}'.format(infile,efile))
                    continue
    
                # Read the ascii files with the number density selections
                if (verbose) : print('Reading {}, {}'.format(infile,efile))
                px,py,pz,vx = np.loadtxt(infile,usecols=(0,1,2,3),unpack=True)
    
                # Read the environmental file                                            
                xx, yy, zz, fenv = np.loadtxt(efile,unpack=True)                         
                env = fenv.astype(int)                                                   
                                                                                         
                # Chech that the coordinates have the same size                          
                if ((len(xx) != len(px)) or                                              
                    (len(yy) != len(py)) or                                              
                    (len(zz) != len(pz))):                                               
                    print('[PROBLEM] Different lengths coordinates: {}\n {}\n'.          
                          format(efile,infile)) ; continue                                
                # Check that the coordinates are ordered in the same way                 
                if ((not np.allclose(xx,px,atol=1e-08,equal_nan=True)) or                
                    (not np.allclose(yy,py,atol=1e-08,equal_nan=True)) or                
                    (not np.allclose(zz,pz,atol=1e-08,equal_nan=True))):                 
                    print('[PROBLEM] Files with different coordinates: {}\n {}\n'.       
                          format(efile,infile)) ; continue   
    
                # Loop over type of environment
                for iie, ienv in enumerate(np.unique(env)):
                    ind = np.where(env == ienv)
                    if (np.shape(ind)[1]<2): 
                        print('Not enough {}'.format(elabels[iie]))
                        continue

                    y = py[ind] ; z = pz[ind]
                    if (space == 'r'): # r-space
                        x = px[ind]
                    else:
                        x  = px[ind] + vx[ind]*(1.+zz)/H(zz)
    
                    # Write the input file for CUTE
                    outfile = inpath+'/iz'+sn+'/'+space+'/'+\
                              elabels[iie]+'_'+cw+'_'+cut+'cut_'+\
                              survey+'_sn'+sn+\
                              '_4cute_'+space+'.dat'
                    if (verbose): print('Output: {}'.format(outfile))

                    #arr_ienv = np.zeros(shape=len(x)) ; arr_ienv.fill(ienv)
                    tofile = np.column_stack((x,y,z))
                    with open(outfile,'w') as outf:
                        np.savetxt(outf, tofile, fmt='%.5f')
    
                    # Write file name in info_file
                    info.write(outfile+' \n')

                # Flush arrays
                px.fill(-999.) ; py.fill(-999.) ; pz.fill(-999.)
                vx.fill(-999.) ; fenv.fill(-999.)
                xx.fill(-999.) ; yy.fill(-999.) ; zz.fill(-999.)
                x.fill(-999.) ; y.fill(-999.) ; z.fill(-999.)

info.close()
print(' - Info file: {}'.format(info_file))
