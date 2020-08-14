import sys, os.path
import subprocess
import numpy as np
from Corrfunc.theory.xi import xi
from Corrfunc.theory.DD import DD
import matplotlib ; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_style
plt.style.use(mpl_style.style1)

def _RR(binfile, lbox, N1, N2=None):
    if N2 is None: N2=N1 #Autocorrelation Function
    vshell = (binfile[1:]**3-binfile[:-1]**3)*4/3.*np.pi
    
    return N1*N2*vshell/float(lbox**3)

Testing = False
verbose = True

space = 'r' #'r' or 'z'  
model = 'gp19'

# Paths
path = '/cosma5/data/durham/violeta/lines/cosmicweb/' 
inpath = path+'selections/'+model
#plotpath = path+'plots/'+model+'selections/corsscorr/' 

# Selections
sn_list = ['39','41']
zz_list = ['0.99','0.83']
surveys1 = ['DEEP2','VVDS-DEEP'] 
surveys2 = ['DESI','eBOSS-SGC']
elabels = ['Voids','Sheets','Filaments','Knots']
cw_list = ['Vweb','Pweb'] 
cut = 'elg'

# Parameters for the monopole calculation
boxsize = 500.
NBin = 30
log_rmin = -2.
log_rmax = np.log10(65.)
nthreads = 16 # 16 cores per node in Cosma5

if Testing: 
    verbose = True
    log_rmax=1.
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI'] 
    sn_list = ['39'] ; zz_list = ['0.99']
    cw_list = ['Vweb']


binfile = np.logspace(log_rmin,log_rmax,NBin+1,True)
bin_middle = (binfile[1:] + binfile[:-1])/2.

# Loop over all the files
for cw in cw_list:
    for iis, sn in enumerate(sn_list):
        zz = zz_list[iis]
        surveys = [surveys1[iis],surveys2[iis]]  

        for survey in surveys: 
            end = cut+'cut_'+survey+'_sn'+sn+'.dat' 
            infile = inpath+'/ascii_files/'+end 
            # xgal 0,ygal 1,zgal 2 (Mpc/h), vxgal 3,vygal 4,vzgal 5 (Km/s), 
            # log10(massh) 6,log10(mass/Msun/h) 7, log10(sfr/Msun/h/Gyr) 8,
            # log10(lum/h^-2 erg/s) 9,log10(lum_ext/h^-2 erg/s) 10,
            # type (0= Centrals; 1,2= Satellites) 11
            efile = path+'env_files/'+model+'/'+cw+'/'+end
            # x,y,z (Mpc/h), environments: 0 Void; 1 sheet; 2 filament; 3 Knots

            # Check if files exist 
            if (not os.path.isfile(infile) or not os.path.isfile(efile)):
                if verbose:
                    print('WARNING: {} or {} not found'.format(infile,efile)) 
                continue

            # Jump files with only the header (1 line) 
            wcl_line = subprocess.check_output(["wc","-l",infile]) 
            wcl = int(wcl_line.split()[0]) 
            ecl_line = subprocess.check_output(["wc","-l",efile]) 
            ecl = int(wcl_line.split()[0]) 
            if (wcl <= 1 or ecl <= 1):  
                if verbose:
                    print('WARNING: {} or {} has too few lines'.format(infile,efile))
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
                      format(infile,efile)) ; continue 
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

                # 3-D pair-counts for the correlation function 
                DD1 = DD(False, nthreads, binfile,
                         x, y, z, 
                         X2=px, Y2=py, Z2=pz,boxsize=boxsize)
                ddc = DD1['npairs']
                
                # Get the cross-correlation functions
                cf = ddc/_RR(binfile,boxsize,len(px),len(x)) - 1

                # For the errors, the subsample DD is needed
                subxi = xi(boxsize, nthreads, binfile, x, y, z)
                dd = subxi['npairs']

                err = (1+cf)/np.sqrt(dd)
            
                # Write output
                outfile = inpath+'/iz'+sn+'/'+space+'/ccf/'+\
                              elabels[iie]+'_'+cw+'_'+cut+'cut_'+\
                              survey+'_sn'+sn+\
                              '_ccf_'+space+'.dat'

                tofile = np.column_stack((bin_middle,cf,err,dd,ddc))
                with open(outfile,'w') as outf:
                    outf.write('# r/h-1Mpc Cross-CF (1+CCF)/sqrt(DD) DD DDtotal \n')
                    np.savetxt(outf, tofile, fmt=['%.5f','%.5f','%.5f','%d', '%d'])

                print('Output: {}'.format(outfile))
