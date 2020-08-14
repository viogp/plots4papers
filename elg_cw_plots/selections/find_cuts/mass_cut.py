import os.path, sys
import numpy as np
from scipy.interpolate import interp1d

inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/selections/'

model = 'gp19/'

snap_list = [41,39]

ngal_vals = [-2.,-3.,-4.2]

# The order of this array needs to match that
# from the files generated with mass_cum.py
surveys = ['All','DEEP2','VVDS-DEEP','eBOSS-SGC','DESI'] 

for zsnap in snap_list:
    # Read the cumulative abundance for the stllar mass
    infile = inpath+model+'mass_cum_sn'+str(zsnap)+'.dat'
    if (not os.path.isfile(infile)):
        print('STOP: {} not found'.format(infile)) ; sys.exit()
    data = np.loadtxt(infile, unpack=True)
    mass = data[0,:] #; print( mass)

    # Write output header
    outfile = inpath+model+'ngal_mass_cuts_sn'+str(zsnap)+'.dat'
    ff = open(outfile,'w') ; print('Outfile: {}'.format(outfile))
    ff.write('# log(ngal), log(M*), Survey \n' )

    # Find mass cuts
    for ngal_val in ngal_vals:
        for ii, survey in enumerate(surveys):
            y = data[ii+1,:]

            if (max(y) >= ngal_val):
                f = interp1d(y,mass)
                cut = f(ngal_val)
            else:
                cut = -999.

            ff.write(' '.join(str(jj) \
                                  for jj in [ngal_val,cut,survey]))
            ff.write(' \n')
    ff.close()
