import os.path, sys
import numpy as np
from scipy.interpolate import interp1d

inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/selections/'

model = 'gp19/'
#model = 'gp19.font/'

snap_list = [41,39]

values = [8.5,10.]

# The order of this array needs to match that 
# from the files generated with mass_cum.py
elgs = ['DEEP2','VVDS-DEEP','eBOSS-SGC','DESI'] 

for zsnap in snap_list:
    infile = inpath+model+'mass_cum_sn'+str(zsnap)+'.dat'
    if (not os.path.isfile(infile)):
        print('STOP: {} not found'.format(infile)) ; sys.exit()
    
    # Read cumulative abundance
    data = np.loadtxt(infile, unpack=True)
    mass = data[0,:] #; print( mass)
    allg = data[1,:]

    # Write output header
    outfile = inpath+model+'mass_cuts_sn'+str(zsnap)+'.dat'
    ff = open(outfile,'w') ; print('Outfile: {}'.format(outfile))
    ff.write('# ')
    ff.write('  '.join([str(ii)+'='+elg for ii,elg in enumerate(elgs)]))
    ff.write(' \n')
    ff.write('# ELG survey, log10(M*) cut value, log(ELG ngal), Corresponding log10(M*) for the whole sample \n' )

    # Find mass cuts
    for value in values:
        for ii, elg in enumerate(elgs):
            y = data[ii+2,:]
            f = interp1d(mass,y)
            ngal_val = f(value)

            f = interp1d(allg,mass)
            mcut = f(ngal_val)

            ff.write(' '.join(str(jj) \
                                  for jj in [ii,value,ngal_val,mcut]))
            ff.write(' \n')
    ff.close()
