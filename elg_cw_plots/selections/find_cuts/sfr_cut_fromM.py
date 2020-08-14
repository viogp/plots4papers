import os.path, sys
import numpy as np
from scipy.interpolate import interp1d

inpath = '/cosma5/data/durham/violeta/lines/cosmicweb/selections/'

model = 'gp19/'

snap_list = [41,39]

mvals = [8.5,10.]

# The order of this array needs to match that
# from the files generated with mass_cum.py
elgs = ['DEEP2','VVDS-DEEP','eBOSS-SGC','DESI'] 

for zsnap in snap_list:
    # Read the number densities for given stellar masses cuts
    infile = inpath+model+'mass_cuts_sn'+str(zsnap)+'.dat'
    if (not os.path.isfile(infile)):
        print('STOP: {} not found'.format(infile)) ; sys.exit()
    cut_val, lngals1 = np.loadtxt(infile, usecols=(1,2), unpack=True)

    # Read the cumulative abundance for the SFR 
    infile = inpath+model+'sfr_cum_sn'+str(zsnap)+'.dat'
    if (not os.path.isfile(infile)):
        print('STOP: {} not found'.format(infile)) ; sys.exit()
    data = np.loadtxt(infile, unpack=True)
    sfr = data[0,:] #; print( sfr)
    allg = data[1,:]

    # Write output header
    outfile = inpath+model+'fromM_sfr_cuts_sn'+str(zsnap)+'.dat'
    ff = open(outfile,'w') ; print('Outfile: {}'.format(outfile))
    ff.write('# ')
    ff.write('  '.join([str(ii)+'='+elg for ii,elg in enumerate(elgs)]))
    ff.write(' \n')
    ff.write('# ELG survey, log10(M*) value, log(ngal), log(SFR_all), log(SFR_ELGs) \n' )

    # Find sfr cuts
    for mval in mvals:
        ind = np.where(cut_val == mval)
        lngals = lngals1[ind] # number densities

        for ii, elg in enumerate(elgs):
            y = data[ii+2,:]
            ngal_val = lngals[ii]

            finterp = interp1d(y,sfr) ; cut_elg = finterp(ngal_val)
            #print('{}, ngal={}, SFR cut={}'.format(elg,ngal_val,cut_elg))

            f = interp1d(allg,sfr) ; cut_all = f(ngal_val)

            ff.write(' '.join(str(jj) \
                                  for jj in [ii,mval,ngal_val,cut_all,cut_elg]))
            ff.write(' \n')
    ff.close()
