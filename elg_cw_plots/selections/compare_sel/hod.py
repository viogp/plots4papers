import sys, os.path
import subprocess
import numpy as np

model = 'gp19/'

sn_list = ['41','39']
surveys1 = ['VVDS-DEEP','DEEP2']
surveys2 = ['eBOSS-SGC','DESI']
nds = ['-2.0','-3.0','-4.2']
#cuts = ['m','sfr']
cuts = ['lo2'] 

verbose = False

# Paths
inpath = '/cosma6/data/dp004/dc-gonz3/lines/cosmicweb/'
ndpath = inpath+'selections/'
hodpath = inpath+'hod/'

# Read the information from the different surveys
for iiz, sn in enumerate(sn_list): 
    # Read the number of haloes
    #hmffile = hodpath+model+'hmf_sn'+sn+'_2vols.txt'
    hmffile = hodpath+model+'hmf_sn'+sn+'.txt' # Full volume
    mhist,mlow,mhigh,nhs = np.loadtxt(hmffile,unpack='True')

    surveys = ['All',surveys1[iiz],surveys2[iiz]]

    for cut in cuts:
        for survey in surveys:
            for nd in nds:
                infile = ndpath+model+'ascii_files/'+cut+\
                         'cut_'+survey+'_nd'+nd+'_sn'+sn+'.dat'

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

                massh,mass,x = np.loadtxt(infile,usecols=(6,7,11),unpack=True)
                gtype = x.astype(int)

                # Initialize array for hod
                hod = np.zeros(shape=(3,len(nhs))) ; hod.fill(-999.)

                for ii,nh in enumerate(nhs):
                    # All
                    ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]))
                    numa = np.shape(ind)[1]
                    if (numa >= 1 and nh > 0):
                        hod[0,ii] = float(numa)/nh 

                    # Centrals
                    ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]) &
                                   (gtype == 0))
                    numc = np.shape(ind)[1]
                    if (numc >= 1 and nh > 0):
                        hod[1,ii] = float(numc)/nh

                    # Satellites
                    ind = np.where((massh>=mlow[ii]) & (massh<mhigh[ii]) &
                                   (gtype > 0))
                    nums = np.shape(ind)[1] 
                    if (nums >= 1 and nh > 0):
                        hod[2,ii] = float(nums)/nh 

                # Take the log of the HODs:
                ind = np.where(hod>0.)
                hod[ind] = np.log10(hod[ind])

                # Output file
                hfile = hodpath+model+cut+\
                         'cut_'+survey+'_nd'+nd+'_sn'+sn+'.dat'
                tofile = np.column_stack((mhist,hod[0,:],hod[1,:],hod[2,:]))
                with open(hfile, 'w') as outf:
                    outf.write('# log10(Mh/Msun/h)_midpoint, log10<Nall>, log10<Ncentrals>, log10<Nsatellites>')
                    np.savetxt(outf,tofile,fmt=('%.5f'))

                print('Output: {}'.format(hfile))


                
