import sys,os.path
import subprocess
import numpy as np

Testing = False

propname = 'mh'

##########################################
model = 'gp19/'
path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
proppath = path+'selections/'+model+'ascii_files/'
volume = np.power(500.,3.)
##########################################

cut = 'elg'

elabels = ['Voids','Sheets','Filaments','Knots']

surveys1 = ['DEEP2','VVDS-DEEP']
surveys2 = ['DESI','eBOSS-SGC']
snaps = ['39','41']
zzs = ['0.99','0.83']
cw_list = ['Vweb','Pweb']

if Testing:
    surveys1 = ['DEEP2'] ; surveys2 = ['DESI']
    snaps = ['39']
    cw_list = ['Vweb'] 

##########################################
# Output percentage summary
envsumfile = path+'env_files/'+model+propname+'_elg.txt'
sumfile = open(envsumfile, 'w')
sumfile.write('# Median host halo mass ELGs, P10, P90, LSE, Survey, Snapshot, Vweb/Pweb \n')
sumfile.close()
sumfile = open(envsumfile, 'a')

# Loop over the different files
for cw in cw_list:
    epath = path+'env_files/'+model+cw+'/'
    plotroot = path+'plots/'+model+'environ/props/'+cw+'/elgs/'
    print('\n Plots: {}* \n'.format(plotroot))

    for iiz, sn in enumerate(snaps):
        surveys = [surveys1[iiz],surveys2[iiz]]

        for survey in surveys:
            end = cut+'cut_'+survey+'_sn'+sn+'.dat'
            efile = epath+end
            pfile = proppath+end

            # Check if files exist and has more than one line
            if (not os.path.isfile(efile) or not os.path.isfile(pfile)):
                if Testing: print('Jumping {} or {}'.format(efile,pfile))
                continue
            wcl_line = subprocess.check_output(["wc", "-l",efile])
            wcl = int(wcl_line.split()[0])
            pcl_line = subprocess.check_output(["wc", "-l",pfile])
            pcl = int(pcl_line.split()[0])
            if (wcl <= 1 or pcl <= 1):
                if Testing: print('Jumping {} or {}'.format(efile,pfile))
                continue

            # Read property file
            # xgal 0, ygal 1, zgal 2 (Mpc/h), vxgal 3,vygal 4,vzgal 5 (Km/s),
            # log10(massh) 6, log10(mass/Msun/h) 7, log10(sfr/Msun/h/Gyr) 8, 
            # log10(lum/h^-2 erg/s) 9, log10(lum_ext/h^-2 erg/s) 10,
            # type 11 (0= Centrals; 1,2= Satellites) 
            px, py, pz, mh = np.loadtxt(pfile, usecols=(0,1,2,6), unpack=True)

            # Write total to summary file
            prop = 10.**mh
            q1, median, q9 = np.percentile(prop, [10, 50, 90])
            lq1 = np.log10(q1)
            lme = np.log10(median)
            lq9 = np.log10(q9)
            sumfile.write('{:.2f} {:.2f} {:.2f} {} {} {} {}\n'.format(
                lme,lq1,lq9,'Total',survey,sn,cw)) 
            
            # Read the environmental file
            xx, yy, zz, fenv = np.loadtxt(efile,unpack=True)
            env = fenv.astype(int)

            # Chech that the coordinates have the same size
            if ((len(xx) != len(px)) or 
                (len(yy) != len(py)) or
                (len(zz) != len(pz))):
                print('[PROBLEM] Different lengths coordinates: {}\n {}\n'.
                      format(efile,pfile)) ; continue
            # Chech that the coordinates are ordered in the same way
            if ((not np.allclose(xx,px,atol=1e-08,equal_nan=True)) or 
                (not np.allclose(yy,py,atol=1e-08,equal_nan=True)) or
                (not np.allclose(zz,pz,atol=1e-08,equal_nan=True))):
                print('[PROBLEM] Files with different coordinates: {}\n {}\n'.
                      format(efile,pfile)) ; continue

            # Loop over type of environment
            for iie, ienv in enumerate(np.unique(env)):
                ind = np.where(env == ienv)
                if(np.shape(ind)[1]<2): continue
                prop = 10.**mh[ind]
                q1, median, q9 = np.percentile(prop, [10, 50, 90])
                lq1 = np.log10(q1)
                lme = np.log10(median)
                lq9 = np.log10(q9)
                sumfile.write('{:.2f} {:.2f} {:.2f} {}\n'.format(lme,lq1,lq9,elabels[iie])) 

print('Summary file: {}'.format(envsumfile))
sumfile.close()
