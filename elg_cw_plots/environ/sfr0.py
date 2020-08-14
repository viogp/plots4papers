import sys,os.path
import subprocess
import numpy as np

Testing = False

propname = 'sfr'

##########################################
model = 'gp19/'
path = '/cosma5/data/durham/violeta/lines/cosmicweb/'
proppath = path+'selections/'+model+'ascii_files/'
hodpath = path+'hod/'
volume = np.power(500.,3.)
##########################################

cut = 'elg'

elabels = ['Voids','Sheets','Filaments','Knots']

cols = ['k','greenyellow','limegreen','forestgreen','darkolivegreen']
lwidth = [4,2,2,2,2] 

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

# Loop over the different files
for cw in cw_list:
    epath = path+'env_files/'+model+cw+'/'
    plotroot = path+'plots/'+model+'environ/props/'+cw+'/elgs/'

    for iiz, sn in enumerate(snaps):
        surveys = [surveys1[iiz],surveys2[iiz]]

        for survey in surveys:
            # Files to read
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
            px, py, pz, mass, sfr = np.loadtxt(pfile, usecols=(0,1,2,7,8), unpack=True)

            # Total values
            tot = len(sfr)
            tot0 = len(sfr[np.where(sfr<=0.)])
            print('Total 0s: {}% ({}) ; {} {} {}'.format(100.*tot0/tot,tot0,cw,survey,sn))

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
                emass = mass[ind]
                esfr = sfr[ind]

                tot = len(esfr)
                tot0 = len(esfr[np.where(esfr<=0.)])
                print('{} 0s: {}% ({})'.format(elabels[iie],100.*tot0/tot,tot0))
