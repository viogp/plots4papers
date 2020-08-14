import sys
import numpy as np

model = 'gp19/'
path = '/cosma5/data/durham/violeta/lines/cosmicweb/env_files/'+model

infile = path+'env_fractions.txt'

sns = ['39','41']
nds = ['-2.0','-3.0','-4.2']
cuts = ['m','sfr','lo2']
nelgs = ['DEEP2','DESI','VVDS-DEEP','eBOSS-SGC']
lse = ['Voids','Sheets','Filaments','Knots']
cw_list = ['Vweb','Pweb']
###########################
dum = -999. ; eps = 1e-8

# Initialize the arrays
alls = np.zeros(shape=(len(cw_list),len(sns),len(nds),len(cuts),len(lse)))
alls.fill(dum)
elgs = np.zeros(shape=(len(cw_list),len(sns),len(nds),len(cuts),len(nelgs),len(lse)))
elgs.fill(dum)

with open(infile) as ff:
    for line in ff:
        if (not line.rstrip() or line[0]=='#'):
            continue

        # Data
        array1 = line.split(':')[1]
        array2 = array1.split('[')[1]
        array3 = array2.split(']')[0]
        array = array3.split()

        # Name
        name0 = line.split(':')[0]
        cw = name0.split('_')[0]
        icw = cw_list.index(cw) 
        name = name0.split('web_')[1]
        cut = name.split('cut')[0]
        ic = cuts.index(cut)
        nd = (name.split('_')[2]).split('nd')[1]
        ind = nds.index(nd)
        sn = (name.split('sn')[1]).split('.')[0]
        isn = sns.index(sn)
        survey = name.split('_')[1]
        if (survey == 'All'):
            alls[icw,isn,ind,ic,:] = array
        else:
            ielg = nelgs.index(survey)
            elgs[icw,isn,ind,ic,ielg,:] = array

# Output summary file with comparissons
envsumfile = path+'fracs_compared.txt'

sumfile = open(envsumfile,'w') 
sumfile.write('#### abs(Difference) (Percentage) ELG compared to the whole population #### \n') 
sumfile.close()
sumfile = open(envsumfile,'a') 

# Loop over arrays, compare mcut with mcut
for icw in range(len(cw_list)):
    sumfile.write('{} ###############\n'.format(cw_list[icw]))
    for ind in range(len(nds)):
        for ilse in range(len(lse)):
            for  ic in range(len(cuts)):
                for ielg in range(len(nelgs)):
                    if (nelgs[ielg] == 'DEEP2' or nelgs[ielg] == 'DESI'):
                        sn = '39' ; isn = sns.index(sn)
                    else:
                        sn = '41' ; isn = sns.index(sn)

                    val = alls[icw,isn,ind,ic,ilse]

                    if (elgs[icw,isn,ind,ic,ielg,ilse]>0. and val>eps):
                        per = abs(val-elgs[icw,isn,ind,ic,ielg,ilse])*100./val
                        diffs = abs(val-elgs[icw,isn,ind,ic,ielg,ilse])

                        sumfile.write('{:.2f} ({:.2f}%) {}; All_frac={:.3f} ### {}, {}cut, nd{}, sn{}\n '.format(diffs,per,nelgs[ielg],val,lse[ilse],cuts[ic],nds[ind],sns[isn]))
                    else:
                        sumfile.write('   not enough data\n ')

sumfile.write('\n ###################\n'.format(cw_list[icw]))

# Loop over arrays, compare ELGs with sfrcut
for icw in range(len(cw_list)):
    sumfile.write('{} ############### ELG vs sfrcut\n'.format(cw_list[icw]))
    for ind in range(len(nds)):
        for ilse in range(len(lse)):
            for  ic in range(len(cuts)):
                for ielg in range(len(nelgs)):
                    if (nelgs[ielg] == 'DEEP2' or nelgs[ielg] == 'DESI'):
                        sn = '39' ; isn = sns.index(sn)
                    else:
                        sn = '41' ; isn = sns.index(sn)

                    val = alls[icw,isn,ind,1,ilse]

                    if (elgs[icw,isn,ind,ic,ielg,ilse]>0. and val>eps):
                        per = abs(val-elgs[icw,isn,ind,ic,ielg,ilse])*100./val
                        diffs = abs(val-elgs[icw,isn,ind,ic,ielg,ilse])

                        sumfile.write('{:.2f} ({:.2f}%) {}; All_sfr_frac={:.3f} ### {}, {}cut, nd{}, sn{}\n '.format(diffs,per,nelgs[ielg],val,lse[ilse],cuts[ic],nds[ind],sns[isn]))
                    else:
                        sumfile.write('   not enough data\n ')


sumfile.write('\n ##### Frac(Sheets+filaments) #####\n')
for icw in range(len(cw_list)):
    sumfile.write('{} ###############\n'.format(cw_list[icw]))
    for ind in range(len(nds)):
        for ielg in range(len(nelgs)):
            if (nelgs[ielg] == 'DEEP2' or nelgs[ielg] == 'DESI'):
                sn = '39' ; isn = sns.index(sn)
            else:
                sn = '41' ; isn = sns.index(sn)

            for  ic in range(len(cuts)):
                addfs = 0.
                for ilse in range(len(lse)):
                    if (lse[ilse] in ['Sheets','Filaments']):
                        if (elgs[icw,isn,ind,ic,ielg,ilse]>0.):
                            addfs = addfs + elgs[icw,isn,ind,ic,ielg,ilse]
                sumfile.write('* {}, {}, {}cut, nd{}, sn{} \n'
                  .format(addfs,nelgs[ielg],cuts[ic],nds[ind],sns[isn]))                

sumfile.close()
print('Comparison of fractions: {}'.format(envsumfile))
