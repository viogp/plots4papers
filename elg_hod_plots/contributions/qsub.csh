#! /bin/tcsh -f

set noms = (deep2 vvdsdeep vvdswide eboss desi)
#set noms = ( eboss)

#set model = 'MillGas/gp14/' 
set model = 'MillGas/gp15newmg.anders/' 

foreach nom ($noms)
    set logpath = /gpfs/data/$user/Junk/contributions_$nom  

    time bsub -P durham -n 1 -q cordelia -J "$nom" -o $logpath.%J.%I run.csh $nom $model
	    
    # Testing
    #./run.csh $nom $model
end

echo 'All submited'
