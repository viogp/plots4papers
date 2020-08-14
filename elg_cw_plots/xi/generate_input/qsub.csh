#! /bin/tcsh -f

set jobname = xi
set logname = /cosma5/data/durham/$user/Junk/$jobname.%A.%a.log

#time bsub -P durham -n 1 -q cordelia -J "$name" -o $logpath.%J.%I run.csh
	    
# Testing
#./run.csh 

set script = run.csh

cat << EOF - ${script} | sbatch 
#!/bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${jobname}
#SBATCH -o ${logname}
#SBATCH -p cordelia
#SBATCH -A durham
#SBATCH -t 72:00:00
#

# Run script follows
EOF


