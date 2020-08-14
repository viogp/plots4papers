#! /bin/tcsh -f

set name    = sel_comp
set logname = /cosma5/data/durham/$user/Junk/$name.%A.%a.log

set script = run.csh

#time bsub -P durham -n 1 -q cordelia -J "$name" -o $logpath.%J.%I run.csh 

cat << EOF - ${script} | sbatch 
#!/bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${name}
#SBATCH -o ${logname}
#SBATCH -p cordelia
#SBATCH -A durham
#SBATCH -t 24:00:00
#

# Run script follows
EOF

