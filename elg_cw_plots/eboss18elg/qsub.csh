#! /bin/tcsh -f

set jobname    = modelplots
set logname = /cosma5/data/durham/$user/Junk/$jobname.%A.%a.log

set script = run.csh

cat << EOF - ${script} | sbatch 
#!/bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${jobname}
#SBATCH -o ${logname}
#SBATCH -p cordelia
#SBATCH -A durham
#SBATCH -t 3:00:00
#

# Run script follows
EOF


