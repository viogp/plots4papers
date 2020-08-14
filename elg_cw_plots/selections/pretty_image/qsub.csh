#! /bin/tcsh -f

set name    = selections
set logname = /cosma5/data/durham/$user/Junk/$name.%A.%a.log

set script = run.csh

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

