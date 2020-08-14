#! /bin/tcsh -f

set name    = hods
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
#SBATCH -t 12:00:00
#

# Run script follows
EOF

