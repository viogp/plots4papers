#! /bin/tcsh -f

#set exec = 'fenv.py'
#set exec = 'fsat.py'
set exec = 'lo2_lf.py'
#set exec = 'hod.py'
#set exec = 'sfr_m.py'
#set exec = 'xi_bias.py'
#set exec = 'mh.py'

set Testing = False

set nom = 'elgs'
set logpath = /cosma5/data/durham/$user/Junk
set logname = ${logpath}/props_$nom.%A.%a.log
set job_file = ${logpath}/${nom}.job

if ($Testing == 'True') then
    echo 'Testing'
    echo $exec
    python3 $exec
else
    cat > $job_file <<EOF
#! /bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${nom}
#SBATCH -o ${logname}
#SBATCH -p cordelia
#SBATCH -A durham
#SBATCH -t 12:00:00
#
echo $exec
python3 $exec
#
EOF

    sbatch $job_file
    rm $job_file
endif

echo 'All submited'
