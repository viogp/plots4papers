#! /bin/tcsh -f

#set exec = 'plot_crosscorr.py'
set exec = 'crosscorr.py'

set Testing = False

set nom = 'elgs'
set logpath = /cosma5/data/durham/$user/Junk
set logname = ${logpath}/crosscor_$nom.%A.%a.log
set job_file = ${logpath}/${nom}.job

if ($Testing == 'True') then
    echo 'Testing'
    echo $exec
    python3 $exec
else
    cat > $job_file <<EOF
#! /bin/tcsh -ef
#
#SBATCH --ntasks 1 --cpus-per-task 16
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
