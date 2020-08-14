#! /bin/tcsh -f

set Testing = False

#set noms = (mass mh sfr ssfr lo2)
set noms = (mass)

set logpath = /cosma5/data/durham/$user/Junk

foreach nom ($noms)
    set logname = ${logpath}/props_$nom.%A.%a.log
    set job_file = ${logpath}/${nom}.job
    set exec =  ${nom}'_env_dis.py'

    if ($Testing == 'True') then
	echo 'Testing'
	echo $exec
	python $exec
    else
	cat > $job_file <<EOF
#! /bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${nom}
#SBATCH -o ${logname}
#SBATCH -p cordelia
#SBATCH -A durham
#SBATCH -t 48:00:00
#
echo $exec
python3 $exec
#
EOF

	sbatch $job_file
	rm $job_file
    endif
end

echo 'All submited'
