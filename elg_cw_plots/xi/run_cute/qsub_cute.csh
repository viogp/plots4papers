#! /bin/tcsh

set zspace  = false
set model   = gp19
#set model   = gp19.font
#set model   = gp19.starvation

echo "Running " $model ", with z-space=" $zspace

set testing = false
echo "Running CUTE with testing =" $testing

set jobname    = cute
set logname = /cosma6/data/dp004/dc-gonz3/Junk/$jobname.%A.%a.log

# Test file: /cosma5/data/durham/violeta/CUTE/CUTE_box/test/params.txt
set file_list = (`awk '{print $1}' in.txt`)

set script = run_cute.csh

#set parfil = 'params.txt'
@ i = 1
foreach file ($file_list)
    if ($zspace == true) then
	set parfil = 'params/params'${i}'_z_'${model}'.txt'
	set outfil = (`echo $file | sed 's#_4cute_#_CUTExi_#'`)
	#set outfil = (`echo $ofile | sed 's#/O#/z/O#'`)
    else
	set parfil = 'params/params'${i}'_'${model}'.txt'
	set outfil = (`echo $file | sed 's#_4cute#_CUTExi#'`)
	#set outfil = (`echo $ofile | sed 's#/O#/r/O#'`)
    endif

    cp example.params.txt $parfil
    ./replace_variable.csh $parfil data_filename $file
    ./replace_variable.csh $parfil output_filename $outfil	    

    if ($testing == true) then
    	./run_test.csh $parfil
    else
    	#time bsub -P durham -n 1 -q cordelia -J "$name" -o $logpath.%J.%I run_cute.csh $parfil

	cat << EOF - ${script} | sbatch
#!/bin/tcsh -ef
#
#SBATCH --ntasks 1
#SBATCH -J ${jobname}
#SBATCH -o ${logname}
#SBATCH -p cosma6
#SBATCH -A dp004
#SBATCH -t 2:00:00
#
# Set parameters
set parfil = ${parfil}

# Run script follows 
EOF
    endif

    @ i += 1
end

echo "The end"
