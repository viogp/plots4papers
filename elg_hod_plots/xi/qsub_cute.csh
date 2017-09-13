#! /bin/tcsh

set testing = false
echo "Running CUTE with testing =" $testing

set zspace  = false #true #true

set name    = cute
set logpath = /gpfs/data/$user/Junk/$name

# Test file: /gpfs/data/violeta/CUTE/CUTE_box/test/params.txt
set file_list = (`awk '{print $1}' in.txt`)

#set parfil = 'params.txt'
@ i = 1
foreach file ($file_list)
    if ($zspace == true) then
	set parfil = 'params/params'${i}'_z_gp15newmg.txt'
	set ofile = (`echo $file | sed 's#_4cute_z.dat#_CUTExi.dat#'`)
	set outfil = (`echo $ofile | sed 's#/O#/z/O#'`)
    else
	set parfil = 'params/params'${i}'_gp15newmg.txt'
	set ofile = (`echo $file | sed 's#_4cute.dat#_CUTExi.dat#'`)
	set outfil = (`echo $ofile | sed 's#/O#/r/O#'`)
    endif
    
    cp example.params.txt $parfil
    ./replace_variable.csh $parfil data_filename $file
    ./replace_variable.csh $parfil output_filename $outfil	    

    if ($testing == true) then
    	./run_cute.csh $parfil
    else
    	echo "Submitting job to Cordelia, param_file =" $parfil
    	time bsub -P durham -n 1 -q cordelia -J "$name" -o $logpath.%J.%I run_cute.csh $parfil
    endif

    @ i += 1
end

echo "The end"
