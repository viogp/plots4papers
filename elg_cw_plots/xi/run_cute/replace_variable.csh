#!/bin/tcsh -f

# USAGE: replace_variable.csh $inputs_file $name $value

if( $#argv != 3) then
    echo ERROR in replace_variable.csh: USAGE replace_variable.csh file name value
    exit 1
endif

set file = $1
set name = $2
set value = $3

# create temporary file
set tempfile = replace_variable_temp$$

# 1st create new file with line containing named variable removed
# the awk command tests whether the 1st string on the line exactly matches $name
set ncute = ${name}'='
set fil1 = (`awk -v var=$ncute '{if($1 == var) print $2}' $file`)

sed -i "s#${fil1}#${value}#" ${file}


