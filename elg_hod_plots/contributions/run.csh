#!/bin/tcsh -f

set nom = $1
set model = $2

set exec =  ${nom}'_contributions.py'

echo $exec $model
python $exec $model


