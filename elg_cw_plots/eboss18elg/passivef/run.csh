#!/bin/tcsh -f

set exec = passive_sn61.py
#set exec = gaea_sn61.py
#set exec = passive_sn39.py
#set exec = passive_sn41.py

python $exec

echo 'The end'

