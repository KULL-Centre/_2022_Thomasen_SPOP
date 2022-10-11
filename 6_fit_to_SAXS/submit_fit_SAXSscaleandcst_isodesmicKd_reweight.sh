#!/bin/bash

python=/storage1/thomasen/software/miniconda3/bin/python3.7

SPOP_type=MATHfree
philimit=$1
Ka=0.42

nohup $python Reweighting.py -i 1000 -sr 0 -os 60 -s 2 -cp ../calc_data_${SPOP_type} -fsl -fk -c 5 10 20 30 40 -k $Ka --theta 100 --thetafinal 0.0000000000000001 --thetadecrease 0.02 --philimit $philimit --mcstepsize 0.1 -tinit 10.0 -tstop 0.1 -tdecr 0.3 &
