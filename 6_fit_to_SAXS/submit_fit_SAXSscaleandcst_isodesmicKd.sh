#!/bin/bash

python=/storage1/thomasen/software/miniconda3/bin/python3.7

SPOP_type=MATHfree
Ka=0.42

nohup $python Reweighting.py -i 500 -sr 501 -os 60 -s 2 -fsl -fk -c 5 10 20 30 40 -cp ../calc_data_${SPOP_type} -k $Ka --mcstepsize 0.1 -tinit 10.0 -tstop 0.1 -tdecr 0.3 &
