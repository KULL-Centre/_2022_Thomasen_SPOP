#!/bin/bash

python=/storage1/thomasen/software/miniconda3/bin/python3.7

SPOP_type=MATHfree
Ka=0.42

nohup $python Reweighting.py -i 1 -sr 2 -os 60 -s 2 -fsl -c 5 10 20 30 40 -k $Ka -cp ../calc_data_${SPOP_type} &
