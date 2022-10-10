#!/bin/bash

for i in $(seq 2 2 12)
do

cd SPOP_${i}mer
cp ../prodrun_grompp.sh .
qsub prodrun_grompp.sh
cd ..

done
