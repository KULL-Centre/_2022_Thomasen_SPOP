#!/bin/bash

for i in $(seq 2 2 12)
do

cd SPOP_${i}mer
cp ../relax_grompp.sh .
qsub relax_grompp.sh
cd ..

done
