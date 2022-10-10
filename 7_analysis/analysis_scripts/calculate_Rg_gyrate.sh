#!/bin/bash

SPOP_type=MATHfree
gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

for i in $(seq 2 2 60)
do

cd ${SPOP_type}_runs/SPOP_${i}mer

$gmx gyrate -f prodrun_AAbackmapped.xtc -s prodrun_AAbackmapped.gro -o Rg_gyrate_AA.xvg <<EOF
Protein
EOF

cd ../..

done


