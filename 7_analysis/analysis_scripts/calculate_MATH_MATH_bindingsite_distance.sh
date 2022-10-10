#!/bin/bash

SPOP_type=$1
oligomersize=$2

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

dir=${SPOP_type}_runs/SPOP_${oligomersize}mer
#rm -r ${dir}/MATH_MATH_bindingsite_distance
mkdir ${dir}/MATH_MATH_bindingsite_distance

last_subunit=$(( ${oligomersize}-1 ))

atoms_per_subunit=5241

ARG45_start=701
ARG45_end=724
TYR62_start=974
TYR62_end=994
SER94_start=1496
SER94_end=1506
TYR98_start=1558
TYR98_end=1578
LYS104_start=1663
PHE108_end=1747

#Loop over subunits in SPOP (except for last subunit)
for subunit in $(seq 1 $last_subunit)
do

MATHright_start=$(( ${MATH_start} + ${atoms_per_subunit} ))

r_ARG45_start=$(( ${ARG45_start} + ${atoms_per_subunit} ))
r_ARG45_end=$(( ${ARG45_end} + ${atoms_per_subunit} ))
r_TYR62_start=$(( ${TYR62_start} + ${atoms_per_subunit} ))
r_TYR62_end=$(( ${TYR62_end} + ${atoms_per_subunit} ))
r_SER94_start=$(( ${SER94_start} + ${atoms_per_subunit} ))
r_SER94_end=$(( ${SER94_end} + ${atoms_per_subunit} ))
r_TYR98_start=$(( ${TYR98_start} + ${atoms_per_subunit} ))
r_TYR98_end=$(( ${TYR98_end} + ${atoms_per_subunit} ))
r_LYS104_start=$(( ${LYS104_start} + ${atoms_per_subunit} ))
r_PHE108_end=$(( ${PHE108_end} + ${atoms_per_subunit} ))


#Make index file for subunit MATH and right neighbouring MATH
$gmx make_ndx -f ${dir}/prodrun_AAbackmapped.gro -o ${dir}/MATH_MATH_bindingsite_distance/MATH_${subunit}_MATHright.ndx <<EOF
a ${ARG45_start}-${ARG45_end} | a ${TYR62_start}-${TYR62_end} | a ${SER94_start}-${SER94_end} | a ${TYR98_start}-${TYR98_end} | a ${LYS104_start}-${PHE108_end} 
a ${r_ARG45_start}-${r_ARG45_end} | a ${r_TYR62_start}-${r_TYR62_end} | a ${r_SER94_start}-${r_SER94_end} | a ${r_TYR98_start}-${r_TYR98_end} | a ${r_LYS104_start}-${r_PHE108_end}
q
EOF

#Calculate MATH and right neighbouring MATH contacts for subunit
$gmx distance -f ${dir}/prodrun_AAbackmapped.xtc -s ${dir}/prodrun_AAbackmapped.gro -n ${dir}/MATH_MATH_bindingsite_distance/MATH_${subunit}_MATHright.ndx -oall ${dir}/MATH_MATH_bindingsite_distance/comdist_MATH_${subunit}_MATHright.xvg -select "com of group 10 plus com of group 11"

${dir}/MATH_MATH_bindingsite_distance/MATH_${subunit}_MATHright.ndx

ARG45_start=$(( ${ARG45_start} + ${atoms_per_subunit} ))
ARG45_end=$(( ${ARG45_end} + ${atoms_per_subunit} ))
TYR62_start=$(( ${TYR62_start} + ${atoms_per_subunit} ))
TYR62_end=$(( ${TYR62_end} + ${atoms_per_subunit} ))
SER94_start=$(( ${SER94_start} + ${atoms_per_subunit} ))
SER94_end=$(( ${SER94_end} + ${atoms_per_subunit} ))
TYR98_start=$(( ${TYR98_start} + ${atoms_per_subunit} ))
TYR98_end=$(( ${TYR98_end} + ${atoms_per_subunit} ))
LYS104_start=$(( ${LYS104_start} + ${atoms_per_subunit} ))
PHE108_end=$(( ${PHE108_end} + ${atoms_per_subunit} ))


done





