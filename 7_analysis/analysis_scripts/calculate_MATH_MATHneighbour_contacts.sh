#!/bin/bash

SPOP_type=$1
oligomersize=$2

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

dir=${SPOP_type}_runs/SPOP_${oligomersize}mer
#rm -r ${dir}/MATH_MATHneighbour_contacts
mkdir ${dir}/MATH_MATHneighbour_contacts

last_subunit=$(( ${oligomersize}-1 ))

atoms_per_subunit=5241
atoms_per_MATH=2245 #set up to (and with) 138VAL

MATH_start=1

#Loop over subunits in SPOP (except for last subunit)
for subunit in $(seq 1 $last_subunit)
do

#Set MATH domain boundaries based on start of current subunit MATH domain
MATH_end=$(( ${MATH_start} + ${atoms_per_MATH} ))
MATHright_start=$(( ${MATH_start} + ${atoms_per_subunit} ))
MATHright_end=$(( ${MATHright_start} + ${atoms_per_MATH} ))

echo MATH subunit: ${subunit}, MATHstart: ${MATH_start}, MATHend: ${MATH_end}, MATHright_start: ${MATHright_start}, MATHright_end: ${MATHright_end}

#Make index file for subunit MATH and right neighbouring MATH
$gmx make_ndx -f ${dir}/prodrun_AAbackmapped.gro -o ${dir}/MATH_MATHneighbour_contacts/MATH_${subunit}_MATHright.ndx <<EOF
a ${MATH_start}-${MATH_end}
a ${MATHright_start}-${MATHright_end}
q
EOF

#Calculate MATH and right neighbouring MATH contacts for subunit
$gmx mindist -f ${dir}/prodrun_AAbackmapped.xtc -s ${dir}/prodrun_AAbackmapped.gro -n ${dir}/MATH_MATHneighbour_contacts/MATH_${subunit}_MATHright.ndx -od ${dir}/MATH_MATHneighbour_contacts/mindist_MATH_${subunit}_MATHright.xvg -on ${dir}/MATH_MATHneighbour_contacts/numcont_MATH_${subunit}_MATHright.xvg -tu us -d 0.7 <<EOF
a_${MATH_start}-${MATH_end}
a_${MATHright_start}-${MATHright_end}
EOF

rm ${dir}/MATH_MATHneighbour_contacts/MATH_${subunit}_MATHright.ndx

#Update start of subunit MATH domain for next loop
MATH_start=$(( ${MATH_start} + ${atoms_per_subunit} ))

done





