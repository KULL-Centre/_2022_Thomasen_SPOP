#!/bin/bash
#

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

MATH_type=MATHfree

atoms_per_MATH=2245 #set up to (and with) 138VAL
atoms_per_subunit=5241
filename=prodrun_AAbackmapped

for i in 60
do

	cd ${MATH_type}_runs/SPOP_${i}mer

	subunit_nr1=$(( ($i/2)-2 ))
	subunit_nr2=$(( ($i/2)-1 ))
	subunit_nr3=$(( ($i/2) ))
	subunit_nr4=$(( ($i/2)+1 ))
	
	BTBBACK1_start=$(( ( $subunit_nr1 * $atoms_per_subunit ) + $atoms_per_MATH ))
	BTBBACK1_stop=$(( ($subunit_nr1+1) * $atoms_per_subunit ))
	BTBBACK2_start=$(( ( $subunit_nr2 * $atoms_per_subunit ) + $atoms_per_MATH ))
	BTBBACK2_stop=$(( ($subunit_nr2+1) * $atoms_per_subunit ))
	BTBBACK3_start=$(( ( $subunit_nr3 * $atoms_per_subunit ) + $atoms_per_MATH ))
	BTBBACK3_stop=$(( ($subunit_nr3+1) * $atoms_per_subunit ))
	BTBBACK4_start=$(( ( $subunit_nr4 * $atoms_per_subunit ) + $atoms_per_MATH ))
        BTBBACK4_stop=$(( ($subunit_nr4+1) * $atoms_per_subunit ))
	
	echo $BTBBACK1_start $BTBBACK1_stop $BTBBACK2_start $BTBBACK2_stop $BTBBACK3_start $BTBBACK3_stop $BTBBACK4_start $BTBBACK4_stop 

	$gmx make_ndx -f ${filename}.gro -o 4_center_domains.ndx <<EOF
a ${BTBBACK1_start}-${BTBBACK1_stop} | a ${BTBBACK2_start}-${BTBBACK2_stop} | a ${BTBBACK3_start}-${BTBBACK3_stop} | a ${BTBBACK4_start}-${BTBBACK4_stop}
q
EOF


	mkdir superposed_frames

	$gmx trjconv -s ${filename}.gro -f ${filename}.xtc -o superposed_frames/AA_frame.pdb -fit rot+trans -n 4_center_domains.ndx -skip 100 -sep <<EOF
a_${BTBBACK1_start}-${BTBBACK1_stop}_a_${BTBBACK2_start}-${BTBBACK2_stop}_a_${BTBBACK3_start}-${BTBBACK3_stop}_a_${BTBBACK4_start}-${BTBBACK4_stop}
1
EOF

	cd ../..
done
