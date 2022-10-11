#!/bin/bash
#
#Script to coarse-grain and run energy minimization of SPOP oligomers using martini


export PATH="/lindorffgrp-isilon/wyong/software/miniconda3/bin:$PATH"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lindorffgrp-isilon/wyong/software/openmpi401/lib

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi
python=/lindorffgrp-isilon/wyong/software/miniconda3/bin/python3

martinize=/lindorffgrp-isilon/wyong/software/miniconda3/bin/martinize2

#wget http://cgmartini.nl/images/tools/insane/insane.py
insane=insane.py

minmdp=minimization.mdp
FF=martini304
ffdir=/lindorffgrp-isilon/wyong/LPMO_Gaston/Martini30b417
posre_script=Write_posre_xline.py
EN_script=choose_rubber_bands_MATHfree.py

dssp=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/mkdssp

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
#export GMX_MAXCONSTRWARN=-1

#Loop through oligomer-sizes
for i in $(seq 6 2 12)
do

#Starting structure
pdb=SPOP_${i}mer.pdb


#Make directory, cp .pdb and insane.py file there and cd there
dir=SPOP_${i}mer_xlineposres0.005
#rm -r $dir
mkdir $dir
cp Structures/$pdb $dir
cp $insane $dir
cp $posre_script $dir
cp $EN_script $dir
cd $dir


#Remove C-terminal OXT-atom from .pdb file
sed -i '/OXT/d' $pdb


#Martinize
$python $martinize -f $pdb -o PRO_topol.top -x PRO_CG.pdb -dssp $dssp -ff $FF -ed -cys auto -elastic -ef 500 -el 0.5 -eu 1.2 -nt -scfix -ff-dir $ffdir/v.3.0.4.17/martini3-protein/force_fields/ -map-dir $ffdir/v.3.0.4.17/martini3-protein/mappings/ 

boxlength=$(( 15+$i*3 ))

#Put protein in box
if [[ "$i" -le "8" ]]
then
	$gmx editconf -f PRO_CG.pdb -o PRO_CG.gro -bt triclinic -box $boxlength 21 21 -center 0 0 0 -princ <<EOF
1
EOF
fi

if [[ "$i" -ge "10" ]]
then
	$gmx editconf -f PRO_CG.pdb -o PRO_CG.gro -bt triclinic -box $boxlength 23 23 -center 0 0 0 -princ <<EOF
1
EOF
fi

#Solvate using insane.py
python2.7 $insane -f PRO_CG.gro -o PRO_SOL_IONS.gro -pbc keep -salt 0.15 -sol W -center -p PRO_topol_SOL_IONS.top


#The next few blocks modify the toplogy file and molecule_0.itp file:

#Remove #include include martini.itp and substitute ion names in topology file
perl -pi -e's/#include "martini.itp"//g' PRO_topol_SOL_IONS.top
perl -pi -e's/NA\+/NA/g' PRO_topol_SOL_IONS.top
perl -pi -e's/CL-/CL/g' PRO_topol_SOL_IONS.top

#Rename molecule_0.itp to PRO.itp and rename "molecule_0" as "Protein" in PRO.itp file
mv molecule_0.itp PRO.itp
perl -pi -e's/molecule_0/Protein/g' PRO.itp

#Add "#include .itp" lines to PRO_topol_SOL_IONS.top
cat <<EOF > others.top
#include "$ffdir/v.3.0.4.17/martini_v3.0.4.itp"
#include "PRO.itp"
#include "$ffdir/v.3.0.4.17/martini_v3.0_phospholipids.itp"
#include "$ffdir/v.3.0.4.17/martini_v3.0_ions.itp"
#include "$ffdir/v.3.0.4.17/martini_v3.0_solvents.itp"
EOF
cat others.top PRO_topol_SOL_IONS.top >a
mv a PRO_topol_SOL_IONS.top

#Write posres
$python $posre_script $i

#Run energy minimization
$gmx grompp -f ../$minmdp -p PRO_topol_SOL_IONS.top -c PRO_SOL_IONS.gro -o min.tpr -pp all_PRO.top -maxwarn 3 -r PRO_SOL_IONS.gro
nohup $gmx mdrun -deffnm min -v -ntomp 8 &

nohup $python $EN_script $i &

rm $insane

cd ..

done
