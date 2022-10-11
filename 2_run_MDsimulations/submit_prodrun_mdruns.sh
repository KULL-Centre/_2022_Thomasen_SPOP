#!/bin/bash

for i in $(seq 2 2 12)
do

cd SPOP_${i}mer
cp ../prodrun_mdrun.sh .

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N ${i}mer_SPOP_md" >> temp
cat prodrun_mdrun.sh >> temp
mv temp prodrun_mdrun.sh

qsub prodrun_mdrun.sh
cd ..

done
