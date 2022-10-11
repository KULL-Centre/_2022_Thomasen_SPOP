#!/bin/bash

for i in $(seq 2 2 12)
do

cd SPOP_${i}mer/Backmapping
cp ../../Backmap.py .

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N ${i}mer_backmap" >> temp
cat ../../Backmap.sh >> temp
mv temp Backmap.sh

qsub Backmap.sh

cd ../..

done

