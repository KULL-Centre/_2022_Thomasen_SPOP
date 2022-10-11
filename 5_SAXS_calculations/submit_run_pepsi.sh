#!/bin/bash

pepsiscript=run_pepsi_constantparams.py

for i in $(seq 2 2 60)
do

dir=SPOP_${i}mer/SAXS
mkdir $dir
cp $pepsiscript $dir
cp run_pepsi.sh $dir
cd $dir

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N ${i}mer_pepsiSAXS" >> temp
cat run_pepsi.sh >> temp
mv temp run_pepsi.sh

qsub run_pepsi.sh -v oligomersize=$i

cd ../..

done
