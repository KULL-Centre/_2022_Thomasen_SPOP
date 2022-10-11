#!/bin/bash

input1=12
input2=12
min_size=14
max_size=$(( ${input1} + ${input2} - 4 ))

for i in $(seq ${min_size} 2 ${max_size})
do

echo "Submitting ${i}-mer"

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N ${i}mer_SPOP_stitch" >> temp
cat Stitch_oligomers.sh >> temp
mv temp Stitch_oligomers_tmp.sh

qsub Stitch_oligomers_tmp.sh -v input1=$input1,input2=$input2,target=$i

rm Stitch_oligomers_tmp.sh

done
