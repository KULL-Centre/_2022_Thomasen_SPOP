#PBS -l nodes=1:ppn=40:thinnode
#PBS -l walltime=1000:00:00
#PBS -l mem=130gb
# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

python=/home/projects/ku_10001/people/fretho/miniconda3/bin/python3.8
pepsiscript=run_pepsi_constantparams.py

$python $pepsiscript $oligomersize
