import subprocess
import os
import pandas as pd
import numpy as np
import mdtraj as md
import time
from jinja2 import Template

proteins = ['SETD2_IDR1', 'SETD2_IDR2', 'SCAF1_IDR1', 'SRC3_IDR1', 'Gli3_IDR1', 'Gli3_IDR2', 'Gli2_IDR2']

submission = Template("""#!/bin/sh
#SBATCH --job-name={{name}}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --nodelist=node592
#SBATCH --partition=sbinlab_ib2
#SBATCH --mem=40GB
#SBATCH -t 1000:00:00
#SBATCH -o {{path}}/out
#SBATCH -e {{path}}/err

source /groups/sbinlab/giulio/.bashrc

conda activate hoomd

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./simulate.py --seq_name {{name}} --path {{path}}""")

for name in proteins:
    os.system(f'rm -r {name}')
    if not os.path.isdir(name):
        os.mkdir(name)
    with open('{:s}.sh'.format(name), 'w') as submit:
        submit.write(submission.render(name=name,path=name))
    subprocess.run(['sbatch','{:s}.sh'.format(name)])
    print(name)
    time.sleep(.6)
