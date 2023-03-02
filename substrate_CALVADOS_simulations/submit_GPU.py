import subprocess
import os
import pandas as pd
import numpy as np
import mdtraj as md
import time
from jinja2 import Template

proteins = ['SCAF1_IDR1', 'SETD2_IDR1']

submission = Template("""#!/bin/sh
#SBATCH --job-name={{name}}
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH -t 1000:00:00
#SBATCH --partition=sbinlab_gpu
#SBATCH -o {{path}}/out
#SBATCH -e {{path}}/err

source /groups/sbinlab/giulio/.bashrc

conda activate hoomd

python ./simulate_GPU.py --seq_name {{name}} --path {{path}}""")

for name in proteins:
    os.system(f'rm -r {name}')
    if not os.path.isdir(name):
        os.mkdir(name)
    with open('{:s}.sh'.format(name), 'w') as submit:
        submit.write(submission.render(name=name,path=name))
    subprocess.run(['sbatch','{:s}.sh'.format(name)])
    print(name)
    time.sleep(.6)
