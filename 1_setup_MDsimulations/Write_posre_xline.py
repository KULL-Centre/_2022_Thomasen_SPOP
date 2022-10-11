import sys
import numpy as np

#Write gromacs position restraints to keep BB-atoms in BTB/BACK of first and last subunit on a line along the x-axis 

oligomersize = int(sys.argv[1])

PROitp = "PRO.itp"

beadtype = np.genfromtxt('PRO_CG.gro', skip_header=2, skip_footer=1, usecols=(1), dtype = str, unpack = True)
atomnum = np.genfromtxt('PRO_CG.gro', skip_header=2, skip_footer=1, usecols=(2), dtype = int, unpack = True)

atomnums = []

atoms_per_subunit = 756

BTBBACK_start = 398
BTBBACK_end = 643

for i in range(len(atomnum)):
    if atomnum[i] in range(BTBBACK_start, BTBBACK_end) or atomnum[i] in range(atoms_per_subunit*(oligomersize-1)+BTBBACK_start, atoms_per_subunit*(oligomersize-1)+BTBBACK_end):
        if beadtype[i] == 'BB':
            atomnums.append(atomnum[i])

with open(PROitp, 'a') as output:
    output.write("[ position_restraints ] \n ;  resid  func.   kx  ky  kz\n")
    for j in atomnums:
        output.write("%i \t1 \t0 \t0.005 \t0.005 \t  ; Restrain to x-axis line\n" % j)

print("Wrote flatbottom posres to PRO.itp file")
