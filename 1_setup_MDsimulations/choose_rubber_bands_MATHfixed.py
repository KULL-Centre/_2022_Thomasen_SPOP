import sys
import numpy as np

subunits=int(sys.argv[1])
topfile="all_PRO.top"
grofile="PRO_CG.gro"
outfile="rubberbands_all_PRO.top"

#Residues for domain boundaries
residues = ["3LYS", "140GLN", "160ARG", "334ALA"]    

#Go through grofile and find atomnumbers of BB beads for the desired residues
#Append to list atompos
resi, beadtype = np.genfromtxt(grofile, usecols = (0,1), dtype=(str), skip_header=2, skip_footer=1, unpack=True)
atom = np.genfromtxt(grofile, usecols = 2, dtype=(int), skip_header=2, skip_footer=1, unpack=True)

atompos = []
for residue in residues:
    atompos_temp = []
    for i in range(len(resi)):
        if resi[i] == residue and beadtype[i] == "BB":
            atompos_temp.append(atom[i])
    atompos.append(atompos_temp)

print("Residue names: " + str(residues), "BB atoms in each subunit: " + str(atompos))

#Find Rubber band lines in topfile
with open(topfile, 'r') as f:
    toplines = f.readlines()

for i in range(len(toplines)):
    if '; Rubber band' in toplines[i]:
        start_line = int(i+1)
        break

for i in range(len(toplines)):
    if ';' in toplines[i] and i>start_line:
        end_line = int(i-1)
        break

print("Rubber bands are from line %s to line %s in %s" % (start_line, end_line, topfile))

#Returns a list of EN lines only between atoms of a protein domain from resi1 to resi2
def get_lines_to_keep_domain(resi1, resi2, subunits, toplines, start_line, end_line, residues):
    
    resi1index = residues.index(resi1)
    resi2index = residues.index(resi2)
    
    lines_to_keep = []
    
    for line in toplines[start_line:end_line]:
        linesplit = line.split()
        for i in range(subunits):
            if             int(linesplit[0]) >= int(atompos[resi1index][i]) and int(linesplit[0]) <= int(atompos[resi2index][i])             and int(linesplit[1]) >= int(atompos[resi1index][i]) and int(linesplit[1]) <= int(atompos[resi2index][i]):
                lines_to_keep.append(line)
    return lines_to_keep

#Returns a list of EN lines only between atoms of a protein domain from resi1 to resi2 and a protein domain from resi3 to resi4
def get_lines_to_keep_interdomain(resi1, resi2, resi3, resi4, subunits, toplines, start_line, end_line, residues):
    
    resi1index = residues.index(resi1)
    resi2index = residues.index(resi2)
    resi3index = residues.index(resi3)
    resi4index = residues.index(resi4)
    
    lines_to_keep = []
    
    for line in toplines[start_line:end_line]:
        linesplit = line.split()
        for i in range(subunits):
            if             int(linesplit[0]) >= int(atompos[resi1index][i]) and int(linesplit[0]) <= int(atompos[resi2index][i])             and int(linesplit[1]) >= int(atompos[resi3index][i]) and int(linesplit[1]) <= int(atompos[resi4index][i]):
                lines_to_keep.append(line)
                
            if             int(linesplit[0]) >= int(atompos[resi3index][i]) and int(linesplit[0]) <= int(atompos[resi4index][i])             and int(linesplit[1]) >= int(atompos[resi1index][i]) and int(linesplit[1]) <= int(atompos[resi2index][i]):
                lines_to_keep.append(line)
    
    return lines_to_keep

#Returns a list of EN lines only between an atom of a protein domain from resi1 to resi2 and the same domain one subunit over
def get_lines_to_keep_intersubunit(resi1, resi2, subunits, toplines, start_line, end_line, residues):
        
    resi1index = residues.index(resi1)
    resi2index = residues.index(resi2)
        
    lines_to_keep = []
        
    for line in toplines[start_line:end_line]:
        linesplit = line.split()
        for i in range(subunits-1):
            if             int(linesplit[0]) >= int(atompos[resi1index][i]) and int(linesplit[0]) <= int(atompos[resi2index][i])             and int(linesplit[1]) >= int(atompos[resi1index][i+1]) and int(linesplit[1]) <= int(atompos[resi2index][i+1]):                
                lines_to_keep.append(line)
            
            if             int(linesplit[0]) >= int(atompos[resi1index][i+1]) and int(linesplit[0]) <= int(atompos[resi2index][i+1])             and int(linesplit[1]) >= int(atompos[resi1index][i]) and int(linesplit[1]) <= int(atompos[resi2index][i]):
                lines_to_keep.append(line)
    
    return lines_to_keep

#Intradomain MATH domain
MATH = get_lines_to_keep_domain("3LYS", "140GLN", subunits, toplines, start_line, end_line, residues)
#Intradomain BTB/BACK domain
BTBBACK = get_lines_to_keep_domain("160ARG", "334ALA", subunits, toplines, start_line, end_line, residues)

#Interdomain MATH BTB/BACK domains
MATH_BTBBACK = get_lines_to_keep_interdomain("3LYS", "140GLN", "160ARG", "334ALA", subunits, toplines, start_line, end_line, residues)

#Intersubunit BTB/BACK domains
BTBBACK_intersubunit = get_lines_to_keep_intersubunit("160ARG", "334ALA", subunits, toplines, start_line, end_line, residues)


#Write to file
with open(outfile, 'a') as f:
    for line in toplines[:start_line]:
        f.write(line)
    
    for line in MATH:
        f.write(line)
    for line in BTBBACK:
        f.write(line)
    for line in MATH_BTBBACK:
        f.write(line)
    for line in BTBBACK_intersubunit:
        f.write(line)
        
    for line in toplines[end_line:]:
        f.write(line)



