import sys
import os

from matplotlib import lines
sys.path.append('/home/yao/scripts/toolkit')
import core

#lattice =(13.00000000,0.00000000,0.00000000,-6.50000000,11.25833025,0.00000000,0.00000000,0.00000000,11.43675000)
if os.path.exists('cp2k.inp'):
    lattice_input = []
    with open('cp2k.inp','r') as f:
        lines = [x.rstrip().split() for x in f]
    for i in lines:
        if i != [] and i[0] == 'ABC':
            for y in [float(x) for x in i[1:]]:
                lattice_input.append(y)
        elif i != [] and i[0] == 'ALPHA_BETA_GAMMA':
            for y in [float(x) for x in i[1:]]:
                lattice_input.append(y)
            break
        elif i != [] and i[0] == 'A':
            for y in [float(x) for x in i[1:]]:
                lattice_input.append(y)
        elif i != [] and i[0] == 'B':
            for y in [float(x) for x in i[1:]]:
                lattice_input.append(y)
        elif i != [] and i[0] == 'C':
            for y in [float(x) for x in i[1:]]:
                lattice_input.append(y)
            break
else:
    lattice_input = False

structure = core.Writefile(sys.argv[1],lattice=lattice_input)
structure.tovasp(outfile_name=sys.argv[2],direct=False)
