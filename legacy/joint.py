import sys
import pandas as pd
import numpy as np
sys.path.append('/home/yao/scripts/toolkit')
import core

try:
    structure1 = core.Structure(sys.argv[1])
    structure2 = core.Structure(sys.argv[2])
except:
    print ('Please put two structure file like: python joint.py POSCAR1 POSCAR2')
joint_axis = input("Input the joint axisas x, y, z: ")
lattice1 = structure1.lattice()
lattice2 = structure2.lattice()
move_dis = float(lattice1[2][2])
layer_dis = 0.5
print (lattice1)
print (lattice2)
coord1, coord2 = structure1.coord(direct=False), structure2.coord(direct=False)
#coord1[joint_axis] = coord1[joint_axis].apply(lambda x: x/2)
#coord2[joint_axis] = coord2[joint_axis].apply(lambda x: x/2+1) #distance between two structures
coord1[joint_axis] = coord1[joint_axis]
coord2[joint_axis] = coord2[joint_axis].apply(lambda x: x+move_dis+layer_dis) #distance between two structures
coord_new = pd.concat([coord1,coord2])
dic= {'x':0, 'y':1, 'z':2}
list = [0,1,2]
list.pop(dic[joint_axis])
if (lattice1[list[0]] == lattice2[list[0]]).all() and \
   (lattice1[list[1]] == lattice2[list[1]]).all():
    pass
else:
    print ('The lattices do not match in the jointing direction')
    ignore = input ('Ignore it?: (y or n)')
    if ignore == 'y':
        pass
    else:
        sys.exit(0)
lattice_new = lattice1
lattice_new[dic[joint_axis]] = lattice1[dic[joint_axis]]+lattice2[dic[joint_axis]]
print (lattice_new)
lattice_new[2][2]=lattice_new[2][2]+layer_dis
print (lattice_new)
#coord_newc = np.dot(coord_new,lattice_new)
coord_newc = pd.DataFrame(coord_new, columns=['x','y','z'], index=coord_new.index)
structure = core.Writefile(coord=coord_newc,lattice=lattice_new)
structure.tovasp(outfile_name='POSCAR_joint',direct=True)
