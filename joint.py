import sys
import pandas as pd
import numpy as np
sys.path.append('/home/yao/scripts/toolkit')
import core

structure1 = core.Structure(sys.argv[1])
structure2 = core.Structure(sys.argv[2])
joint_axis = input("Input the joint axisas x, y, z: ")
coord1, coord2 = structure1.coord(direct=True), structure2.coord(direct=True)
coord1[joint_axis] = coord1[joint_axis].apply(lambda x: x/2)
coord2[joint_axis] = coord2[joint_axis].apply(lambda x: x/2+0.5)
coord_new = pd.concat([coord1,coord2])
dic= {'x':0, 'y':1, 'z':2}
lattice1 = structure1.lattice()
lattice2 = structure1.lattice()
lattice_new = lattice1
lattice_new[dic[joint_axis]] = lattice1[dic[joint_axis]]+lattice2[dic[joint_axis]]
coord_newc = np.dot(coord_new,lattice_new)
coord_newc = pd.DataFrame(coord_newc, columns=['x','y','z'], index=coord_new.index)
structure = core.Writefile(coord=coord_newc,lattice=lattice_new)
structure.tovasp(outfile_name='POSCAR_joint',direct=True)