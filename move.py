import sys
sys.path.append('/home/yao/scripts/toolkit')
import core

structure = core.Structure(sys.argv[1])
move_coord = structure.move(num_inp=False,dis_inp=False)
core.Writefile(coord=move_coord,lattice=structure.lattice()).tovasp(outfile_name='POSCAR_moved',direct=False)