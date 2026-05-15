import sys
import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from ase.io import read
from ase.visualize.plot import plot_atoms
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

#get surface atoms by range
def getsurf_range(atoms, range_inp):
    pos = atoms.get_positions()
    return [
        i for i in range(len(atoms))\
        if min(range_inp) < pos[i][2] < max(range_inp)
    ]

#get surface atoms by mesh and search, assume adsorption of H atom
#compare the distance between mesh point with the atoms, scale counts a lot 
def getsurf_grid(atoms, xypos, scale=1.2):
    radiiAll = scale * np.array([
        covalent_radii[i] for i in atoms.get_atomic_numbers()])
    z_test = max(atoms.get_positions()[:, 2])+1.5
    surfIndex = []
    while len(surfIndex) == 0:
        if z_test<0:
            break
        testAtom = Atoms(['H'], [[xypos[0], xypos[1], z_test]])
        tmpAtoms = atoms.copy()
        tmpAtoms.extend(testAtom)
        dist_test = list(tmpAtoms.get_distances(-1, list(range(len(atoms)))) - radiiAll)
        for i in dist_test:
            if i < 0: surfIndex.append(dist_test.index(i))
        z_test -= 0.1
    return surfIndex
#get the mesh in XY plane, mesh counts a lot 
def get_surf_grid(atoms, mesh=20):
    cell = atoms.get_cell()
    grids = np.linspace(0, 1, mesh, endpoint=False)
    surfList = []
    for i in grids:
        for j in grids:
            surfList += getsurf_grid(atoms,[cell[0][0]*i+cell[1][0]*j, cell[0][2]*i+cell[1][1]*j])
    # for i,j in zip(grids, grids):
    #     surfList += getsurf_grid(atoms,\
    #         [cell[0][0]*i+cell[1][0]*j, cell[0][2]*i+cell[1][1]*j])
    return sorted(list(set(surfList)))

#prepare for calculating hollow and bridge site with period
def get_extended_atoms(atoms):
    tmpAtoms = atoms.copy()
    tmpAtoms = atoms*[3,3,1]
    tmpAtoms.set_positions(tmpAtoms.get_positions()\
        -atoms.get_cell()[0]-atoms.get_cell()[1])
    return tmpAtoms

#delet the atoms outside od the cell
def del_outside(list,cell):
    coord_d = np.dot(list,(np.linalg.inv(cell)))
    del_index=[]
    for i in range(len(list)):
        if coord_d[i].min()<0 or coord_d[i].max()>1:
            del_index.append(i)
    list=np.array([i for num,i in enumerate(list) if num not in del_index])
    return list

if __name__ == '__main__':
    #get input and surface atoms
    f_inp = sys.argv[1]
    if "." in f_inp:
        f_name = f_inp.split('.')[0]
    else:
        f_name = sys.argv[1]
    surf = read(f_inp)
    input1 = input('Input 1 or 2 for the searching method\n'
                  '1:search surface atoms by range\n'
                  '2:search by mesh and covalent radii\n')
    if input1 == '1':
        range_inp = input('Input the range like:9 12\n')
        range_inp = list(map(float, range_inp.split(' ')))
        surfList = getsurf_range(surf, range_inp)
    elif input1 == '2':
        surfList = get_surf_grid(surf)
    del surf[[i for i in range(len(surf)) if i not in surfList]]

    #consider the periodicity
    surf_ext = get_extended_atoms(surf)
    pos_ext = surf_ext.get_positions()
    top_atom = surf.get_positions()
    #hollow got through delaunay of the top sites
    tri = Delaunay(pos_ext[:, :2])
    pos_nodes = pos_ext[tri.simplices]
    hollow = np.array([t.sum(axis=0)/3 for t in pos_nodes])
    #bridge site got through half of the node
    bridge = []
    for i in pos_nodes:
        bridge.append((i[0]+i[1])/2)
        bridge.append((i[0]+i[2])/2)
        bridge.append((i[1]+i[2])/2)
    bridge = np.array(bridge)
    #delete the atoms outside
    cell_org = surf.get_cell()
    hollow = del_outside(hollow,cell_org)
    bridge = del_outside(bridge,cell_org)

    #draw graph 
    fig,ax=plt.subplots(2,1)
    ax[0].triplot(pos_ext[:,0], pos_ext[:,1], triangles=tri.simplices, color='grey',)
    ax[0].plot(hollow[:,0],  hollow[:,1],  'ok', label='Hollow')
    ax[0].plot(bridge[:,0],  bridge[:,1],  'or', label='Bridge')
    ax[0].plot(top_atom[:,0], top_atom[:,1], 'ob', label='Top')
    ax[0].set_xlim([cell_org[1][0],cell_org[0][0]])
    ax[0].set_ylim([cell_org[0][1],cell_org[1][1]])
    ax[0].legend(loc='lower left')
    plot_atoms(surf,ax=ax[1],radii=0.5,show_unit_cell=2)
    plt.savefig(f_name+'.png',  bbox_inches = "tight", facecolor=fig.get_facecolor(), transparent=True)