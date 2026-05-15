import sys
import pandas as pd
import numpy as np
import math
import core

try:
    structure = core.Structure(sys.argv[1])
except:
    print ('Please put two structure file like: python weight.py POSCAR')

ext_lattice,coord = structure.extend(1,1,1)
print (coord)
#coord['z'] = coord['z'].map(lambda x: x-0.5)
ele_list = coord.index.tolist()
ele_type = []
[ele_type.append(i) for i in ele_list if not i in ele_type ]
weight = coord.mean().tolist()
def distance(array):
    v = array-np.array(weight)
    return (math.hypot(v[0],v[1],v[2]))
dic={}
for i in ele_type:
    dis_list = []
    num = coord.loc[i].shape[0]
    for j in range(num):
        dis_list.append(distance(np.array(coord.loc[i].iloc[j])))
    std = np.std(dis_list,ddof=1)
    dic[i]=std
print ('standard deviation \n',dic)
