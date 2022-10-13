#!/usr/bin/env python
import numpy as np 
import pandas as pd
import sys
import core
import os
import time
import matplotlib.pyplot as plt

#input_file = sys.argv[1]
#lattice for pt
#lattice_input = (13.55300000,0.00000000,0.00000000,6.77999952,11.73500057,0.00000000,4.21399883,2.43199956,17.70789335)
#lattice for 1t
#lattice_input = (13,13,14.625975,90,90,120)
#lattice for 1t_vac
#lattice_input = (13,0,0,-6.5,11.25833025,0,0,0,30)
#lattice for 1t_112
#lattice_input = (13,13,26.138826,90,90,120)
#lattice for 1t_64
# lattice_input = (13,0,0,-6.5,11.25833025,0,0,0,16.139999)

st = time.process_time()
#read lattice from cp2k.inp file
lattice_input = []
if os.path.exists('cp2k.inp'):
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
print ('Read lattice information from .inp file:',lattice_input)
#default setting for density anlysis
minz = 0
maxz = 27
intervalz = float(sys.argv[3])
water_num = 64

a=1
start_num = int(sys.argv[1])
end_num = int(sys.argv[2])
all_c,all_d = pd.Series(core.Analysis('pos_'+str(start_num).zfill(4)+'.xyz',lattice=lattice_input).density(interval=intervalz,min=minz,max=maxz,axis='z',element='O'))
for i in np.arange(start_num+1,end_num+1):
    input_file='pos_'+str(i).zfill(4)+'.xyz'
    density_c,density_d = core.Analysis(input_file,lattice=lattice_input).density(interval=intervalz,min=minz,max=maxz,axis='z',element='O')
    a+=1
    if density_c.sum() == water_num:
        print (input_file)
        all_d = all_d + density_d
        all_c = all_c + density_c
    else:
        print (input_file,'There might be some atoms missed')
all_d = all_d/a*18/6.02*10
print (all_d)
all_d.to_csv('density_'+str(start_num).zfill(4)+'_'+str(end_num).zfill(4)+'_'+str(intervalz)+'.csv',index_label=False,header=None)
all_d.plot(kind='line',marker='o',color='r')
plt.savefig('density_'+str(start_num).zfill(4)+'_'+str(end_num).zfill(4)+'_'+str(intervalz)+'.png')
plt.show()
et = time.process_time()
print ('Consumed time: ',et-st)