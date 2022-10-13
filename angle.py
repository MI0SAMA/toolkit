#!/usr/bin/env python
import time
import pandas as pd
import sys
sys.path.append('/home/yao/scripts/toolkit')
import core
import os

st = time.process_time()
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
print (lattice_input)
# lattice_input = (13.55300000,0.00000000,0.00000000,6.77999952,11.73500057,0.00000000,4.21399883,2.43199956,17.70789335)
normal, abnormal = pd.DataFrame(), pd.DataFrame()
start,end = int(sys.argv[1]),int(sys.argv[2])
for i in range (start,end+1):
    num=str(i).zfill(4)
    file_inp='pos_'+num+'.xyz'
    a, b =core.Analysis(file_inp,lattice=lattice_input).angle(period='xy')
    normal = pd.concat((normal,a),axis=0)
    abnormal = pd.concat((abnormal,b),axis=0)
    print (file_inp)
normal.index.name,abnormal.index.name='atmo_index','atmo_index'
file1, file2= 'angle_normal_'+str(start)+'_'+str(end)+'.csv','angle_abnormal_'+str(start)+'_'+str(end)+'.csv'
normal.to_csv(file1)
abnormal.to_csv(file2)
# print ('LOOP'+str(j+1)+'DONE')
et = time.process_time()
print (et-st)