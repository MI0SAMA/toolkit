import numpy as np
import pandas as pd
import os

#read type
type_dic = {}
with open ('../type_map.raw') as f:
    lines = [x.rstrip().split() for x in f]
    for i in lines:
        type_dic[i[0]]=0
with open('../type.raw') as f:
    lines = [x.rstrip().split() for x in f]
    for i in lines:
        ele = list(type_dic.keys())[int(i[0])]
        type_dic[ele] += 1
force = np.load('force.npy')
ele_list=list(type_dic.keys())
force_dic=type_dic
print (type_dic)
for i in force:
    for j in range(len(ele_list)):
        if j == 0:
            sta=0
            end=type_dic[ele_list[j]]*3-1
        else:
            sta=type_dic[ele_list[j-1]]*3
            end=(type_dic[ele_list[j]]+type_dic[ele_list[j-1]])*3-1
        print (list(i[sta:end]))
    break