import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import matplotlib.cm as cm
import seaborn as sns

file_inp = sys.argv[1]
df_angle = pd.read_csv(file_inp,index_col=False)
name_list = ['Bond(deg)','Bond_Surface(deg)','Mid_Surface(deg)']
# print (df_angle[name_list[0]])
def cal_poss(ang_list):
    # ang_list = ang_list1.tolist()
    pos_list = []
    for i in range(-90,90,10):
        num = 0
        for j in range(len(ang_list)):
            if float(ang_list[j]) > float(i) and float(ang_list[j]) < float(i+10):
                # ang_list.remove(ang_list[j])
                num += 1
        num = float(num)/len(ang_list)
        pos_list.append(num)
    return pos_list
y1,y2=cal_poss(df_angle[name_list[1]]),cal_poss(df_angle[name_list[2]])
x1 = range(-90,90,10)
l1=plt.plot(x1,y1,'r--',label=name_list[1])
l2=plt.plot(x1,y2,'g--',label=name_list[2])
plt.plot(x1,y1,'ro-',x1,y2,'g+-')
plt.title('Possbility distribution')
plt.xlabel('Angle (deg)')
plt.ylabel('Possbility (a.u.)')
plt.legend()
plt.savefig('pos_'+sys.argv[1].split('.')[0])
plt.show()