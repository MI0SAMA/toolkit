import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import matplotlib.cm as cm
import seaborn as sns

list1=[-90,-20,20,90]
list2=[-90,-50,50,90]
file_inp = sys.argv[1]
df_angle = pd.read_csv(file_inp,index_col=False)
df_heat = pd.DataFrame(np.zeros((3,3)))
for row in df_angle.itertuples():
    if float(sys.argv[2])<float(row[7].strip('[').strip(']').split()[2])<float(sys.argv[3]):
        x = float(row[5])
        y = float(row[6])
        flag = True
        for a in range(0,3):
            if x>list1[a] and x<list1[a+1]:
                for b in range(0,3):
                    if y>list2[b] and y<list2[b+1]:
                        df_heat.iloc[a][b]=df_heat.iloc[a][b]+1
                        flag = False
                        break
            if not flag:
                break
# if df_heat.sum().sum()==df_angle.shape[0]:
df_heat=df_heat
plt.figure(figsize=(8, 6))
g = sns.heatmap(df_heat, annot=True, fmt='g', cmap='Blues', cbar=False, robust=True
               ,xticklabels=False, yticklabels=False, annot_kws={'fontsize': 22})
plt.savefig(file_inp.split(".")[0]+'_near_3.png')
plt.show()
s1,s2,s3,s4=df_heat.iloc[0][1]+df_heat.iloc[2][1],df_heat.iloc[2][2],df_heat.iloc[0][0],df_heat.iloc[1][0]+df_heat.iloc[1][1]+df_heat.iloc[1][2]
total=df_heat.sum().sum()
if s1+s2+s3+s4==total:
    s1,s2,s3,s4=s1/total,s2/total,s3/total,s4/total
    print (s1,s2,s3,s4)
else:
    print (df_heat,s1/total,s2/total,s3/total,s4/total)
#  print(getattr(row, 'Main_position').strip('[').strip(']').split( )[2])
# df1 = df.loc[:, ['Bond_Surface(deg)', 'Mid_Surface(deg)']]
# print (df1)
# df.loc[:, ['Bond_Surface(deg)']].hist(bins=100, color = "blue", grid =True, label = 'bond1', alpha=0.6)
# df.loc[:, ['Mid_Surface(deg)']].hist(bins=100, color = "red", grid =True, label = 'bond2', alpha=0.6)
# plt.show()