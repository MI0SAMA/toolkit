import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import matplotlib.cm as cm
import seaborn as sns

file_inp = sys.argv[1]
df_angle = pd.read_csv(file_inp,index_col=False)
df_heat = pd.DataFrame(np.zeros((18,18)))
for row in df_angle.itertuples():
    x = float(row[5])
    y = float(row[6])
    flag = True
    for a in np.arange(0,18):
        if x>(a*10-90) and x<(a*10-80):
            for b in np.arange(0,18):
                if y>(b*10-90) and y<(b*10-80):
                    df_heat.iloc[a][b]=df_heat.iloc[a][b]+1
                    flag = False
                    break
        if not flag:
            break
if df_heat.sum().sum()==df_angle.shape[0]:
    df_heat=df_heat
    plt.figure(figsize=(8, 6))
    g = sns.heatmap(df_heat, vmin=0, annot=True,  fmt='g', cmap='Blues', cbar=False, robust=True)
    plt.savefig(file_inp.split(".")[0]+'.png')
# plt.show()
#  print(getattr(row, 'Main_position').strip('[').strip(']').split( )[2])
# df1 = df.loc[:, ['Bond_Surface(deg)', 'Mid_Surface(deg)']]
# print (df1)
# df.loc[:, ['Bond_Surface(deg)']].hist(bins=100, color = "blue", grid =True, label = 'bond1', alpha=0.6)
# df.loc[:, ['Mid_Surface(deg)']].hist(bins=100, color = "red", grid =True, label = 'bond2', alpha=0.6)
# plt.show()