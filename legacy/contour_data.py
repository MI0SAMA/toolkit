import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import os
#get energy log
if os.path.exists('energy.log'):
    pass
else:
    print(os.system("for i in {1..49}; do cd opt_$i; \
a=`grep without OUTCAR | tail -1 | awk '{print $5}'`; cd ..; echo $a >> ./energy.log; done"))

ene=[]
with open('energy.log','r') as f:
    lines = [x.rstrip().split() for x in f]
    for i in lines:
         ene.append(float(i[0]))
ene = np.array(ene)
ene = ene - ene.min()

#rgenerate posotion
grid=np.sqrt(len(ene))-1
vec = []
a=np.array([1/grid,0])
b=np.array([0,1/grid])
for i in range(int(grid+1)):
    a1 = i*a
    for j in range(int(grid+1)):
        b1 = j*b
        vec.append((a1+b1).tolist())
#read enegy of each position
#generate data of z
xi=np.linspace(0,1,500)
yi=np.linspace(0,1,500)
xi,yi=np.meshgrid(xi,yi)
zi=griddata(vec,ene.tolist(),(xi,yi),method='cubic')#linear nearest
#position = np.vstack([xi.ravel(), yi.ravel()])
with open('energy_surface.dat','w') as f:
    f.truncate()
    for i in np.arange(0,500,1):
        for j in np.arange(0,500,1):
            z = zi[i][j]
            # x = xi[i][j]+0.1-1 if (xi[i][j]+0.1) > 1 else xi[i][j]+0.1
            # y = yi[i][j]+0.1-1 if (yi[i][j]+0.1) > 1 else yi[i][j]+0.1
            x,y = xi[i][j],yi[i][j]
            f.write(str(x)+' '+str(y)+' '+str(z)+' '+'\n')
