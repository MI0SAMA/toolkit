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
#draw contor with image and line
contor1=plt.contour(xi,yi,zi,colors='black',linewidths=1)
contor2=plt.imshow(zi, cmap=plt.cm.jet, extent=(0,1,0,1), origin='lower')
plt.xlabel('Glide Position (a.u.)',fontproperties = 'Times New Roman',fontsize=24)
plt.ylabel('Glide Position (a.u.)',fontproperties = 'Times New Roman',fontsize=24)
plt.yticks(fontproperties = 'Times New Roman', size = 18)
plt.xticks(fontproperties = 'Times New Roman', size = 18)
plt.clabel(contor1, inline=1, fontsize=12,colors='black')
cb = plt.colorbar(contor2)
cb.set_label('Energy difference (eV)', fontproperties = 'Times New Roman',fontsize=24)
cb.ax.tick_params(labelsize=18)
plt.savefig("contour.jpg", dpi=300, pad_inches=0,bbox_inches='tight')
#scatter
plt.figure()
x,y=[],[]
for i in vec:
	x.append(i[0])
	y.append(i[1])
color = [plt.get_cmap("seismic",36)(int(float(i-ene.min())/(ene.max()-ene.min())*36)) for i in ene]
plt.set_cmap(plt.get_cmap("seismic", 36))
im=plt.scatter(x, y, s=800,c=color,marker='.')
cb = plt.colorbar(im, format=matplotlib.ticker.FuncFormatter(lambda x,pos:int(x*(ene.max()-ene.min()))+ene.min()))
# set grid
# plt.grid(True, color = 'r')
cb.set_label('Energy (eV)', fontproperties = 'Times New Roman',fontsize=24)
plt.xlabel('Glide Position (a.u.)',fontproperties = 'Times New Roman',fontsize=24)
plt.ylabel('Glide Position (eV)',fontproperties = 'Times New Roman',fontsize=24)
plt.xticks(np.arange(0, 1.1, 1/grid),fontproperties = 'Times New Roman',fontsize=18)
plt.yticks(np.arange(0, 1.1, 1/grid),fontproperties = 'Times New Roman',fontsize=18)
cb.ax.tick_params(labelsize=18)
plt.grid(True, color = 'black')
plt.savefig("scatter.jpg", dpi=300, bbox_inches='tight')
#3d energy surface
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(xi,yi,zi,cmap='rainbow') 
plt.xlabel('Relative position (a.u.)',fontsize=14)
plt.ylabel('Relative position (a.u.)',fontsize=14)
plt.savefig("energy_surface.jpg", dpi=300, pad_inches=0)
#plt.show()


