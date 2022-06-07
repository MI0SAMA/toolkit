#!/usr/bin/env python
import numpy as np 
import pandas as pd
from scipy import spatial
import itertools
import sys

"""
USAGE
object can be a xyz,vasp,cif structure file
"""
class Structure(object):
    #initial
    def __init__(self,filename,lattice=False):
        self._filename = str(filename)
        if "." in self._filename:
            self._filetype = self._filename.split('.')[-1].lower()
        elif 'poscar' in self._filename.lower():
            self._filetype = 'vasp'
        with open(self._filename,'r') as f:
            self._filecontent = [x.rstrip().split() for x in f]
        self.lattice_input = lattice

    @property
    def filecontent(self) -> list:
        return self._filecontent
    @property
    def filename(self) -> str:
        return self._filename
    @property
    def filetype(self) -> str:
        if self._filetype in ['vasp','cif','xyz']:
            return self._filetype
        else:
            raise ValueError('Sturcture input wrong')
    
    #universal method
    def f2l(self,list,f):
        return [f(x) for x in list]
    
    def is_num(self,s):
        try:
            float(s)
            return True
        except ValueError:
            pass
        try:
            import unicodedata
            unicodedata.numeric(s)
            return True
        except (TypeError, ValueError):
            pass
        return False
    
    def cell_p2v(self,A,B,C,alpha,beta,gamma):
        alpha,beta,gamma = self.f2l([alpha,beta,gamma],lambda x: np.deg2rad(x))
        lx = A
        xy = B * np.cos(gamma)
        xz = C * np.cos(beta)
        ly = B* np.sin(gamma)
        if not ly > 0.1:
            raise RuntimeError("ly:=B* np.sin(gamma)=={}, must be greater than 0.1",format(ly))
        yz = (B*C*np.cos(alpha)-xy*xz)/ly
        if not C**2-xz**2-yz**2 > 0.01:
            raise RuntimeError("lz^2:=C**2-xz**2-yz**2=={}, must be greater than 0.01",format(C**2-xz**2-yz**2))
        lz = np.sqrt(C**2-xz**2-yz**2)
        lattice = np.array([[lx, 0 , 0],
        [xy, ly, 0 ],
        [xz, yz, lz]],dtype=float)
        return lattice

    #method
    """
    USAGE
    get lattice information as a 3*3 nparray
    for xyz file, the information should be input by fuction parameter as:
    A,B,C,alpha,beta,gamma or xx,xy,xz,yx,yy,yz,zx,zy,zz
    """
    def lattice(self):#get lattice as a nparray
        if self.filetype == 'vasp':
            lattice = np.array(self.filecontent[2:5],dtype=float)
        elif self.filetype == 'cif':
            dic = {}
            for i in self.filecontent:#Improve
                if len(i)==2:
                    dic[i[0]] = i[1]
            lg = np.array([dic['_cell_length_a'],dic['_cell_length_b'],dic['_cell_length_c']],dtype=float)
            ag = np.array([dic['_cell_angle_alpha'],dic['_cell_angle_beta'],dic['_cell_angle_gamma']],dtype=float)
            lattice = self.cell_p2v(lg[0],lg[1],lg[2],ag[0],ag[1],ag[2])#Improve 
        elif self.filetype == 'xyz' and self.lattice_input:
            if len(self.lattice_input) == 6:
                A,B,C,alpha,beta,gamma=self.lattice_input
                lattice = self.cell_p2v(A,B,C,alpha,beta,gamma) 
            elif len(self.lattice_input) == 9:
                lattice = np.zeros((3,3))
                for i in range(3):
                    lattice[i][:] = self.lattice_input[3*i:3*i+3]
        else:
            raise ValueError('Input of lattice information wrong')
        lattice[np.logical_and(lattice<0.0001,lattice>0)]=0
        return lattice

    """
    USAGE
    get absolute coordinate as pdframe
    """
    def coord(self):
        if self.filetype == 'vasp':
            at_num = list(map(int,self.filecontent[6]))
            at_type = dict(zip(self.filecontent[5],at_num))
            for i in 8,9:
                if self.is_num(self.filecontent[i][0]):
                    coord = np.array(self.filecontent[i:i+np.sum(at_num)])[:,:3].astype(float)
                    if self.filecontent[i-1][0][0] == 'D':
                        coord = np.dot(coord,self.lattice()).astype(float)
                        break
                    else:
                        break
                else:
                    pass
            coord = pd.DataFrame(coord,columns=['x','y','z'],index=[val for val in at_type.keys() for i in range(at_type[val])])    
        elif self.filetype == 'xyz':
            at_num = int(self.filecontent[0][0])
            coord = np.array(self.filecontent[2:2+at_num])[:,0:4]
            coord = pd.DataFrame(coord[:,1:4].astype(float),columns=['x','y','z'],index=[coord[:,0]])
        elif self.filetype == 'cif':
            count = 0
            for i in self.filecontent:
                if i !=[] and i[0] == '_atom_site_type_symbol':
                    break
                else:
                    count += 1
            coord = np.array(self.filecontent[count+1:])
            coord = pd.DataFrame(np.dot(coord[:,2:5].astype(float),self.lattice()),columns=['x','y','z'],index=[coord[:,-1]])
        return coord

    """
    USAGE
    x,y,z is the extend size,*args could be lattice information
    returen lattice and abosolute coord information after extend
    """
    def extend(self,x,y,z):
        org_coord = self.coord()
        org_lattice = self.lattice()
        ext_lattice = org_lattice * [x,y,z]
        def list_gen(x,y,z):
            alllist = []
            for max in (x,y,z):
                i = 0
                list = []
                while (i < max):
                    list.append(i+1)
                    i+=1
                alllist.append(list)
            return alllist
        extend_list = np.array(list(itertools.product(list_gen(x,y,z)[0],list_gen(x,y,z)[1],list_gen(x,y,z)[2])))
        ext_coord = pd.DataFrame()
        for i in extend_list:
            con_coord = org_coord + np.dot((i-[1,1,1]),org_lattice)#important
            ext_coord = pd.concat([ext_coord,con_coord])
        return ext_lattice,ext_coord

"""
use coord lattice or structure to generate POSCAR or xyz file
"""
class Writefile(object):
    def __init__(self, filename=False, coord=False, lattice=False):
        if filename and not coord and not lattice:
            self._filename = str(filename)
            self._coord = Structure(filename).coord()
            self._lattice = Structure(filename).lattice()
        elif filename and not coord:
            self._filename = str(filename)
            self._coord = Structure(filename).coord()
            if isinstance(lattice,np.ndarray):
                self._lattice = lattice
            else:
                self._lattice = Structure(filename,lattice=lattice).lattice()
        elif isinstance(coord,pd.DataFrame):
            self._coord = coord
            if isinstance(lattice,np.ndarray):
                self._lattice = lattice
            else:
                self._lattice = Structure(filename,lattice=lattice).lattice()    
        else:
            raise ValueError('Input of sturcture information wrong')
    
    @property
    def filename(self) -> str:
        if self._filename:
            return self._filename
        else:
            return False
    @property
    def coord(self) -> pd.DataFrame:
        return self._coord.sort_index() #the origin coord will be sorted by element type
    @property
    def lattice(self) -> np.array:
        return self._lattice
    
    """
    USAGE
    from structure file or cooed and lattice information to generate POSCAR
    the fix module will supported latter
    """
    def tovasp(self,outfile_name=False,direct=True):
        ele_type=self.coord.index.unique()
        ele_dic={}
        for i in ele_type:
            ele_dic[i[0]] = len(self.coord.loc[[i[0]]])
        if outfile_name == False:
            file_out = 'POSCAR'
        elif outfile_name:
            file_out = str(outfile_name)
        f_out = open(file_out, 'w')
        f_out.truncate()
        f_out.write(
        'Generated by toolkit@yao'+'\n'
        +'1.00000000000000'+'\n'
        +'\n'.join('\t'.join('%0.6f' %x for x in y) for y in self.lattice)+'\n'
        +'\t'.join(ele_dic.keys())+'\n'
        +'\t'.join(str(i) for i in ele_dic.values())+'\n'
        )
        if direct:
            direct_coord = np.dot(self.coord,(np.linalg.inv(self.lattice)))
            f_out.write(
            'D'+'\n'+
            '\n'.join('\t'.join('%0.6f' %x for x in y) for y in direct_coord)
            )
        elif not direct:
            f_out.write(
            'C'+'\n'+
            '\n'.join('\t'.join('%0.6f' %x for x in y) for y in self.coord.to_numpy())
            )
        f_out.close()

    def toxyz(self,outfile_name=False):
        if outfile_name == False and self.filename:
            if "." in self.filename:
                file_out = str(self.filename.split('.')[1]).strip('\\')+'.xyz'
        elif outfile_name:
            file_out = str(outfile_name)
        else:
            file_out = 'structure.xyz'
        f_out = open(file_out, 'w')
        f_out.truncate()
        content_df = self.coord.applymap(lambda x:'%0.6f' %x)
        f_out.write(
        str(len(self.coord))+'\n'+
        'Generated by toolkit@yao'+'\n'
        +'\n'.join('\t'.join(str(x) for x in y) for y in content_df.reset_index().to_numpy())
        )
        f_out.close()
        
class Analysis(Writefile):
    def __init__(self, filename=False, coord=False, lattice=False):
        super().__init__(filename, coord, lattice)
    
    def density(self,interval=0.2,axis='z',element='H') -> pd.Series:
        lattice = self.lattice
        coord = self.coord.loc[[element],axis].sort_values()
        min,max = np.floor(coord.min()),np.ceil(coord.max())
        num_dic = {}
        for j in np.arange(min,max,interval):
            count = 0
            for i in coord:
                if i>j and i<(j+interval):
                    coord.drop(i)
                    count += 1 
                    num_dic[str(np.around(j,3))+'~'+str(np.around(j+interval,3))] = count
                else:
                    continue
        axis_dic = {'x':[1,2],'y':[0,2],'z':[0,1]}
        area = np.linalg.norm(np.cross(lattice[axis_dic[axis][0]],lattice[axis_dic[axis][1]]))
        density_s = pd.Series(num_dic)/area/interval
        return density_s

    def angle(self,period=False,surface=True,main_ele='O',sub_ele='H',cutoff=1.1) -> pd.DataFrame:
        def list_gen(x,y,z):
            alllist = []
            for max in (x,y,z):
                i = 0
                list = []
                while (i < max):
                    list.append(i+1)
                    i+=1
                alllist.append(list)
            return alllist
        def df_area(df,a,b):
            for i in range(len(period)):
                df = df[(df[period[i]]>a) & (df[period[i]]<b)]               
            return df
        def vec_angle(vec1,vec2):
            a=np.inner(np.array(vec1),np.array(vec2))
            b=np.linalg.norm(np.array(vec1))*np.linalg.norm(np.array(vec2))
            return np.degrees(np.arccos(a/b))
        def surface_angle(angle):
            if angle < 90:
                angle = 90 - angle
            elif angle > 90:
                angle = angle - 90
            return angle
        org_coord = self.coord.loc[[main_ele,sub_ele]]
        org_lattice = self.lattice
        dic_xyz={'xy':(3,3,1),'yz':(1,3,3),'xz':(3,1,3),'xyz':(3,3,3)}
        if period:
            x,y,z = dic_xyz[period]
            extend_list = np.array(list(itertools.product(list_gen(x,y,z)[0],list_gen(x,y,z)[1],list_gen(x,y,z)[2])))
            ext_coord = pd.DataFrame()
            ext_lattice = org_lattice * [x,y,z]
            for i in extend_list:
                con_coord = org_coord + np.dot((i-[1,1,1]),org_lattice)#important
                ext_coord = pd.concat([ext_coord,con_coord])
            ext_coord = ext_coord.dot(np.linalg.inv(ext_lattice))
            ext_coord.columns=['x','y','z']
            main_coord,sub_coord = ext_coord.loc[main_ele],ext_coord.loc[sub_ele]
            main_np = np.dot(df_area(main_coord,1/3,2/3).to_numpy(),ext_lattice)
            sub_np = np.dot(df_area(sub_coord,0.3,0.7).to_numpy(),ext_lattice)
            if len(main_np) == len(org_coord.loc[[main_ele]]):
                pass
            else:
                raise RuntimeError('Error when extract the main element')
            angle_dic1,angle_dic2 = \
                {'distence1':[],'distence2':[],'angle1':[],'angle2':[],'angle3':[],'main_ele_pos':[]}, \
                {'main_ele_pos':[],'sub_ele_pos':[]}
            # def dist():
            for i in main_np:
                dist = spatial.distance.cdist(sub_np,[i],'euclidean')
                dist = pd.DataFrame(dist,columns=['distence'])
                dist = dist[dist['distence']<cutoff]
                if len(dist.index) == 2:
                    angle_dic1['distence1'].append(dist.iloc[0][0])
                    angle_dic1['distence2'].append(dist.iloc[1][0])
                    angle1 = vec_angle(sub_np[dist.index[0]]-i,sub_np[dist.index[1]]-i) #angle of sub-main-sub element
                    angle_dic1['angle1'].append(angle1)
                    if surface:
                        surface_dic = {'xy':[0,0,1],'yz':[0,1,0],'xz':(1,0,0)}
                        angle2 = surface_angle(vec_angle(sub_np[dist.index[0]]-i,surface_dic[period])) #angle of bond and surface
                        angle_dic1['angle2'].append(angle2)
                        angle3 = surface_angle(vec_angle((sub_np[dist.index[0]]+sub_np[dist.index[1]])/2-i,surface_dic[period])) #angle of mid_vector and surface
                        angle_dic1['angle3'].append(angle3)
                    angle_dic1['main_ele_pos'].append(i)
                elif len(dist.index) != 2:
                    angle_dic2['main_ele_pos'].append(i)
                    angle_dic2['sub_ele_pos'].append([sub_np[dist.index]])
                sub_np = np.delete(sub_np,dist.index,0)#delete the element already selected
        elif not period:
            print ('no period wiil be supported soon, use xyz instead')
        info_normal,info_abnormal = pd.DataFrame(angle_dic1),pd.DataFrame(angle_dic2)
        info_normal.columns=['Distence1(Å)','Distence2(Å)','Bond(deg)','Bond_Surface(deg)','Mid_Surface(deg)','Main_position']
        info_abnormal.columns=['Main_position','Sub_position']
        return info_normal,info_abnormal



if __name__ == '__main__':
    lattice_input=(13.00000000,0.00000000,0.00000000,-6.50000000,11.25833025,0.00000000,0.00000000,0.00000000,16.13999900)
    #ext_lattice, ext_coord = Structure(sys.argv[1],lattice=lattice_input).extend(1,2,2)
    #Writefile(coord=ext_coord,lattice=ext_lattice).tovasp(direct=True)
    #print (org_lattice, ext_coord.sort_index(),a)
    #lattice=(13.00000000,0.00000000,0.00000000,-6.50000000,11.25833025,0.00000000,0.00000000,0.00000000,16.13999900)
    a, b =Analysis(sys.argv[1],lattice=lattice_input).angle(period='xy')    
    print (a)