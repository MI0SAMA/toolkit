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
        elif self.filetype == 'xyz':
            if self.lattice_input:
                if len(self.lattice_input) == 6:
                    A,B,C,alpha,beta,gamma=self.lattice_input
                    lattice = self.cell_p2v(A,B,C,alpha,beta,gamma) 
                elif len(self.lattice_input) == 9:
                    lattice = np.zeros((3,3))
                    for i in range(3):
                        lattice[i][:] = self.lattice_input[3*i:3*i+3]
            elif not self.lattice_input:
                raise ValueError('Pleace input the lattice information for .xyz file')
        else:
            raise ValueError('Input structure file type wrong')
        lattice[np.logical_and(lattice<0.0001,lattice>0)]=0
        return lattice

    """
    USAGE
    get absolute coordinate as pdframe
    """
    def coord(self,direct=False):
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
        if direct == True:
            coord_d = np.dot(coord,(np.linalg.inv(self.lattice())))
            coord_d = pd.DataFrame(coord_d, columns=['x','y','z'], index=coord.index)
            return coord_d
        elif direct == False:
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
    USAGE
    """
    def move(self,num_inp=False,dis_inp=False,direct=False):
        org_coord = self.coord(direct=direct)
        num, dis = [], []
        #input information
        if num_inp == False:
            num_inp = input("Input the atom number to move, eg:(1 2-10): ")  
        elif isinstance(num_inp,str):
            pass
        else:
            print ('Please correct the input atom numbers')

        if dis_inp == False:
            dis_inp = input("Input the distance to move, eg:(1 1 1): ")  
        elif isinstance(dis_inp,str):
            pass
        else:
            print ('Please correct the input atom numbers')
        #moved atoms number and move distance 
        for i in num_inp.split(' '):
            if '-' in i:
                for i in range(int(i.split('-')[0]),int(i.split('-')[1])+1):
                    num.append(int(i)-1)
            else:
                num.append(int(i)-1)
        for i in dis_inp.split(' '):
            dis.append(float(i))
        #move atoms
        move_coord = org_coord
        for i in num:
            move_coord.iloc[i] = org_coord.iloc[i]+dis
        if direct:
            move_coord = np.dot(move_coord,self.lattice())
            move_coord = pd.DataFrame(move_coord, columns=['x','y','z'], index=org_coord.index)
        elif not direct:
            pass
        return move_coord
"""    
    def joint(self,coord1,coord2,lattice1,direction='x'):
        joint_axis = input("Input the joint axisas x, y, z: ")
        coord1[joint_axis] = coord1[joint_axis].apply(lambda x: x/2)
        coord2[joint_axis] = coord2[joint_axis].apply(lambda x: x/2+0.5)
        coord_new = pd.concat([coord1,coord2])
        dic= {'x':0, 'y':1, 'z':2}
        lattice1 = structure1.lattice()
        lattice2 = structure1.lattice()
        lattice_new = lattice1
        lattice_new[dic[joint_axis]] = lattice1[dic[joint_axis]]+lattice2[dic[joint_axis]]
        coord_newc = np.dot(coord_new,lattice_new)
        coord_newc = pd.DataFrame(coord_newc, columns=['x','y','z'], index=coord_new.index)
        structure = core.Writefile(coord=coord_newc,lattice=lattice_new)
        structure.tovasp(outfile_name='POSCAR_joint',direct=True)
"""

"""
Use coord lattice or structure to generate POSCAR or xyz file
"""
class Writefile(object):
    def __init__(self, filename=False, coord=False, lattice=False):
        if filename and coord==False and lattice==False:
            self._filename = str(filename)
            self._coord = Structure(filename).coord()
            self._lattice = Structure(filename).lattice()
        elif filename and lattice and coord==False:
            self._filename = str(filename)
            self._coord = Structure(filename).coord()
            if isinstance(lattice,np.ndarray):
                self._lattice = lattice
            else:
                self._lattice = Structure(filename,lattice=lattice).lattice()
        elif filename==False:
            self._coord = coord
            if isinstance(lattice,np.ndarray):
                self._lattice = lattice
            else:
                self._lattice = Structure(filename,lattice=lattice).lattice()    
        else:
            raise ValueError('Input sturcture information wrong')
    
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
            ele_dic[i] = len(self.coord.loc[[i]])
        if not outfile_name:
            file_out = 'POSCAR'
        elif outfile_name:
            file_out = str(outfile_name)
        f_out = open(file_out, 'w')
        f_out.truncate()
        f_out.write(
        'Generated by toolkit@yao'+'\n'
        +'1.00000000000000'+'\n'
        +'\n'.join('\t'.join('%0.6f' %x for x in y) for y in self.lattice)+'\n'
        +' '.join(ele_dic.keys())+'\n'
        +' '.join(str(i) for i in ele_dic.values())+'\n'
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
        
    # def tolmp(self,outfile_name=False):
    #     if outfile_name == False and self.filename:
    #         if "." in self.filename:
    #             file_out = str(self.filename.split('.')[1]).strip('\\')+'.lmp'
    #     elif outfile_name:
    #         file_out = str(outfile_name)
    #     f_out = open(file_out, 'w')
    #     f_out.truncate()    

class Analysis(Writefile):
    def __init__(self, filename=False, coord=False, lattice=False):
        super().__init__(filename, coord, lattice)
    
    def density(self,interval=0.2,min=False,max=False,axis='z',element='H') -> pd.Series:
        lattice = self.lattice
        coord_one = self.coord.loc[[element],axis].sort_values()
        if min == False and max == False:
            min,max = np.floor(coord_one.min()),np.ceil(coord_one.max())
        else:
            pass
        num_dic = {}
        for j in np.arange(min,max,interval):
            count = 0
            num_dic[str(np.around(j,3))+'~'+str(np.around(j+interval,3))] = count
            for i in coord_one:
                if i>j and i<(j+interval):
                    coord_one.drop(i)#have bug 
                    count += 1 
                    num_dic[str(np.around(j,3))+'~'+str(np.around(j+interval,3))] = count
                else:
                    continue
        density_c = pd.Series(num_dic)
        axis_dic = {'x':[1,2],'y':[0,2],'z':[0,1]}
        area = np.linalg.norm(np.cross(lattice[axis_dic[axis][0]],lattice[axis_dic[axis][1]]))
        density_d = density_c/area/interval
        return density_c,density_d
        
    def angle(self,period=False,main_ele='O',sub_ele='H',cutoff=1.1) -> pd.DataFrame:
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
            return 90 - angle
        def del_empt(dic):
            for key in list(dic.keys()):
                if dic.get(key) == []:
                    del dic[key]
            return dic
        def half_coord(i,period):
            dic_half={'xy':(1/2,1/2,1),'yz':(1,1/2,1/2),'xz':(1/2,1,1/2),'xyz':(1/2,1/2,1/2)}
            i = np.dot(i,(np.linalg.inv(ext_lattice)))
            i = i*np.array([dic_half[period]])
            i = np.dot(i,ext_lattice)
            return i
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
                {'Distence1(A)':[],'Distence2(A)':[],'Bond(deg)':[],'Bond_Surface(deg)':[],'Mid_Surface(deg)':[],'Main_position':[]}, \
                {'Main_position':[],'Sub_position':[]}
            # def dist():
            for i in main_np:
                dist = spatial.distance.cdist(sub_np,[i],'euclidean')
                dist = pd.DataFrame(dist,columns=['distence'])
                dist = dist[dist['distence']<cutoff]
                if len(dist.index) == 2:
                    angle_dic1['Distence1(A)'].append(dist.iloc[0][0])
                    angle_dic1['Distence2(A)'].append(dist.iloc[1][0])
                    angle1 = vec_angle(sub_np[dist.index[0]]-i,sub_np[dist.index[1]]-i) #angle of sub-main-sub element
                    angle_dic1['Bond(deg)'].append(angle1)
                    if period != 'xyz':
                        surface_dic = {'xy':[0,0,1],'yz':[0,1,0],'xz':(1,0,0)}
                        angle2 = surface_angle(vec_angle(sub_np[dist.index[0]]-i,surface_dic[period])) #angle of bond and surface
                        angle_dic1['Bond_Surface(deg)'].append(angle2)
                        angle3 = surface_angle(vec_angle((sub_np[dist.index[0]]+sub_np[dist.index[1]])/2-i,surface_dic[period])) #angle of mid_vector and surface
                        angle_dic1['Mid_Surface(deg)'].append(angle3)
                    i = half_coord(i,period)
                    angle_dic1['Main_position'].append(i)
                elif len(dist.index) != 2:
                    i = half_coord(i,period)
                    angle_dic2['Main_position'].append(i)
                    for j in sub_np[dist.index]:
                        save_list=[]
                        j = half_coord(j,period)
                        save_list.append(j)
                    angle_dic2['Sub_position'].append(save_list)
                sub_np = np.delete(sub_np,dist.index,0)#delete the element already selected
            if period == 'xyz':
                angle_dic1 = del_empt(angle_dic1)
        elif not period:
            print ('Pleas input the structure periodicity')
        info_normal,info_abnormal = pd.DataFrame(angle_dic1),pd.DataFrame(angle_dic2)
        info_normal.iloc[:,-1] = info_normal.iloc[:,-1].map((lambda x: x[0]))
        info_abnormal.iloc[:,0] = info_abnormal.iloc[:,0].map((lambda x: x[0]))
        info_abnormal.iloc[:,1] = info_abnormal.iloc[:,1].map((lambda x: x[0]))
        return info_normal,info_abnormal



if __name__ == '__main__':
    lattice_input=(13.00000000,0.00000000,0.00000000,-6.50000000,11.25833025,0.00000000,0.00000000,0.00000000,16.13999900)
    # structure = Structure(sys.argv[1])
    # move_coord = structure.move(num_inp='1 23-24',dis_inp='1 1 1')
    # print(structure.coord()-move_coord)
    # print(move_coord)
    #ext_lattice, ext_coord = Structure(sys.argv[1],lattice=lattice_input).extend(1,2,2)