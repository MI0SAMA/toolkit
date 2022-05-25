#!/usr/bin/env python
import numpy as np 
import pandas as pd
import sys 

#read information from stucture file
class Structure(object):
    #initial
    def __init__(self,filename):
        self._filename = str(filename)
        if "." in self._filename:
            self._filetype = self._filename.split('.')[-1].lower()
        elif 'poscar' in self._filename.lower():
            self._filetype = 'vasp'
        with open(self._filename,'r') as f:
            self._filecontent = [x.rstrip().split() for x in f]

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
    def lattice(self,*args):#get lattice as a nparray
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
        elif self.filetype == 'xyz' and args != ():
            if len(args) == 6:
                A,B,C,alpha,beta,gamma=args
                lattice = self.cell_p2v(A,B,C,alpha,beta,gamma) 
            elif len(args) == 9:
                lattice = np.zeros((3,3))
                for i in range(3):
                    lattice[i][:] = args[3*i:3*i+3]
        else:
            raise ValueError('Input of lattice information wrong')
        lattice[np.logical_and(lattice<0.0001,lattice>0)]=0
        return lattice

    def coord(self):#get coordinate as pdframe
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

# def main():
#     a=Structure(sys.argv[1])
#     print(a.coord().type())

# if __name__ == '__main__':
#     main()