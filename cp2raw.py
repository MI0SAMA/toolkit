#!/usr/bin/env python
import os

outfile = input("Input the outfile name, Eneter for default: cp2k.out ")
if outfile=="":
    outfile = "cp2k.out"
else:
    pass
with open(outfile,'r') as f:
    for x in f:
        line = x.rstrip().split()
        if 'Project' in line and 'name' in line:
            project_name = line[-1]

step_inp = input("Input the chosen step: eg(2000 3000), Enter for defalut: all ")
if step_inp=="":
    step_end = len(open(project_name + '-1.ener','r').readlines())
    step = [0,step_end-1]
else:
    step_stat, step_end = step_inp.rstrip().split()
    step = [int(step_stat), int(step_end)]

def type(project_name):
    count = 0
    element = []
    with open (project_name + '-pos-1.xyz') as f:
        for x in f:
            line = x.rstrip().split()
            if count == 0:
                element_num = int(line[0])
                count += 1
            elif count == element_num+2:
                break
            else:
                if len(line) != 4:
                    count += 1
                if len(line) == 4:
                    element.append(line[0])
                    count += 1
    element_type = list(dict.fromkeys(element))
    element_type = dict(zip(element_type,(list(range(len(element_type))))))
    return element_type,element_num,element

def box():
    box = []
    with open(outfile,'r') as f:
        for x in f:
            line = x.rstrip().split()
            if 'CELL_TOP|' in line:
                if 'Vector' in line:
                    for i in [4,5,6]:
                        box.append(float(line[i]))
                if 'Angle' in line:
                    break
    return box

def energy(project_name):
    energy = []
    with open (project_name + '-1.ener') as f:
        for x in f:
            line = x.rstrip().split()
            try:
                energy.append(float(line[-2])*27.211386)
            except ValueError:
                pass
    return energy

#force = [[x1,y1,x2,y2],[x1,y1,x2,y2]]
def force(project_name,element_num): 
    force = []
    def con(a):
        a = float(a)*27.211386/0.529177
        return a
    if os.path.isfile(project_name+'-force-1.xyz'):
        pass
    else:
        forcefile = input("Input the force file name: ")
        print(os.system(
        'cd '+forcefile+' && cat *-force-1_* > ../'+project_name+'-force-1.xyz && cd ..'))
    with open (project_name + '-force-1.xyz') as f_in:
        i = 0
        force1 = []
        for x in f_in:
            line = x.rstrip().split()
            if len(line) == 6:
                force1.extend((con(line[3]),con(line[4]),con(line[5])))
                i += 1
                if i == element_num:
                    force.append(force1)
                    force1 = []
                    i = 0
    return (force)
                
def coord(project_name,element_num):
    def con(a):
        return (float(a))
    coord = []
    with open (project_name + '-pos-1.xyz') as f_in:
        i=0
        coord1=[]
        for x in f_in:
            line = x.rstrip().split()
            if len(line) == 4:
                coord1.extend((con(line[1]),con(line[2]),con(line[3])))
                i += 1
                if i == element_num:
                    coord.append(coord1)
                    i = 0 
                    coord1=[]      
    return coord

def all2raw(box,energy,coord,force,step,element_type,element):
    if len(energy)==len(force)==len(coord):
        step_num = step[1]-step[0]
        f_box,f_ene,f_for,f_coo=open('box.raw','w'),open('energy.raw','w'),open('force.raw','w'),open('coord.raw','w')
        for i in range(step_num):
            print (str(step[0]+i))
            f_box.write(' '.join(str(j) for j in box)+'\n')
            f_ene.write(str(energy[step[0]+i])+'\n')
            f_for.write(' '.join(str(j) for j in force[step[0]+i])+'\n')
            f_coo.write(' '.join(str(j) for j in coord[step[0]+i])+'\n')
        f1 = open('type.raw','w')
        f1.truncate()
        for i in element:
            f1.write(str(element_type[i])+'\n')
        f1.close()
        f2 = open('type_map.raw','w')
        f2.truncate()
        for i in element_type.keys():
            f2.write(i+'\n')
        f2.close()
        f_box.close()
        f_ene.close()
        f_for.close()
        f_coo.close()
    else:
        print ("Somthing wrong when reading the information")
        print(len(energy),len(force),len(coord))

if __name__ == '__main__':
    element_type, element_num, element= type(project_name)
    box,energy,force,coord=box(),energy(project_name),force(project_name,element_num),coord(project_name,element_num)
    all2raw(box,energy,coord,force,step,element_type,element)
    print(os.system("mkdir raw && mv *.raw raw/"))