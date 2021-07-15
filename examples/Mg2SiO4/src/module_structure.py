import numpy as np
import time
import random
import shutil
import math
import yaml    
import subprocess
from subprocess import check_output
import os
from os.path import isfile, join

from module_crossover import *

pi=3.14159265358979
rad2deg=180.0/pi
 
def asin(a):
        return math.atan2(a,math.sqrt(1.0-a*a))
def acos(a):
        return pi/2.0-asin(a)
def norm(x):
        return (math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))
def dotprod (x,y):
        return ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] )
def angle (v1,v2):
        myacos = dotprod(v1,v2)/norm(v1)/norm(v2)
        if myacos>1.0 :
            myacos = 1.0
        if myacos<-1.0:
            myacos = -1.0
        return(acos(myacos)*180.0/3.14159265358979)
        gc.collect()

def random_gen(Vmin, inp_file, atomnumlist, atomnamelist, tolerance_matrix, filename, gen, pop):

    t1 = time.time()

    sub_lattice_ratio = inp_file['random_condition']['sublattice_generation'] 
    max_sub_lattice = inp_file['random_condition']['max_sub_lattice'] 

    # sublattice generation

    atomnumlist_temp = [atomnumlist[i] for i in range(len(atomnumlist))]

    if random.random() < sub_lattice_ratio:
        sub_lattice_gen = True 
    else:
        sub_lattice_gen = False

    if sub_lattice_gen == True:

        if max_sub_lattice == 4:
            if random.random() < 0.5:
                for i in range(len(atomnumlist)):
                    atomnumlist_temp[i] = int(atomnumlist[i] / 4)
                Vmin = float(Vmin / 4.0)
                sublattice = 4
            else:
                for i in range(len(atomnumlist)):
                    atomnumlist_temp[i] = int(atomnumlist[i] / 2)
                Vmin = float(Vmin / 2.0)
                sublattice = 2
 
        if max_sub_lattice == 2:
            for i in range(len(atomnumlist)):
                atomnumlist_temp[i] = int(atomnumlist[i] / 2)
            Vmin = float(Vmin / 2.0)
            sublattice = 2

    # extract features
    wyckoff =  inp_file['random_condition']['force_general_Wyckoff_site']
    attempts     = inp_file['random_condition']['maximum_attempts_for_one_space_group_and_volume']
    scale_factor = inp_file['random_condition']['scale_factor']
    
    # make tolerance array
    atomtypenum = len(atomnumlist)
    tolerance   = [] 
    for i in range(1,atomtypenum+1):
        for j in range(i,atomtypenum+1):
            tolerance.append([[i,j],tolerance_matrix[i-1][j-1]])
    
    for i in range(len(tolerance)):
        tolerance[i][1] = tolerance[i][1] * scale_factor

    # preprocessing
    atom_array = []
    for i in range(len(atomnumlist_temp)):
        atom_array += [i+1]*atomnumlist_temp[i]

    # generate random_structure
    temp = 0
    while 1:

        V = (random.random()*0.6+0.7)*Vmin

        lmin = V**(1.0/3.0) * 0.3
        lmax = V**(1.0/3.0) * 3.0

        temp += 1
        spg = random.randint(1,230)
        import pyrandspg

        pymin = pyrandspg.LatticeStruct(lmin, lmin, lmin, 60.0, 60.0, 60.0)
        pymax = pyrandspg.LatticeStruct(lmax, lmax, lmax, 120.0, 120.0, 120.0)

        input_ = pyrandspg.RandSpgInput(spg, atom_array, pymin, pymax, 0.0, V*0.95, V*1.05, attempts, tolerance, wyckoff)

        crystal = pyrandspg.RandSpg.randSpgCrystal(input_)

        structure_= crystal.getPOSCARString()
        del pyrandspg
        
        if 'nan' not in structure_:
            del crystal, input_, pymin, pymax
            break
        else:
            del spg, crystal, input_, pymax, pymin, structure_

    structure = structure_.split('\n')
    elementorder = list(map(int,structure[5].split()))
    elementnumber= list(map(int,structure[6].split()))
    atomtypenum  = len(elementnumber)

    written_array = []
    for i in range(len(elementnumber)):
        written_array += [elementorder[i]]*elementnumber[i]

    written_coor = []
    tot_atom_num = sum(atomnumlist_temp)
    for i in range(tot_atom_num):
        written_coor.append(list(map(float,structure[8+i].split())))

    coor_dir = []
    for i in range(len(atomnumlist_temp)):
        for j in range(tot_atom_num):
            if written_array[j] == i+1:
                coor_dir.append(written_coor[j])

    latt = []
    for i in range(2,5):
        latt.append(list(map(float,structure[i].split())))

    latt     = np.array(latt)
    coor_dir = np.array(coor_dir)

    coor = coor_dir.dot(latt)
    
    # make supercell
    if sub_lattice_gen == True:
        coor,latt = make_into_supercell(coor,latt,sublattice)

    atomtype = []
    for i in range(len(atomnumlist)):
        atomtype += [i+1]*atomnumlist[i]

    # convert into lammps format
    p_a= math.sqrt(latt[0][0]**2.0 + latt[0][1]**2.0 + latt[0][2]**2.0);
    p_b= math.sqrt(latt[1][0]**2.0 + latt[1][1]**2.0 + latt[1][2]**2.0);
    p_c= math.sqrt(latt[2][0]**2.0 + latt[2][1]**2.0 + latt[2][2]**2.0);
    alpha=  angle(latt[1],latt[2]);  beta=  angle(latt[0],latt[2]); gamma=  angle(latt[0],latt[1]);  # Angles in degree
    alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians

    lx=   p_a
    p_xy= p_b * math.cos(gammar)
    p_xz= p_c * math.cos(betar)
    ly=   math.sqrt(p_b**2.0 - p_xy**2.0)
    p_yz= (p_b*p_c*math.cos(alphar)-p_xy*p_xz)/(ly)
    lz=   math.sqrt(p_c**2.0 - p_xz**2.0 - p_yz**2.0)

    with open("coo_"+filename,"w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write(str(sum(atomnumlist))+" atoms\n")
        fw.write(str(len(atomnumlist))+" atom types\n")
        fw.write("\n")
        fw.write("0.000000 %10.6f xlo xhi\n"%(lx))
        fw.write("0.000000 %10.6f ylo yhi\n"%(ly))
        fw.write("0.000000 %10.6f zlo zhi\n"%(lz))
        fw.write("\n")
        fw.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy, p_xz, p_yz))
        fw.write("\n")
        fw.write("Atoms\n")
        fw.write("\n")
        for i in range(sum(atomnumlist)):
            fw.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, atomtype[i], coor[i][0], coor[i][1], coor[i][2]))

    #spring
    tot_atom_num = sum(atomnumlist)

    del sub_lattice_ratio, max_sub_lattice, wyckoff, attempts, scale_factor, atomtypenum, tolerance, lmin, lmax, atom_array, structure, structure_, tot_atom_num

    del written_array, written_coor, coor_dir, latt, atomtype, elementorder, elementnumber,atomnumlist_temp 

    t2 = time.time()

    with open("../random_structure_log","a") as f:
        f.write("  %d    %d      %.2f s    %d        %d\n"%(gen, pop, t2-t1, temp, spg))

    del spg 

    return coor, [lx,ly,lz,p_xy,p_xz,p_yz]

    gc.collect()
    del coor, temp

def con_permutation(parent, output_path,parent_data, atomnamelist, atomnumlist, inp_file):

    tot_atom_num = sum(atomnumlist)
 
    coor = np.array(parent_data['coordinate'])
    lattold = np.array(parent_data['lattice'])

    latt = [[0.0,0.0,0.0] for i in range(3)]
    latt = np.array(latt)
    latt[0][0] = lattold[0]
    latt[1][1] = lattold[1]
    latt[2][2] = lattold[2]
    latt[1][0] = lattold[3]
    latt[2][0] = lattold[4]
    latt[2][1] = lattold[5]

    # making number dictionary : example) for In4Te8O20F : atom_to_num['In'] = 0, atom_to_num['Te'] = 1, atom_to_num['O'] = 2, ...
    atom_to_num = {}
    for i in range(len(atomnumlist)):
        atom_to_num[atomnamelist[i]] = i

    # making group : example)  [[0,1,2,3],[4,5,6,7],[8,9]] for A4B4C2
    group = []
    group.append([i for i in range(0,atomnumlist[0])])
    for j in range(1,len(atomnumlist)-1):
        group.append([i for i in range(sum(atomnumlist[0:j]),sum(atomnumlist[0:j+1]))])
    group.append([i for i in range(sum(atomnumlist[0:-1]),sum(atomnumlist))])

    # making permutation group lists
    permutation_groups = []
    if 'material' not in inp_file:
        permutation_groups.append(group)

    else:
        for permutationgroup in inp_file['permutation_allow']:
            atomarray = inp_file['permutation_allow'][permutationgroup]
            #print (atom_to_num[atomarray[1]])
            permutation_groups.append([group[atom_to_num[atomarray[j]]] for j in range(len(atomarray))])

    # pick what group to change
    random_choice_number = random.randint(1,2**len(permutation_groups)-1)
    random_choice = str(bin(random_choice_number))
    random_choice = str(bin(random_choice_number))[2:len(str(bin(random_choice_number)))]
    random_choice = [0]*(len(permutation_groups)-len(random_choice)) + list(map(int,random_choice))

    # change atoms
    for i in range(len(permutation_groups)):
        if random_choice[i] == 1:
            # pick subgroups to change
            change_groups = random.sample([j for j in range(len(permutation_groups[i]))],2)
            # pick atoms to change
            atomnum1 =  random.sample(permutation_groups[i][change_groups[0]],1)[0]
            atomnum2 =  random.sample(permutation_groups[i][change_groups[1]],1)[0]

            # exchange atoms
            temp = [0,0,0]
            for j in range(3):
                temp[j] = coor[atomnum1][j]
            for j in range(3):
                coor[atomnum1][j] = coor[atomnum2][j]
            for j in range(3):
                coor[atomnum2][j] = temp[j]
            del temp

    # convert into lammps format
    p_a= math.sqrt(latt[0][0]**2.0 + latt[0][1]**2.0 + latt[0][2]**2.0);
    p_b= math.sqrt(latt[1][0]**2.0 + latt[1][1]**2.0 + latt[1][2]**2.0);
    p_c= math.sqrt(latt[2][0]**2.0 + latt[2][1]**2.0 + latt[2][2]**2.0);
    alpha=  angle(latt[1],latt[2]);  beta=  angle(latt[0],latt[2]); gamma=  angle(latt[0],latt[1]);  # Angles in degree
    alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians

    lx=   p_a
    p_xy= p_b * math.cos(gammar)
    p_xz= p_c * math.cos(betar)
    ly=   math.sqrt(p_b**2.0 - p_xy**2.0)
    p_yz= (p_b*p_c*math.cos(alphar)-p_xy*p_xz)/(ly)
    lz=   math.sqrt(p_c**2.0 - p_xz**2.0 - p_yz**2.0)

    latt = [lx,ly,lz,p_xy,p_xz,p_yz]
    atomtype = []
    for i in range(len(atomnumlist)):
        atomtype += [i+1]*atomnumlist[i]

    with open("coo","w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write(str(sum(atomnumlist))+" atoms\n")
        fw.write(str(len(atomnumlist))+" atom types\n")
        fw.write("\n")
        fw.write("0.000000 %10.6f xlo xhi\n"%(lx))
        fw.write("0.000000 %10.6f ylo yhi\n"%(ly))
        fw.write("0.000000 %10.6f zlo zhi\n"%(lz))
        fw.write("\n")
        fw.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy, p_xz, p_yz))
        fw.write("\n")
        fw.write("Atoms\n")
        fw.write("\n")
        for i in range(sum(atomnumlist)):
            fw.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, atomtype[i], coor[i][0], coor[i][1], coor[i][2]))


    del atom_to_num, atomtype, permutation_groups, random_choice, change_groups, permutationgroup, lattold
    return coor, latt 
    gc.collect()

def all_permutation(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file):

    tot_atom_num = sum(atomnumlist)    
    
    coor = np.array(parent_data['coordinate'])
    lattold = np.array(parent_data['lattice'])

    latt = [[0.0,0.0,0.0] for i in range(3)]
    latt = np.array(latt)
    latt[0][0] = lattold[0]
    latt[1][1] = lattold[1]
    latt[2][2] = lattold[2]
    latt[1][0] = lattold[3]
    latt[2][0] = lattold[4]
    latt[2][1] = lattold[5]

    # making number dictionary : example) for In4Te8O20F : atom_to_num['In'] = 0, atom_to_num['Te'] = 1, atom_to_num['O'] = 2, ...
    atom_to_num = {}
    for i in range(len(atomnumlist)):
        atom_to_num[atomnamelist[i]] = i

    # making group : example)  [[0,1,2,3],[4,5,6,7],[8,9]] for A4B4C2
    group = []
    group.append([i for i in range(0,atomnumlist[0])])
    for j in range(1,len(atomnumlist)-1):
        group.append([i for i in range(sum(atomnumlist[0:j]),sum(atomnumlist[0:j+1]))])
    group.append([i for i in range(sum(atomnumlist[0:-1]),sum(atomnumlist))])

    # making permutation group lists
    permutation_groups = []
    permutation_groups.append(group)

    # pick what group to change
    random_choice_number = random.randint(1,2**len(permutation_groups)-1)
    random_choice = str(bin(random_choice_number))
    random_choice = str(bin(random_choice_number))[2:len(str(bin(random_choice_number)))]
    random_choice = [0]*(len(permutation_groups)-len(random_choice)) + list(map(int,random_choice))

    # change atoms

    for i in range(len(permutation_groups)):
        if random_choice[i] == 1:
            # pick subgroups to change
            change_groups = random.sample([j for j in range(len(permutation_groups[i]))],2)

            # decide permute number
            max_num = min(len(permutation_groups[i][change_groups[0]]),len(permutation_groups[i][change_groups[1]]))
            max_num = math.floor(max_num/2)

            if max_num < 1:
                max_num = 1
            how_many_permut = random.randint(1,max_num)

            atomnum1s =  random.sample(permutation_groups[i][change_groups[0]],how_many_permut)
            atomnum2s =  random.sample(permutation_groups[i][change_groups[1]],how_many_permut)

            # exchange atoms
            for k in range(how_many_permut):

                atomnum1 = atomnum1s[k]
                atomnum2 = atomnum2s[k]
    
                #with open("test","a") as fw:
                #    fw.write(str(atomnum1))
                #    fw.write("\n")
                #    fw.write(str(atomnum2))
                #    fw.write("\n")

                temp = [0,0,0]
                for j in range(3):
                    temp[j] = coor[atomnum1][j]
                for j in range(3):
                    coor[atomnum1][j] = coor[atomnum2][j]
                for j in range(3):
                    coor[atomnum2][j] = temp[j]

    # convert into lammps format
    p_a= math.sqrt(latt[0][0]**2.0 + latt[0][1]**2.0 + latt[0][2]**2.0);
    p_b= math.sqrt(latt[1][0]**2.0 + latt[1][1]**2.0 + latt[1][2]**2.0);
    p_c= math.sqrt(latt[2][0]**2.0 + latt[2][1]**2.0 + latt[2][2]**2.0);
    alpha=  angle(latt[1],latt[2]);  beta=  angle(latt[0],latt[2]); gamma=  angle(latt[0],latt[1]);  # Angles in degree
    alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians

    lx=   p_a
    p_xy= p_b * math.cos(gammar)
    p_xz= p_c * math.cos(betar)
    ly=   math.sqrt(p_b**2.0 - p_xy**2.0)
    p_yz= (p_b*p_c*math.cos(alphar)-p_xy*p_xz)/(ly)
    lz=   math.sqrt(p_c**2.0 - p_xz**2.0 - p_yz**2.0)

    latt = [lx,ly,lz,p_xy,p_xz,p_yz]

    atomtype = []
    for i in range(len(atomnumlist)):
        atomtype += [i+1]*atomnumlist[i]

    with open("coo","w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write(str(sum(atomnumlist))+" atoms\n")
        fw.write(str(len(atomnumlist))+" atom types\n")
        fw.write("\n")
        fw.write("0.000000 %10.6f xlo xhi\n"%(lx))
        fw.write("0.000000 %10.6f ylo yhi\n"%(ly))
        fw.write("0.000000 %10.6f zlo zhi\n"%(lz))
        fw.write("\n")
        fw.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy, p_xz, p_yz))
        fw.write("\n")
        fw.write("Atoms\n")
        fw.write("\n")
        for i in range(sum(atomnumlist)):
            fw.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, atomtype[i], coor[i][0], coor[i][1], coor[i][2]))
    
    del atom_to_num, atomtype, permutation_groups, random_choice, change_groups, lattold, temp
    return coor, latt 
    gc.collect()

def latticemutation(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file):
    

    tot_atom_num = sum(atomnumlist)
    
    coor = np.array(parent_data['coordinate'])
    lattold = np.array(parent_data['lattice'])

    latt = [[0.0,0.0,0.0] for i in range(3)]
    latt = np.array(latt)
    latt[0][0] = lattold[0]
    latt[1][1] = lattold[1]
    latt[2][2] = lattold[2]
    latt[1][0] = lattold[3]
    latt[2][0] = lattold[4]
    latt[2][1] = lattold[5]

    
    # conversion from cartesian to direct
    coor = coor.dot(np.linalg.inv(latt))

    # lattice mutation
    for i in range(3):
        for j in range(i,3):
            randomvalue = 100

            if i == j:
                while randomvalue > 0.5 or randomvalue < -0.5:
                    randomvalue = np.random.normal(scale=0.1)
                latt[i][i] = latt[i][i]*(1.0 + randomvalue)
            else:
                while randomvalue > 1.0 or randomvalue < -1.0:
                    randomvalue = np.random.normal(scale=0.5)
                latt[j][i] = latt[j][i] + latt[j][j]*(randomvalue/2.0)

    
    coor = coor.dot(latt)
    
    # convert into lammps format
    p_a= math.sqrt(latt[0][0]**2.0 + latt[0][1]**2.0 + latt[0][2]**2.0);
    p_b= math.sqrt(latt[1][0]**2.0 + latt[1][1]**2.0 + latt[1][2]**2.0);
    p_c= math.sqrt(latt[2][0]**2.0 + latt[2][1]**2.0 + latt[2][2]**2.0);
    alpha=  angle(latt[1],latt[2]);  beta=  angle(latt[0],latt[2]); gamma=  angle(latt[0],latt[1]);  # Angles in degree
    alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians

    lx=   p_a
    p_xy= p_b * math.cos(gammar)
    p_xz= p_c * math.cos(betar)
    ly=   math.sqrt(p_b**2.0 - p_xy**2.0)
    p_yz= (p_b*p_c*math.cos(alphar)-p_xy*p_xz)/(ly)
    lz=   math.sqrt(p_c**2.0 - p_xz**2.0 - p_yz**2.0)

    # print
    atomtype = []
    for i in range(len(atomnumlist)):
        atomtype += [i+1]*atomnumlist[i]
    latt = [lx,ly,lz,p_xy,p_xz,p_yz]

    with open("coo","w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write(str(sum(atomnumlist))+" atoms\n")
        fw.write(str(len(atomnumlist))+" atom types\n")
        fw.write("\n")
        fw.write("0.000000 %10.6f xlo xhi\n"%(lx))
        fw.write("0.000000 %10.6f ylo yhi\n"%(ly))
        fw.write("0.000000 %10.6f zlo zhi\n"%(lz))
        fw.write("\n")
        fw.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy, p_xz, p_yz))
        fw.write("\n")
        fw.write("Atoms\n")
        fw.write("\n")
        for i in range(sum(atomnumlist)):
            fw.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, atomtype[i], coor[i][0], coor[i][1], coor[i][2]))
    
    del atomtype, lattold

    return coor, latt 
    gc.collect()


def kept_best(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file):
    
    tot_atom_num = sum(atomnumlist)
 
    coor = np.array(parent_data['coordinate'])
    lattold = np.array(parent_data['lattice'])

    latt = [[0.0,0.0,0.0] for i in range(3)]
    latt = np.array(latt)
    latt[0][0] = lattold[0]
    latt[1][1] = lattold[1]
    latt[2][2] = lattold[2]
    latt[1][0] = lattold[3]
    latt[2][0] = lattold[4]
    latt[2][1] = lattold[5]

    # convert into lammps format
    p_a= math.sqrt(latt[0][0]**2.0 + latt[0][1]**2.0 + latt[0][2]**2.0);
    p_b= math.sqrt(latt[1][0]**2.0 + latt[1][1]**2.0 + latt[1][2]**2.0);
    p_c= math.sqrt(latt[2][0]**2.0 + latt[2][1]**2.0 + latt[2][2]**2.0);
    alpha=  angle(latt[1],latt[2]);  beta=  angle(latt[0],latt[2]); gamma=  angle(latt[0],latt[1]);  # Angles in degree
    alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians

    lx=   p_a
    p_xy= p_b * math.cos(gammar)
    p_xz= p_c * math.cos(betar)
    ly=   math.sqrt(p_b**2.0 - p_xy**2.0)
    p_yz= (p_b*p_c*math.cos(alphar)-p_xy*p_xz)/(ly)
    lz=   math.sqrt(p_c**2.0 - p_xz**2.0 - p_yz**2.0)

    # print
    latt = [lx,ly,lz,p_xy,p_xz,p_yz]
    atomtype = []
    for i in range(len(atomnumlist)):
        atomtype += [i+1]*atomnumlist[i]
    
    with open("coo","w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write(str(sum(atomnumlist))+" atoms\n")
        fw.write(str(len(atomnumlist))+" atom types\n")
        fw.write("\n")
        fw.write("0.000000 %10.6f xlo xhi\n"%(lx))
        fw.write("0.000000 %10.6f ylo yhi\n"%(ly))
        fw.write("0.000000 %10.6f zlo zhi\n"%(lz))
        fw.write("\n")
        fw.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy, p_xz, p_yz))
        fw.write("\n")
        fw.write("Atoms\n")
        fw.write("\n")
        for i in range(sum(atomnumlist)):
            fw.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, atomtype[i], coor[i][0], coor[i][1], coor[i][2]))
    
    del atomtype, lattold
    return coor, latt 
    gc.collect()


def crossover(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file, average_atom_e, tolerance_matrix):
    
    tot_atom_num = sum(atomnumlist)
    
    #info of structure 1
    coor1 = np.array(parent_data['coordinate'][0])
    latt1old = np.array(parent_data['lattice'][0])

    latt1 = [[0.0,0.0,0.0] for i in range(3)]
    latt1 = np.array(latt1)
    latt1[0][0] = latt1old[0]
    latt1[1][1] = latt1old[1]
    latt1[2][2] = latt1old[2]
    latt1[1][0] = latt1old[3]
    latt1[2][0] = latt1old[4]
    latt1[2][1] = latt1old[5]

    ae1 = parent_data['atom_e'][0]

    #info of structure 2
    coor2 = np.array(parent_data['coordinate'][1])
    latt2old = np.array(parent_data['lattice'][1])

    latt2 = [[0.0,0.0,0.0] for i in range(3)]
    latt2 = np.array(latt2)
    latt2[0][0] = latt2old[0]
    latt2[1][1] = latt2old[1]
    latt2[2][2] = latt2old[2]
    latt2[1][0] = latt2old[3]
    latt2[2][0] = latt2old[4]
    latt2[2][1] = latt2old[5]
    
    ae2 = parent_data['atom_e'][1]
    
    # conversion from cartesian to direct
    coor1 = coor1.dot(np.linalg.inv(latt1))
    coor2 = coor2.dot(np.linalg.inv(latt2))

    ######## Genetic Algorithm! ######
    
    # making slabs
    dir1, E1, slab1 = cut_the_slab(coor1, latt1, ae1, atomnumlist, average_atom_e, -1, inp_file)
    dir2, E2, slab2 = cut_the_slab(coor2, latt2, ae2, atomnumlist, average_atom_e, -1, inp_file)

    # remake slab of higher-energy structure fit to low-energy structure
    if E1 > E2:
        latt_direction = decide_direction(latt1,latt2,dir2)
        a_temp, E_temp, slab1          = cut_the_slab(coor1, latt1, ae1, atomnumlist, average_atom_e, latt_direction, inp_file)
        dominance = 2
    else:
        latt_direction = decide_direction(latt2,latt1,dir1)
        a_temp, E_temp, slab2          = cut_the_slab(coor2, latt2, ae2, atomnumlist, average_atom_e, latt_direction, inp_file)
        dominance = 1

    # integrate slabs
    
    integrated_structure, interface = integrate_slabs(slab1,slab2,dominance,atomnumlist,atomnamelist, average_atom_e, inp_file)
    latt = integrated_structure['latt']

    # remove atoms
    coor_new, type_new, atomnumlist_difference, vacant_sites  = remove_redundant_atoms(integrated_structure, atomnamelist, atomnumlist, average_atom_e, inp_file)

#    print (coor_new)

    # add atoms
    coor_final = add_atoms(coor_new, latt, atomnumlist_difference, atomnumlist, atomnamelist, tolerance_matrix, inp_file, interface)
    
    # convert into lammps format
    p_a= math.sqrt(latt[0][0]**2.0 + latt[0][1]**2.0 + latt[0][2]**2.0);
    p_b= math.sqrt(latt[1][0]**2.0 + latt[1][1]**2.0 + latt[1][2]**2.0);
    p_c= math.sqrt(latt[2][0]**2.0 + latt[2][1]**2.0 + latt[2][2]**2.0);
    alpha=  angle(latt[1],latt[2]);  beta=  angle(latt[0],latt[2]); gamma=  angle(latt[0],latt[1]);  # Angles in degree
    alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians

    lx=   p_a
    p_xy= p_b * math.cos(gammar)
    p_xz= p_c * math.cos(betar)
    ly=   math.sqrt(p_b**2.0 - p_xy**2.0)
    p_yz= (p_b*p_c*math.cos(alphar)-p_xy*p_xz)/(ly)
    lz=   math.sqrt(p_c**2.0 - p_xz**2.0 - p_yz**2.0)

    atomtype = []
    for i in range(len(atomnumlist)):
        atomtype += [i+1]*atomnumlist[i]
    latt = [lx,ly,lz,p_xy,p_xz,p_yz]

    # SPRING CONSTANT CALCULATION
    latt, coor_final = spring(latt, coor_final, tolerance_matrix,atomnumlist,atomnamelist,atomtype,tot_atom_num)

    # write structure
    with open("coo","w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write(str(sum(atomnumlist))+" atoms\n")
        fw.write(str(len(atomnumlist))+" atom types\n")
        fw.write("\n")
        fw.write("0.000000 %10.6f xlo xhi\n"%(lx))
        fw.write("0.000000 %10.6f ylo yhi\n"%(ly))
        fw.write("0.000000 %10.6f zlo zhi\n"%(lz))
        fw.write("\n")
        fw.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy, p_xz, p_yz))
        fw.write("\n")
        fw.write("Atoms\n")
        fw.write("\n")
        for i in range(sum(atomnumlist)):
            fw.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, atomtype[i], coor_final[i][0], coor_final[i][1], coor_final[i][2]))

    return coor_final, latt 
    del coor_new, type_new, atomnumlist_difference, vacant_sites, integrated_structure, interface, slab1, slab2
    del coor1, coor2, latt1, latt2, ae1, ae2, latt1old, latt2old
    gc.collect()


def make_into_supercell(coor,latt,sublattice):

    if sublattice == 2:
        lattice_direction = [random.randint(0,2)]
    
    if sublattice == 4:
        lattice_direction = [random.randint(0,2), random.randint(0,2)]
    
    for latt_dir in lattice_direction: 
        coor,latt = supercell_double(coor,latt,latt_dir)

    return coor, latt
    gc.collect()

def supercell_double(coor,latt,latt_dir):
    
    for i in range(3):
        latt[latt_dir][i] = latt[latt_dir][i]*2

    totnum = -1
    COOR = np.array([[0.0,0.0,0.0] for i in range(2*2*len(coor))])
    for i in range(len(coor)):
        for j in range(2):
            totnum += 1
            COOR[totnum][0] = coor[i][0] + float(j*latt[latt_dir][0]/2)
            COOR[totnum][1] = coor[i][1] + float(j*latt[latt_dir][1]/2)
            COOR[totnum][2] = coor[i][2] + float(j*latt[latt_dir][2]/2)
          
    return COOR, latt
    gc.collect()
 
def read_coo_file(tot_atom_num):
    with open("coo","r") as f:
        for i in range(5):
            f.readline()
        latt = []
        for i in range(3):
            latt.append(float(f.readline().split()[1]))
        f.readline()
        latt += list(map(float,f.readline().split()[0:3]))
        for i in range(3):
            f.readline()
        coor = []
        for i in range(tot_atom_num):
            coor.append(list(map(float,f.readline().split()[2:5])))

    return coor, latt

'''
def write_structure(coor,latt,atomnumlist,atomnamelist):
    atomtype = []
    for i in range(len(atomnumlist)):
        atomtype = atomtype + [i+1 for j in range(atomnumlist[i])]
    lx = latt[0]
    ly = latt[1]
    lz = latt[2]
    p_xy = latt[3]
    p_xz = latt[4]
    p_yz = latt[5]
    with open("coo","w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write(str(sum(atomnumlist))+" atoms\n")
        fw.write(str(len(atomnumlist))+" atom types\n")
        fw.write("\n")
        fw.write("0.000000 %10.6f xlo xhi\n"%(lx))
        fw.write("0.000000 %10.6f ylo yhi\n"%(ly))
        fw.write("0.000000 %10.6f zlo zhi\n"%(lz))
        fw.write("\n")
        fw.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy, p_xz, p_yz))
        fw.write("\n")
        fw.write("Atoms\n")
        fw.write("\n")
        for i in range(sum(atomnumlist)):
            fw.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, atomtype[i], coor[i][0], coor[i][1], coor[i][2]))
    gc.collect()
'''
