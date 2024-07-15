import pyrandspg
import time
import random
import numpy as np

atomnumlist = [16,8,16]
tolerance  = [[[1,1],2.2] , [[1,2],2.00], [[1,3],1.30], [[2,2], 1.70], [[2,3], 1.40], [[3,3],2.20]]

#atomnumlist = [16,8,8,8]
#tolerance  = [[[1,1],2.2] , [[1,2],2.00], [[1,3],1.30], [[1,4],1.5], [[2,2], 1.70], [[2,3], 1.40], [[2,4],1.5], [[3,3],2.20], [[3,4], 1.50], [[4,4], 1.50]]

scale_for_random = 1.0

Vmin = 300.0

####
T1 = time.time()

lmin = (Vmin*0.7)**(1.0/3.0) * 0.5
lmax = (Vmin*1.3)**(1.0/3.0) * 2.0

for i in range(len(tolerance)):
    tolerance[i][1] = tolerance[i][1] * scale_for_random

atom_array = []
for i in range(len(atomnumlist)):
    atom_array += [i+1]*atomnumlist[i]


# generate random_structure
while 1:
    spg = random.randint(1,231)
    input_ = pyrandspg.RandSpgInput(spg, atom_array, pyrandspg.LatticeStruct(lmin, lmin, lmin, 60.0, 60.0, 60.0), pyrandspg.LatticeStruct(lmax, lmax, lmax, 120.0, 120.0, 120.0), 0.7, Vmin*0.7, Vmin*1.3, 100, tolerance, False)

    crystal = pyrandspg.RandSpg.randSpgCrystal(input_)

    structure_= crystal.getPOSCARString()
    structure = structure_.split('\n')
    
    if 'nan' not in structure_:
        break

T2 = time.time()
    
elementorder = list(map(int,structure[5].split()))
elementnumber= list(map(int,structure[6].split()))
atomtypenum  = len(elementnumber)

written_array = []
for i in range(len(elementnumber)):
    written_array += [elementorder[i]]*elementnumber[i]

written_coor = []
tot_atom_num = sum(atomnumlist)
for i in range(tot_atom_num):
    written_coor.append(list(map(float,structure[8+i].split())))

coor_dir = []
for i in range(len(atomnumlist)):
    for j in range(tot_atom_num):
        if written_array[j] == i+1:
            coor_dir.append(written_coor[j])

latt = []
for i in range(2,5):
    latt.append(list(map(float,structure[i].split())))

latt     = np.array(latt)
coor_dir = np.array(coor_dir)

coor = coor_dir.dot(latt)
print (T2-T1)

print ("1.000")
for i in range(3):
    for j in range(3):
        print (latt[i][j],end =' ' )
    print ()
print ("Li Hf N")
print ("16 8 16")
print ("Cartesian")

for i in range(tot_atom_num):
    for j in range(3):
        print (coor[i][j],end =' ' )
    print ()
