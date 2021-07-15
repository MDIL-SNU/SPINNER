import numpy as np
import bisect
import random
import math
from module_lammps import *
import lammps
from lammps import lammps
import gc

##### convert to lammps function ####
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


###### choose function ######
def cdf(weights):
    total = sum(weights)
    result = []
    cumsum = 0
    for w in weights:
        cumsum += w
        result.append(cumsum / total)
    return result
    gc.collect()

def choice(population, weights):
    assert len(population) == len(weights)
    cdf_vals = cdf(weights)
    x = random.random()
    idx = bisect.bisect(cdf_vals, x)
    return population[idx]

    del cdf_vals
    gc.collect()

##### main codes #####
def cut_the_slab(coor, latt, ae, atomnumlist, average_atom_e, latt_direction, inp_file):

    # make atom type array   
    atom_type_array = []
    atomtype = 0
    for i in atomnumlist:
        atom_type_array += [atomtype]*i
        atomtype += 1

    # make super coordinate
    tot_atom_num = sum(atomnumlist)
    COOR = []
    AE   = []    
    TYPE = []
    for i in range(tot_atom_num):
        
        element_type = atom_type_array[i]

        COOR.append([coor[i][0], coor[i][1], coor[i][2]])
        COOR.append([coor[i][0]+1.0, coor[i][1]+0.0, coor[i][2]+0.0])
        COOR.append([coor[i][0]+0.0, coor[i][1]+1.0, coor[i][2]+0.0])
        COOR.append([coor[i][0]+0.0, coor[i][1]+0.0, coor[i][2]+1.0])
        COOR.append([coor[i][0]+1.0, coor[i][1]+1.0, coor[i][2]+0.0])
        COOR.append([coor[i][0]+1.0, coor[i][1]+0.0, coor[i][2]+1.0])
        COOR.append([coor[i][0]+0.0, coor[i][1]+1.0, coor[i][2]+1.0])
        COOR.append([coor[i][0]+1.0, coor[i][1]+1.0, coor[i][2]+1.0])
        AE.append(ae[i] - average_atom_e[element_type])
        AE.append(ae[i] - average_atom_e[element_type])
        AE.append(ae[i] - average_atom_e[element_type])
        AE.append(ae[i] - average_atom_e[element_type])
        AE.append(ae[i] - average_atom_e[element_type])
        AE.append(ae[i] - average_atom_e[element_type])
        AE.append(ae[i] - average_atom_e[element_type])
        AE.append(ae[i] - average_atom_e[element_type])
        TYPE.append(element_type)
        TYPE.append(element_type)
        TYPE.append(element_type)
        TYPE.append(element_type)
        TYPE.append(element_type)
        TYPE.append(element_type)
        TYPE.append(element_type)
        TYPE.append(element_type)

    # make nx,ny,nz lists
    N = inp_file['crossover_condition']['num_of_grid_for_cut']
    dE= inp_file['crossover_condition']['energy_range_for_cut_select']

    E = []
    weights = []
    n_arrays = []

    for axis in range(3):
        for n in range(N):

            l = float(n/N)
            n_arrays.append([axis,l])

            E_temp = 0.0
            atom_num_in_slab = 0

            for atom in range(8*tot_atom_num):
                in_the_range = 1
                for k in range(3):
                    if k == axis:
                        if l >= COOR[atom][k] or COOR[atom][k] > l+0.5:
                            in_the_range = 0
                    else:
                        if 0.0 > COOR[atom][k] or COOR[atom][k] > 1.0:
                            in_the_range = 0

                if in_the_range == 1:
                    E_temp += AE[atom]
                    atom_num_in_slab += 1

            if (latt_direction != -1 and axis != latt_direction) or atom_num_in_slab == 0:
                E.append(100.0)
                case = 1
                weights.append(0.0000)
            else:
                E.append(E_temp)
                case = 0
                weights.append(float(np.exp(-E_temp/dE/float(atom_num_in_slab))))



    # choose one nx, ny, nz value 
    population = [i for i in range(3*N)]

    if sum(weights) == 0.0:
        with open("/home/kang1717/cryspnn/cryspnn_pyrandspg_odin/error_message","a") as fw:
           fw.write(str(case))
#          for i in range(len(COOR)):
#              fw.write("%f %f %f\n"%(COOR[i][0],COOR[i][1],COOR[i][2]))

    chosen_num = choice(population, weights)

    selected_axis = n_arrays[chosen_num][0]
    l             = n_arrays[chosen_num][1]
    E_selected    = E[chosen_num]     

    selected_coor = []
    selected_type = []

    # select atoms from nx, ny, nz value
    for atom in range(8*tot_atom_num):

        in_the_range = 1
        for k in range(3):
            if k == selected_axis:
                if l >= COOR[atom][k] or COOR[atom][k] > l+0.5:
                   in_the_range = 0
            else: 
                if 0.0 > COOR[atom][k] or COOR[atom][k] > 1.0:
                   in_the_range = 0

        if in_the_range == 1:
            selected_coor.append([COOR[atom][0],COOR[atom][1],COOR[atom][2]])
            selected_type.append(TYPE[atom])

    # rotate lattice
    selected_atom_num = len(selected_coor)
    selected_latt = []
    
    latt_temp = [[0.0 for i in range(3)] for j in range(3)]

    for i in range(3):
        for j in range(3):
            latt_temp[i][j] = latt[i][j]

    if selected_axis != 0:
        for i in range(3):  # replace 0 and selected_axis
            latt_temp[0][i]             = latt[selected_axis][i]
            latt_temp[selected_axis][i] = latt[0][i]

        for i in range(selected_atom_num):
            temp                            = selected_coor[i][selected_axis]
            selected_coor[i][selected_axis] = selected_coor[i][0]
            selected_coor[i][0]             = temp

    else:
        for i in range(3):
            latt_temp.append([latt[i][0],latt[i][1],latt[i][2]])

    latt_temp = np.array(latt_temp)

    # change y z if necessary (always ly > lz)
    ly = np.linalg.norm(latt_temp[1])
    lz = np.linalg.norm(latt_temp[2])

    if ly < lz:
        for i in range(3):
            temp            = latt_temp[1][i]
            latt_temp[1][i] = latt_temp[2][i]
            latt_temp[2][i] = temp

    # convert into lammps format: [[xx,0,0],[xy,yy,0],[xz,yz,zz]]
    p_a= math.sqrt(latt_temp[0][0]**2.0 + latt_temp[0][1]**2.0 + latt_temp[0][2]**2.0);
    p_b= math.sqrt(latt_temp[1][0]**2.0 + latt_temp[1][1]**2.0 + latt_temp[1][2]**2.0);
    p_c= math.sqrt(latt_temp[2][0]**2.0 + latt_temp[2][1]**2.0 + latt_temp[2][2]**2.0);
    alpha=  angle(latt_temp[1],latt_temp[2]);  beta=  angle(latt_temp[0],latt_temp[2]); gamma=  angle(latt_temp[0],latt_temp[1]);  # Angles in degree
    alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians

    lx=   p_a
    p_xy= p_b * math.cos(gammar)
    p_xz= p_c * math.cos(betar)
    ly=   math.sqrt(p_b**2.0 - p_xy**2.0)
    p_yz= (p_b*p_c*math.cos(alphar)-p_xy*p_xz)/(ly)
    lz=   math.sqrt(p_c**2.0 - p_xz**2.0 - p_yz**2.0)

    selected_latt = [[lx,0.0,0.0],[p_xy,ly,0.0],[p_xz,p_yz,lz]] 

    # shift atoms to the center
    for i in range(selected_atom_num):
        selected_coor[i][0] += 0.25 - l

    output = {}
    output['latt'] = selected_latt
    output['coor'] = selected_coor
    output['type'] = selected_type

    return selected_axis, E_selected, output
    del atom_type_array, COOR, AE, TYPE, E, weights, n_array, population, chosen_num
    del selected_axis, selected_coor, selected_type, selected_latt, latt_temp
    gc.collect()

# choose lattice direction that fits best in to the reference lattice
def decide_direction(latt_change, latt_ref,dir0):
    
    latt_change = np.array(latt_change)
    latt_ref    = np.array(latt_ref)

    dir1 = (dir0+1)%3
    dir2 = (dir0+2)%3

    l1_ref      = np.linalg.norm(latt_ref[dir1])
    l2_ref      = np.linalg.norm(latt_ref[dir2])

    if l2_ref > l1_ref:
        temp = l2_ref
        l2_ref = l1_ref
        l1_ref = temp

    A = [[0,0.0],[1,0.0],[2,0.0]]
    A = np.array(A)

    for i in range(3):
        dir1 = (i+1)%3
        dir2 = (i+2)%3
        l1_change = np.linalg.norm(latt_change[dir1])
        l2_change = np.linalg.norm(latt_change[dir2])

        if l2_change > l1_change:
            temp = l2_change
            l2_change = l1_change
            l1_change = temp

        A[i][1] = (l1_change - l1_ref)**2.0 + (l2_change - l2_ref)**2.0

    A_sorted = sorted(A, key = lambda x : x[1])

    if random.random() < 0.6:
        return int(A_sorted[0][0])
    else:
        return random.randint(0,2)

    gc.collect()
    del A_sorted
    del A

def integrate_slabs(slab1, slab2, winner, atomnumlist, atomnamelist, average_atom_e, inp_file):

    # integrate lattice
    lattice = [[0.0 for i in range(3)] for j in range(3)]

    if winner == 1:
         for i in range(3):
             lattice[1][i] = slab1['latt'][1][i]
             lattice[2][i] = slab1['latt'][2][i]

    else:
         for i in range(3):
             lattice[1][i] = slab2['latt'][1][i]
             lattice[2][i] = slab2['latt'][2][i]

    for i in range(3):
        lattice[0][i] = float((slab1['latt'][0][i]+slab2['latt'][0][i])/2.0)

    info = []
    N = inp_file['crossover_condition']['grid_for_shift']

    lmp = lammps()
    create_transient_lmp(lattice, atomnamelist, lmp)

    # main loop
    COORS = []
    temp = 0

    shifty = float(random.random()/N)
    shiftz = float(random.random()/N)
    for xy_mirror in range(2):
        for xz_mirror in range(2):
            for yz_mirror in range(2):
                for y_shift in range(N):
                    for z_shift in range(N):
                        
                        dy = float(y_shift/N)+shifty
                        dz = float(z_shift/N)+shiftz

                        # integrate coordinates
                        COOR, TYPE = make_coor(slab1,slab2,xy_mirror, xz_mirror, yz_mirror, dy, dz, atomnamelist, average_atom_e)

                        E = evaluate_energy(lattice, COOR, TYPE, atomnumlist, atomnamelist, lmp, temp)
                        COORS.append([COOR,E,xy_mirror,xz_mirror,yz_mirror,y_shift,z_shift])
                        temp += 1                 
    COORS_sorted = sorted(COORS, key = lambda x : x[1])
    COOR_selected = COORS_sorted[0][0]

    # fine search
    COORS_fine = []
    for y_shift in range(-round(N/2),round(N/2)):
        for z_shift in range(-round(N/2),round(N/2)):

            xy_mirror = COORS_sorted[0][2]
            xz_mirror = COORS_sorted[0][3]
            yz_mirror = COORS_sorted[0][4]

            dy = float(y_shift/N/N) + COORS_sorted[0][5] + shifty
            dz = float(z_shift/N/N) + COORS_sorted[0][6] + shiftz

            COOR, TYPE = make_coor(slab1,slab2,xy_mirror, xz_mirror, yz_mirror, dy, dz, atomnamelist, average_atom_e)
            E = evaluate_energy(lattice, COOR, TYPE, atomnumlist, atomnamelist, lmp, temp)
            temp += 1
                        
            COORS_fine.append([COOR,E,xy_mirror,xz_mirror,yz_mirror,y_shift,z_shift])

    COORS_sorted = sorted(COORS_fine, key = lambda x : x[1])
    COOR_selected = COORS_sorted[0][0]

    lmp.close()
    del lmp

    # save interface region
    latt1 = np.array(slab1['latt'])
    latt2 = np.array(slab2['latt'])

    x1 = np.linalg.norm(latt1)
    x2 = np.linalg.norm(latt2)

    interface = x2/(x1+x2) 

    # calculate atomic energy
    atom_e = atomic_energy(lattice, COOR, TYPE, average_atom_e,atomnamelist,atomnumlist)
    
    integrated_structure = {}
    integrated_structure['coor']   = COOR
    integrated_structure['type']   = TYPE
    integrated_structure['atom_e'] = atom_e
    integrated_structure['latt']   = lattice

    return integrated_structure, interface
    gc.collect()
    del lattice, info, COORS, COOR_selected, COORS_fine, COORS_sorted, COOR, TYPE    
  
def create_transient_lmp(latt,atomnamelist,lmp):
    
    lx   = latt[0][0]
    ly   = latt[1][1]
    lz   = latt[2][2]
    p_xy = latt[1][0]
    p_xz = latt[2][0]
    p_yz = latt[2][1]
    
    with open("coo_temp","w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write("1 atoms\n")
        fw.write(str(len(atomnamelist))+" atom types\n")
        fw.write("\n")
        fw.write("0.000000 %10.6f xlo xhi\n"%(lx))
        fw.write("0.000000 %10.6f ylo yhi\n"%(ly))
        fw.write("0.000000 %10.6f zlo zhi\n"%(lz))
        fw.write("\n")
        fw.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy, p_xz, p_yz))
        fw.write("\n")
        fw.write("Atoms\n")
        fw.write("\n1 1 0 0 0")
    
    # start lammps
    lmp.command("units           metal")
    lmp.command("newton          on")
    lmp.command("dimension       3")
    lmp.command("boundary        p p p")
    lmp.command("atom_style 	atomic ")
    lmp.command("atom_modify     map yes")
    lmp.command("box tilt large")
    lmp.command("read_data coo_temp  ")
    lmp.command("pair_style nn   ")
    sentence = "pair_coeff * *  potential "
    for atom in atomnamelist:
        sentence +=  atom + " "
    lmp.command(sentence)
    for i in range(len(atomnamelist)):
        sentence = "mass {} 15.0".format(i+1)
        lmp.command(sentence)
    lmp.command("neighbor 2 bin")
    lmp.command("neigh_modify 	every 1 delay 0 check yes")
    lmp.command("compute 	_rg all gyration ")
    lmp.command("group one id 1")
    lmp.command("delete_atoms group one")
    gc.collect()

def atomic_energy(lattice, COOR, TYPE, average_atom_e,atomnamelist, atomnumlist):
    tot_atom_num = len(TYPE)
    lmp = lammps()
    # make coo
    
    lx   = lattice[0][0]
    ly   = lattice[1][1]
    lz   = lattice[2][2]
    p_xy = lattice[1][0]
    p_xz = lattice[2][0]
    p_yz = lattice[2][1]
    
    with open("coo_temp","w") as fw:
        fw.write(" POSCAR to lmp\n")
        fw.write("\n")
        fw.write(str(len(TYPE))+" atoms\n")
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
        for i in range(len(TYPE)):
            fw.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, TYPE[i]+1, COOR[i][0], COOR[i][1], COOR[i][2]))
    
    # start lammps
    lmp.command("units           metal")
    lmp.command("newton          on")
    lmp.command("dimension       3")
    lmp.command("boundary        p p p")
    lmp.command("atom_style 	atomic ")
    lmp.command("atom_modify     map yes")
    lmp.command("box tilt large")
    lmp.command("read_data coo_temp  ")
    lmp.command("pair_style nn   ")
    sentence = "pair_coeff * *  potential "
    for atom in atomnamelist:
        sentence +=  atom + " "
    lmp.command(sentence)
    for i in range(len(atomnamelist)):
        sentence = "mass {} 15.0".format(i+1)
        lmp.command(sentence)
    lmp.command("neighbor 2 bin")
    lmp.command("neigh_modify 	every 1 delay 0 check yes")
    lmp.command("compute 	_rg all gyration ")

    # atomic energy
    lmp.command("compute atompe all pe/atom")
    lmp.command("minimize 0 0 0 0")
    atom_e = lmp.extract_compute("atompe",1,1) 
    atom_e_0 = atom_e[0:len(TYPE)]
    lmp.close()
    del lmp
    atom_E = [0.0 for i in range(len(TYPE))]
 
    for i in range(len(TYPE)):
        atom_E[i] = atom_e_0[i] - average_atom_e[TYPE[i]]

    return atom_E
    gc.collect()

def make_coor(slab1,slab2, xy_mirror, xz_mirror, yz_mirror, y_shift, z_shift, atomnumlist, average_atom_e):

    coor1 = [[0.0, 0.0, 0.0] for i in range(len(slab1['coor']))]
    coor2 = [[0.0, 0.0, 0.0] for i in range(len(slab2['coor']))]

    N1 = len(coor1)
    for i in range(N1):
        coor1[i][0] = slab1['coor'][i][0]
        coor1[i][1] = slab1['coor'][i][1]
        coor1[i][2] = slab1['coor'][i][2]

    N2 = len(coor2)
    for i in range(N2):
        coor2[i][0] = slab2['coor'][i][0]
        coor2[i][1] = slab2['coor'][i][1]
        coor2[i][2] = slab2['coor'][i][2]
   
    latt1 = np.array(slab1['latt'])
    latt2 = np.array(slab2['latt'])

    x1 = np.linalg.norm(latt1)
    x2 = np.linalg.norm(latt2)

    interface = x2/(x1+x2) 

    for i in range(N1):
        if xz_mirror == 1:
            coor1[i][1] = 1.0 - coor1[i][1]     
        if xy_mirror == 1:
            coor1[i][2] = 1.0 - coor1[i][2]     
        if yz_mirror == 1:
            coor1[i][0] = 1.0 - coor1[i][0]

        coor1[i][0] = (coor1[i][0]-0.5)*2*x1/(x1+x2) + 0.5 + (0.5*x1 + x2)/(x1+x2) - 0.5
        coor1[i][1] = coor1[i][1] + y_shift
        coor1[i][2] = coor1[i][2] + z_shift
    
        
    for i in range(N2):
        coor2[i][0] = (coor2[i][0]-0.5)*2*x2/(x1+x2) + 0.5 + (0.5*x2)/(x1+x2) - 0.5

    for j in range(3):
        for i in range(N1):
            if coor1[i][j] > 1.0:
                coor1[i][j] -= 1.0

    COOR = []
    TYPE = []

    for i in range(len(atomnumlist)):
        for j in range(N1):
            if slab1['type'][j] == i:
                COOR.append([coor1[j][0],coor1[j][1],coor1[j][2]])
                TYPE.append(i)
        for j in range(N2):
            if slab2['type'][j] == i:
                COOR.append([coor2[j][0],coor2[j][1],coor2[j][2]])
                TYPE.append(i)

    return COOR, TYPE
    gc.collect()
    del coor1, coor2, latt1, latt2 

def evaluate_energy(lattice, COOR, TYPE, atomnumlist, atomnamelist, lmp, t):
    
    # change into the cartesian_coordinate
    COOR = np.array(COOR)
    lattice = np.array(lattice)
 
    coor_car = COOR.dot(lattice)
    
    # create atoms
    for i in range(len(COOR)):
        lmp.command("create_atoms {} single {} {} {}".format(TYPE[i]+1,coor_car[i][0],coor_car[i][1],coor_car[i][2]))        
    

    # evaluate energy
    lmp.command("run 1")
    lmp.command("variable e equal etotal")
    lmp.command("variable v equal vol")
    e = lmp.extract_variable("e",0,0)

    # delete all atoms
    lmp.command("delete_atoms group all")

    return e
    gc.collect()
    del coor_car

def remove_redundant_atoms(integrated_structure, atomnamelist, atomnumlist, average_atom_e, inp_file):

    COOR    = integrated_structure['coor']
    TYPE    = integrated_structure['type']
    lattice = integrated_structure['latt']

    tot_atom_num = len(COOR)

    atomnumlist_present = [0 for i in range(len(atomnumlist))]
 
    for i in range(len(atomnumlist)):
        for j in range(len(TYPE)):
            if i == TYPE[j]:
                atomnumlist_present[i] += 1
    
    atomnumlist_difference = [0 for i in range(len(atomnumlist))]

    
    # calculate atomic E
    atom_e = atomic_energy(lattice, COOR, TYPE, average_atom_e,atomnamelist,atomnumlist)

    # eliminate atoms
    vacant_sites = [] # save vacant sites

    COOR_NEW = []
    TYPE_NEW = []

    for i in range(len(atomnumlist)):
        atomnumlist_difference[i] = atomnumlist_present[i] - atomnumlist[i]
        
        info = []

        if atomnumlist_difference[i] > 0: # if element is redundant
            for j in range(tot_atom_num):
                if TYPE[j] == i:
                    info.append([[COOR[j][0],COOR[j][1],COOR[j][2]], atom_e[j]])

            info_sorted = sorted(info, key = lambda x : x[1])

            for j in range(atomnumlist[i]):
                COOR_NEW.append([info_sorted[j][0][0], info_sorted[j][0][1], info_sorted[j][0][2]])
                TYPE_NEW.append(i)
           
            for j in range(atomnumlist[i],atomnumlist_present[i]):
                vacant_sites.append([info_sorted[j][0][0], info_sorted[j][0][1], info_sorted[j][0][2]])

            atomnumlist_difference[i] = 0

        else:
            for j in range(tot_atom_num):
                if TYPE[j] == i:
                    TYPE_NEW.append(i)
                    COOR_NEW.append([COOR[j][0],COOR[j][1],COOR[j][2]])

    return COOR_NEW, TYPE_NEW, atomnumlist_difference, vacant_sites
    gc.collect()
    del COOR, TYPE, lattice, atomnumlist_present, atom_e, info

def add_atoms(coor, latt, atomnumlist_difference, atomnumlist, atomnamelist, tolerance_matrix, inp_file, interface):
    
    coor = np.array(coor)
    latt = np.array(latt)

    coor = coor.dot(latt)  

    shuffle_list = [i for i in range(len(atomnumlist_difference))]

    N = inp_file['crossover_condition']['iteration_for_add_atoms']

    lmp = lammps()
    create_transient_lmp(latt, atomnamelist, lmp)
    
    info = []
    for n in range(N):

        random.shuffle(shuffle_list)
        
        COOR = initialize_coor(coor,atomnumlist,atomnumlist_difference)

        for i in shuffle_list:
    
            if atomnumlist_difference[i] < 0:

                add_num = -atomnumlist_difference[i]
                success_num = 0

                atomtype = i

                for j in range(add_num):
                        COOR = add_one_atom(COOR, latt, atomtype, tolerance_matrix, interface)

        COOR_l, TYPE = convert_to_std_format(COOR)

        E = evaluate_energy(latt, COOR_l, TYPE, atomnumlist, atomnamelist, lmp, n)

        info.append([COOR_l, E])

    info_sorted = sorted(info, key = lambda x : x[1])
    
    COOR_return = info_sorted[0][0]

    lmp.close()
    del lmp
    return COOR_return

    gc.collect()
    del shuffle_list, info, COOR_l, TYPE, COOR, info_sorted

def convert_to_std_format(COOR):

    coor_return = []
    type_return = []

    for i in range(len(COOR)):
        for j in range(len(COOR[i])):
            coor_return.append([COOR[i][j][0],COOR[i][j][1],COOR[i][j][2]])
            type_return.append(i)

    return coor_return, type_return
    gc.collect()


def add_one_atom(COORS, latt, atomtype, tolerance_matrix, interface):

    while 1:
    
        coor_new = random_put(latt, interface)

        if tolerance_satisfy(coor_new, COORS, atomtype, tolerance_matrix) == 1:
            break

    COORS[atomtype].append([coor_new[0],coor_new[1],coor_new[2]])
   
    return COORS
    gc.collect()

def tolerance_satisfy(coor_new, COORS, atomtype, tolerance_matrix):

    satisfy = 1

    for i in range(len(COORS)):

        for j in range(len(COORS[i])):

            for x in range(-1,2):
                for y in range(-1,2):
                    for z in range(-1,2):

                        coor0 = np.array(COORS[i][j])
                        
                        d = np.linalg.norm(coor_new-coor0)

                        if d < tolerance_matrix[atomtype][i]:

                            satisfy = 0

    return satisfy
    gc.collect()

def random_put(latt, interface):

    if random.random() < 0.8:  # put atom in the interface

        if random.random() < 0.5: # 0.45 ~ 0.55

            a = random.random()*0.1+interface-0.05

        else:                     # 0.95 ~ 0.05

            if random.random() < 0.5:

                a = random.random()*0.05+0.95

            else:
                a = random.random()*0.05

    else:
        a = random.random()

    b = random.random()
    c = random.random()

    coor_new = np.array([a,b,c])
    latt     = np.array(latt)

    coor_new = coor_new.dot(latt)

    return coor_new
    gc.collect()

def initialize_coor(coor_ref, atomnumlist, atomnumlist_difference):

    coor = []

    global_atom_num = 0

    for i in range(len(atomnumlist_difference)):
            
        coor.append([])

        atomnum = atomnumlist[i]
        diff    = atomnumlist_difference[i]

        for j in range(atomnum + diff):
           
            coor[i].append([coor_ref[global_atom_num][0],coor_ref[global_atom_num][1],coor_ref[global_atom_num][2]])

            global_atom_num += 1

    return coor
    gc.collect()

def spring(latt,coor,tolerance_matrix,atomnumlist,atomnamelist,atomtype,tot_atom_num):

    lx =   latt[0]
    ly =   latt[1]
    lz =   latt[2]
    p_xy = latt[3]
    p_xz = latt[4]
    p_yz = latt[5]

    with open("coo0","w") as fw:
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
    
    lmp = lammps()
    lmp.command("units           metal")
    lmp.command("newton          on")
    lmp.command("dimension       3")
    lmp.command("boundary        p p p")
    lmp.command("atom_style 	atomic ")
    lmp.command("atom_modify     map yes")
    lmp.command("box tilt large")
    lmp.command("read_data coo0 ")
    lmp.command("pair_style zero 6.0 ")
    lmp.command("pair_coeff * * ")
    for i in range(len(atomnamelist)):
        sentence = "mass {} 15.0".format(i+1)
        lmp.command(sentence)
    lmp.command("neighbor 2 bin")
    lmp.command("neigh_modify 	every 1 delay 0 check yes")
    lmp.command("compute 	_rg all gyration ")               # All atoms will be involved in the simulation http://lammps.sandia.gov/doc/compute.html 

    fix_holdem0(atomnumlist, tolerance_matrix, lmp)
    relax_with_unfixed_lattice(tot_atom_num,lmp,4)
    latt, coor = get_latt_coor(tot_atom_num,lmp)
    lattr = [latt[0][0],latt[1][1],latt[2][2],latt[1][0],latt[2][0],latt[2][1]]
    lmp.close()
    del lmp
    return lattr, coor
    gc.collect()

def fix_holdem0(atomnumlist, tolerance_matrix,lmp):
    atomnum_array = []
    for i in range(len(atomnumlist)):
        atomnum_array += [i]*atomnumlist[i]

    number = 0
    for i in range(len(atomnum_array)):
        for j in range(i,len(atomnum_array)):
            number += 1
            sentence = "fix holdem"+ str(number) + " all restrain lbound " + str(i+1) + "  " + str(j+1) + " 0.1 0.1 " + str(tolerance_matrix[int(atomnum_array[i])][int(atomnum_array[j])])
            lmp.command(sentence) 
            sentence = "fix_modify holdem"+str(number)+" energy yes"
            lmp.command(sentence) 
    gc.collect()
    del atomnum_array
