import os
import sys
import shutil
from module_input import *
import random
import subprocess
import numpy as np
import math
from module_potcar import *

input_name  = sys.argv[1]
mpi_command = sys.argv[2]
NPROC       = sys.argv[3]

inp_file, error_message = input_yaml(input_name)

ENCUT = inp_file['DFT_setting']['ENCUT']
NPAR  = inp_file['DFT_setting']['NPAR']
KPAR  = inp_file['DFT_setting']['KPAR']
vasp_gam = inp_file['DFT_setting']['vasp_gam']
vasp_std = inp_file['DFT_setting']['vasp_std']

pot_dir=inp_file['DFT_setting']['potential_directory']


folder= inp_file['output_dir']

k_spacing = inp_file['DFT_setting']['k-spacing']

# getting info from output folder
os.chdir(folder)

f = open("BestResults","r")
breadlines = f.readlines()
f.close()

with open("CONTCARs","r") as f:
    poscars = f.readlines()

atomnamelist = []
atomnumlist  = []
with open("structure","r") as f:
    atomnamelist = list(f.readline().split())
    atomnumlist  = list(map(int,f.readline().split()))

tot_atom_num = sum(atomnumlist)


with open("totalGen","r") as f:
    freadlines = f.readlines()
    num_struct = len(freadlines) - 2

dat = []
for i in range(2,len(breadlines)):
    gen = int(breadlines[i].split()[0])
    pop = int(breadlines[i].split()[1])
    anti= float(breadlines[i].split()[5])
    E   = float(breadlines[i].split()[2])
    dat.append([gen,pop,anti,E])

num_struct  = int(len(poscars)/(8+tot_atom_num))

# make dft folder inside the output folder
os.mkdir("one_shot")
os.chdir("one_shot")

# make incar


# make k-point module
def make_kpt(k_spacing):
    with open("POSCAR","r") as f:
        f.readline()
        f.readline()
        latt = []
        latt.append(list(map(float,f.readline().split()[0:3])))
        latt.append(list(map(float,f.readline().split()[0:3])))
        latt.append(list(map(float,f.readline().split()[0:3])))
        latt = np.array(latt)

        k0 = math.ceil(1.0/(np.linalg.norm(latt[0])*k_spacing))
        k1 = math.ceil(1.0/(np.linalg.norm(latt[1])*k_spacing))
        k2 = math.ceil(1.0/(np.linalg.norm(latt[2])*k_spacing))

    with open("KPOINTS","w") as f:
        f.write("Auto k-point\n")
        f.write(" 0\n")
        f.write("Monk-horst\n") 
        f.write("%d  %d  %d\n"%(k0,k1,k2))
        f.write("0 0 0")

    return [k0,k1,k2]
 
# run VASP
lmp_E = []
gens  = []
pops  = []

# make INCAR ISYM = 1

shutil.copy("../../src/INCAR","./")
with open("INCAR","a") as f:
    f.write("   ENCUT = "+str(ENCUT)+"\n")
    f.write("   NPAR  = "+str(NPAR)+"\n")
    f.write("   KPAR  = "+str(KPAR)+"\n")
    f.write("   ISYM  =  1\n")

# run vasp with ISYM = 1
i_struct = 0
for i in range(num_struct):
    for j in range(len(dat)):
        if "generation: {} / population: {}\n".format(dat[j][0],dat[j][1]) in poscars[i*(8+tot_atom_num)]:
            lmp_E.append(dat[j][3])
            gens.append(dat[j][0])
            pops.append(dat[j][1])

            i_struct += 1


            # make POSCAR and kpt
            with open("POSCAR_"+str(i_struct),'w') as f:
                for k in range((tot_atom_num+8)*i,(tot_atom_num+8)*(i+1)):
                    f.write(poscars[k])

            shutil.copy("POSCAR_"+str(i_struct),"POSCAR")
            kptarray = make_kpt(k_spacing)
            shutil.copy("KPOINTS","KPOINTS"+str(i_struct))


            # make POTCAR 
            if i_struct == 1:
                make_potcar(pot_dir)
            
            # run VASP
            if kptarray == [1,1,1]:
                with open("stdout.x"+str(i_struct),"w") as f:
                    subprocess.call([mpi_command,'-np',NPROC,vasp_gam],stdout=f)

            else:            
                with open("stdout.x"+str(i_struct),"w") as f:
                    subprocess.call([mpi_command,'-np',NPROC,vasp_std],stdout=f)

            os.rename('OUTCAR','OUTCAR'+str(i_struct))

# run vasp with ISYM = 0
os.remove("INCAR")
shutil.copy("../../src/INCAR","./")

with open("INCAR","a") as f:
    f.write("   ENCUT = "+str(ENCUT)+"\n")
    f.write("   NPAR  = "+str(NPAR)+"\n")
    f.write("   KPAR  = "+str(KPAR)+"\n")
    f.write("   ISYM  =  0\n")
    f.write("   SYMPREC  =  1e-10\n")

vasp_E = []
for i in range(1,i_struct+1):
    try:
        check = subprocess.check_output(['grep',"free  ",'OUTCAR'+str(i)])
    except:
        check = ''

    if check =='':

        os.remove("OUTCAR"+str(i))

        shutil.copy("POSCAR_"+str(i),"POSCAR")
        kptarray = make_kpt(k_spacing)
   
        if kptarray == [1,1,1]:
            with open("stdout.x"+str(i),"w") as f:
                subprocess.call([mpi_command,'-np',NPROC,vasp_gam],stdout=f)

        else:            
            with open("stdout.x"+str(i),"w") as f:
                subprocess.call([mpi_command,'-np',NPROC,vasp_std],stdout=f)

        shutil.copy('OUTCAR','OUTCAR'+str(i))
    
        freeE = float(subprocess.check_output(['grep',"free  ",'OUTCAR'+str(i)]).split()[4])

        vasp_E.append(freeE)

    else:
        vasp_E.append(float(check.split()[4]))

# write correlation
os.chdir("../")

with open("correlation","w") as f:
    f.write("gen   pop   vasp_E  lammps_E\n")
    for i in range(i_struct):
        f.write("%d   %d   %.3f   %.3f\n"%(gens[i],pops[i],vasp_E[i],lmp_E[i]))


