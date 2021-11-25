import os
import sys
import shutil
from module_input import *
import random
E_limit = 0.1
 
# input
input_name  = sys.argv[1]
mpi_command = sys.argv[2]
iteration   = sys.argv[3]

inp_file, error_message = input_yaml(input_name)

encut = inp_file['DFT_setting']['ENCUT']
folder= inp_file['output_dir']

num_sample   = int(inp_file['sampling_trainingset']['num_of_samples'])

pot_dir = inp_file['DFT_setting']['potential_directory']
vasp_gam= inp_file['DFT_setting']['vasp_gam']
vasp_std= inp_file['DFT_setting']['vasp_std']
vasp_ncl= inp_file['DFT_setting']['vasp_ncl']
NPAR    = inp_file['DFT_setting']['NPAR']
KPAR    = inp_file['DFT_setting']['KPAR']

# getting info from output folder
os.chdir(folder)

absolute_path = os.getcwd() 

f = open("correlation","r")
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

N_high = 0

dat = []
dat_tot = []
for i in range(1,len(breadlines)):
    gen  = int(breadlines[i].split()[0])
    pop  = int(breadlines[i].split()[1])
    lmp  = float(breadlines[i].split()[2])
    vasp = float(breadlines[i].split()[3])

    if lmp - vasp > float(E_limit*tot_atom_num):
        N_high += 1
    else:
        dat.append([gen,pop,vasp,i,lmp])
        
    dat_tot.append([gen,pop,vasp,i,lmp])

num_struct  = int(len(poscars)/(8+tot_atom_num))

# make dft folder inside the output folder
if os.path.isdir('dft') == False:
    os.mkdir("dft")

os.chdir("dft")
shutil.copytree('../../src/amp2_normal', 'src')
shutil.copytree('../../src/amp2_normal0','src0')

# make config.yaml
with open("config.yaml","w") as f:
  f.write("directory:\n")
  f.write("  submit: ./Submit\n")
  f.write("  output: ./Output\n")
  f.write("  done: ./Done\n")
  f.write("  error: ./ERROR\n")
  f.write("  src_path: ./src\n")
  f.write("  pot_path_gga: "+pot_dir+"\n")
  f.write("  pot_path_lda: "+pot_dir+"\n")
  f.write("\n")

  f.write("program:\n")
  f.write("  vasp_std: "+vasp_std+"\n") 
  f.write("  vasp_gam: "+vasp_gam+"\n") 
  f.write("  vasp_ncl: "+vasp_ncl+"\n") 
  f.write("  mpi_command: "+mpi_command+"\n")
  f.write("\n")

  f.write("vasp_parallel:\n")
  f.write("  NPAR:  "+str(NPAR)+"\n")
  f.write("  KPAR:  "+str(KPAR)+"\n")
  f.write("\n")

  f.write("calculation:\n")
  f.write("  kp_test: 1\n")
  f.write("  encut_test: 0\n")
  f.write("  relaxation: 1\n")
  f.write("  band:  0\n")
  f.write("  density_of_states: 0\n")
  f.write("  hse_oneshot: 0\n")
  f.write("  dielectric: 0\n")
  f.write("  effective_mass: 0\n")
  f.write("  plot: 0\n")
  f.write("  magnetic_ordering: 0\n")
  f.write("  large_off: -1\n")
  f.write("\n")

  f.write("relaxation:\n")
  f.write("  converged_inonic_step: 3\n")
  f.write("  new_kp_percent: 20\n")
  f.write("  incar:\n")
  f.write("    NSW: 10\n")

  f.write("cif2vasp:\n")
  f.write("  pot_name:\n")
  f.write("    GGA:\n")
  f.write("      Cs: Cs_sv_GW\n")
  f.write("      Ca: Ca_sv\n")

with open("config0.yaml","w") as f:
  f.write("directory:\n")
  f.write("  submit: ./Submit\n")
  f.write("  output: ./Output\n")
  f.write("  done: ./Done\n")
  f.write("  error: ./ERROR\n")
  f.write("  src_path: ./src0\n")
  f.write("  pot_path_gga: "+pot_dir+"\n")
  f.write("  pot_path_lda: "+pot_dir+"\n")
  f.write("\n")

  f.write("program:\n")
  f.write("  vasp_std: "+vasp_std+"\n") 
  f.write("  vasp_gam: "+vasp_gam+"\n") 
  f.write("  vasp_ncl: "+vasp_ncl+"\n") 
  f.write("  mpi_command: "+mpi_command+"\n")
  f.write("\n")

  f.write("vasp_parallel:\n")
  f.write("  NPAR:  "+str(NPAR)+"\n")
  f.write("  KPAR:  "+str(KPAR)+"\n")
  f.write("\n")

  f.write("calculation:\n")
  f.write("  kp_test: 1\n")
  f.write("  encut_test: 0\n")
  f.write("  relaxation: 1\n")
  f.write("  band:  0\n")
  f.write("  density_of_states: 0\n")
  f.write("  hse_oneshot: 0\n")
  f.write("  dielectric: 0\n")
  f.write("  effective_mass: 0\n")
  f.write("  plot: 0\n")
  f.write("  magnetic_ordering: 0\n")
  f.write("  large_off: -1\n")
  f.write("\n")

  f.write("relaxation:\n")
  f.write("  converged_inonic_step: 3\n")
  f.write("  new_kp_percent: 20\n")
  f.write("  incar:\n")
  f.write("    NSW: 10\n")

  f.write("cif2vasp:\n")
  f.write("  pot_name:\n")
  f.write("    GGA:\n")
  f.write("      Cs: Cs_sv_GW\n")
  f.write("      Ca: Ca_sv\n")

# make structure in Submit folder
os.mkdir('Submit')
os.chdir('Submit')
i_struct = 0



num_sample = num_sample - N_high

 

if len(dat) < num_sample: 
  for i in range(num_struct):
    for j in range(len(dat)):
        if "generation: {} / population: {}\n".format(dat[j][0],dat[j][1]) in poscars[i*(8+tot_atom_num)]:
            i_struct = i_struct + 1            
            with open("POSCAR_"+str(i_struct),'w') as f:
                for k in range((tot_atom_num+8)*i,(tot_atom_num+8)*(i+1)):
                    f.write(poscars[k])


else:
  dat_sorted = sorted(dat, key = lambda x : x[2])
  for j in range(num_sample):
    for i in range(num_struct):
        if "generation: {} / population: {}\n".format(dat_sorted[j][0],dat_sorted[j][1]) in poscars[i*(8+tot_atom_num)]:
            print ("generation: {} / population: {} ".format(dat_sorted[j][0],dat_sorted[j][1]))
            i_struct = i_struct + 1            
            with open("POSCAR_"+str(i_struct),'w') as f:
                for k in range((tot_atom_num+8)*i,(tot_atom_num+8)*(i+1)):
                    f.write(poscars[k])


fstr = open("../../str_list"+iteration,"w")
if N_high != len(dat_tot):
    fstr.write("[final" + iteration+"-low: 2.0]\n")
for i in range(len(dat_tot)):
    if dat_tot[i][4] - dat_tot[i][2] <= float(E_limit*tot_atom_num):
        fstr.write(absolute_path+"/one_shot"+iteration+"/OUTCAR"+str(dat_tot[i][3])+"  ::1000000\n")

if N_high != 0:
    fstr.write("[final" + iteration+"-high: 0.2]\n")
for i in range(len(dat_tot)):
    if dat_tot[i][4] - dat_tot[i][2] > float(E_limit*tot_atom_num):
        fstr.write(absolute_path+"/one_shot"+iteration+"/OUTCAR"+str(dat_tot[i][3])+"  ::1000000\n")

 
fstr.close()
