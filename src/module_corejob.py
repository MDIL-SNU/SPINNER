import os
import subprocess
import shutil
from module_structure import *
from module_rdf import *
import time
import module_lammps
from module_lammps import *
import math

def core_job(gen, rank, pop, parent, struct_type,Vmin, inp_file, src_path, output_path, atomnamelist, atomnumlist, parent_data, continue_, Emin, tolerance_matrix, average_atom_e, startnum):
   
    rerelax_best      = inp_file['structure']['re-relax_best']
    tot_atom_num = inp_file['tot_atom_num']
    atom_num = len(inp_file['material'])

    accurate_potential = inp_file['relax_condition']['further_calculate_with_accurate_potential']    

    #make folder if it does not exist
    if os.path.isdir(output_path+"/pop"+str(rank)) == False:
        os.mkdir(output_path+"/pop"+str(rank))
        shutil.copy(output_path+"/input.yaml", output_path+"/pop"+str(rank))
        shutil.copy(output_path+"/potential",  output_path+"/pop"+str(rank))
        shutil.copy(output_path+"/potential10",  output_path+"/pop"+str(rank))
        shutil.copy(output_path+"/structure",  output_path+"/pop"+str(rank))
        if accurate_potential == True:
            shutil.copy(output_path+'/potential_accurate',  output_path+"/pop"+str(rank))

    os.chdir(output_path+"/pop"+str(rank))

    # make structure

    t1 = time.time()

    if struct_type == 0:
        shutil.copy("../random_structure/coo_"+str(pop),"coo")
        if gen == startnum:
            coor0 = [[parent_data['coordinate'][i][j] for j in range(3)] for i in range(len(parent_data['coordinate']))]
            latt0 = [parent_data['lattice'][j] for j in range(6)]
        else:
            coor0, latt0 = read_coo_file(tot_atom_num)

    if struct_type == 1:
        coor0, latt0 = con_permutation(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file)
    if struct_type == 2:
        coor0, latt0 = all_permutation(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file) 
    if struct_type == 3:
        coor0, latt0 = latticemutation(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file)
    if struct_type == 4:
        coor0, latt0 = kept_best(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file)
    if struct_type == 5:
        coor0, latt0 = crossover(parent, output_path, parent_data, atomnamelist, atomnumlist, inp_file, average_atom_e, tolerance_matrix)
       
    attempt = 0
 
    t2 = time.time()

    # run calculation
    if rerelax_best == False:
        if struct_type == 4 and gen > continue_+1:
            E = parent_data['E']
            V = parent_data['V']
            latt = parent_data['lattice']
            coor = parent_data['coordinate']
            status = ""
        else:
            E, V, latt, coor, rdf0, atom_e, status = run_lammps(Emin, tolerance_matrix, atomnamelist, atomnumlist, inp_file, accurate_potential,Vmin)

    else:
        E, V, latt, coor, rdf0, atom_e, status = run_lammps(Emin, tolerance_matrix, atomnamelist, atomnumlist, inp_file,accurate_potential,Vmin)

    t3 = time.time()

    # write rdfs
    if rerelax_best == False:
        if struct_type == 4 and gen > continue_+1:
            rdfs = parent_data['rdf']
        else:
            rdfs = read_rdf_from_lammps(atomnumlist, inp_file, rdf0)
    else: 
        rdfs = read_rdf_from_lammps(atomnumlist, inp_file, rdf0)
    t4 = time.time()

    # write POSCAR CONTCAR
    poscar = "generation: "+str(gen)+" / population: "+str(pop)+"\n"
    poscar += "1.0000\n"
    poscar += str(latt0[0]) + " 0.0 0.0\n"
    poscar += str(latt0[3]) + " " + str(latt0[1]) + " 0.0\n"
    poscar += str(latt0[4]) + " " + str(latt0[5]) + " " + str(latt0[2]) + "\n"

    for j in range(len(atomnamelist)):
        poscar += atomnamelist[j] + " "
    poscar += "\n"
    for j in range(len(atomnumlist)):
        poscar += str(atomnumlist[j]) + " "
    poscar += "\n"

    poscar += "Cartesian\n"
    for j in range(tot_atom_num):
        poscar += str(coor0[j][0])+" "
        poscar += str(coor0[j][1])+" "
        poscar += str(coor0[j][2])+" "
        poscar += "\n"

    #write CONTCARs
    contcar = "generation: "+str(gen)+" / population: "+str(pop)+"\n"
    contcar += "1.000\n"
    contcar += str(latt[0]) + " 0.0 0.0\n"
    contcar += str(latt[3]) + " " + str(latt[1]) + " 0.0\n"
    contcar += str(latt[4]) + " " + str(latt[5]) + " " + str(latt[2]) + "\n"

    for j in range(len(atomnamelist)):
        contcar += atomnamelist[j] + " "
    contcar += "\n"
    for j in range(len(atomnumlist)):
        contcar += str(atomnumlist[j]) + " "
    contcar += "\n"

    contcar += "Cartesian\n"
    for j in range(tot_atom_num):
        contcar += str(coor[j][0])+" "
        contcar += str(coor[j][1])+" "
        contcar += str(coor[j][2])+" "
        contcar += "\n"

    t5 = time.time()

    del coor0, latt0, rdf0

    return t2-t1, t3-t2, t4-t3, t5-t4, E, V, poscar, contcar, rdfs, latt, coor, atom_e, status, attempt
