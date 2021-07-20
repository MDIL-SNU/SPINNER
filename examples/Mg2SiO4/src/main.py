###########################################
### Date: 2021-07-09                    ###
### Sungwoo Kang                        ###
###########################################
import os, sys, subprocess, yaml, shutil, time, shutil
from subprocess import check_output
import random
from time import sleep
from module_input  import *
#from module_log import * 
from module_pre import *
from module_post import *
from module_continue import *
from module_time  import *
from module_clean import *
from module_corejob import *
from module_structure import *

import mpi4py
from mpi4py import MPI
import gc


input_file = str(sys.argv[1])

# Getting calculation information
inp_file, error_message = input_yaml("../"+input_file)

if error_message != "":
    with open("../ERROR_"+input_file,"w") as fw:
        fw.write(error_message)
    sys.exit()

# Append input directory
src_path = os.getcwd()
input_dir = inp_file['input_dir']

# check if there is any missing input files


# MPI setting
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
corenum = comm.Get_size()

num_of_atom = len(inp_file['material'])

# info that need in the main code
generation        = inp_file['structure']['generation']
population        = inp_file['structure']['population']
i_population      = inp_file['structure']['i_population']
num_of_best_input = inp_file['structure']['num_of_best']
rerelax_best      = inp_file['structure']['re-relax_best']
tot_atom_num      = inp_file['tot_atom_num']
#AutoFrac          = inp_file['operator']['Auto_Frac']
output_dir        = inp_file['output_dir']
continue_         = inp_file['continue']['continue_num']
dis_type          = inp_file['similarity_metric']['type']
accurate_potential = inp_file['relax_condition']['further_calculate_with_accurate_potential']
 
output_path = src_path+'/../'+output_dir

rdf_grid = inp_file['similarity_metric']['rdf_grid']

# check if the same-name directory is in the current folder

if rank == 0:
  if os.path.isdir(output_path):
    with open("../ERROR_"+input_file,"w") as fw:
        fw.write("output_dir already exists in the current folder. rename it")
    for i in range(1,corenum):
        comm.send(0, dest = i, tag = 0)
    sys.exit()
  else:
    for i in range(1,corenum):
        comm.send(1, dest = i, tag = 0)
else:
  if_there_is_same_directory = comm.recv(source = 0, tag=0)
  if if_there_is_same_directory == 0:
    sys.exit()

#main core
if rank == 0:    
    #making folders    
    os.mkdir(output_path)
    shutil.copy(src_path+'/../'+input_file,              output_path)
    with open(output_path+"/input.yaml","w") as f:
        yaml.dump(inp_file,f)
    shutil.copy(src_path+'/../'+input_dir+'/potential',  output_path)
    if accurate_potential == True:
        shutil.copy(src_path+'/../'+input_dir+'/potential_accurate',  output_path)
    os.mkdir(output_path+"/random_structure")
    with open(output_path+"/random_structure_log","w") as f:
        f.write(" gen  pop      time    Attempt Spacegroup\n")    

    make_potcar_cutoff_10(output_path)   
 
    # write structure info
    fstructure = open(output_path+'/structure','w')
    f_input = open(output_path+"/input.yaml","r")
    f_input_readlines = f_input.readlines()
    f_input.close()
    for i in range(len(f_input_readlines)):
        if 'material:' in f_input_readlines[i]:
            material_num = i
            break;
    for i in range(material_num+1,material_num+len(inp_file['material'])+1):
        fstructure.write("  ")
        fstructure.write(f_input_readlines[i].split()[0][0:-1])
    fstructure.write("\n")
    for i in range(material_num+1,material_num+len(inp_file['material'])+1):
        fstructure.write("  ")
        fstructure.write(f_input_readlines[i].split()[1])
    fstructure.close()


    #write start log
    with open(output_path+'/totalGen','w') as f_individual:
        f_individual.write("     Generation  Foldernum    Energy (eV)    Volume (A3)          Type          Parents    keptBest(generation,population)     \n")
        f_individual.write("-------------------------------------------------------------------------------------------------------------------------------\n")
    with open(output_path+'/totalbest','w') as f_best:
        f_best.write("     Generation  Foldernum    Energy (eV)    Volume (A3)          Type          Parents    keptBest(generation,population)     \n")
        f_best.write("-------------------------------------------------------------------------------------------------------------------------------\n")
        

    # pre-processing if continue or not
    if continue_ == 0:
        results = {}    # results for previous generation
        results['E']  = [0. for i in range(i_population)]
        results['V']  = [0. for i in range(i_population)]
        results['coordinate'] = [[[0.0,0.0,0.0] for i in range(tot_atom_num)] for k in range(i_population)]
        results['lattice']    = [[0.0 for i in range(6)] for k in range(i_population)]
        results['rdf']        = [[0.0 for i in range(int(rdf_grid*num_of_atom*(num_of_atom+1)/2))] for k in range(i_population)]
        results['atom_e']     = [[0.0 for i in range(tot_atom_num)] for k in range(i_population)]
        oldresults = {}    # results for previous generation
        oldresults['E']  = [0. for i in range(i_population)]
        oldresults['V']  = [0. for i in range(i_population)]
        oldresults['coordinate'] = [[[0.0,0.0,0.0] for i in range(tot_atom_num)] for k in range(i_population)]
        oldresults['lattice']    = [[0.0 for i in range(6)] for k in range(i_population)]
        oldresults['rdf']        = [[0.0 for i in range(int(rdf_grid*num_of_atom*(num_of_atom+1)/2))] for k in range(i_population)]
        oldresults['atom_e']     = [[0.0 for i in range(tot_atom_num)] for k in range(i_population)]

#    results['atomic_energy'] = [0. for i in range(tot_atom_num*i_population)]
    #results['keptbest'] = [0 for i in range(i_population)]
    #making folder & energy V log
   
        history = {} 

        if num_of_best_input == 0:
            num_of_best = 0
        else:
            num_of_best = num_of_best_input

        startnum = 1

        Emin = 100.0

        global_best_results = []

        # write log
        with open(output_path+"/timelog","w") as ftimelog:
            ftimelog.write("---------------------\n")
            ftimelog.write("        START        \n")
            ftimelog.write("---------------------\n")
            ftimelog.write("\n")
        with open(output_path+"/specific_time","a") as f:
            f.write("                             Struct Lammps  RDF  Gets Commun Core#\n")
    else:
        original_dir         = inp_file['continue']['original_dir']
        original_path        = src_path+"/../"+original_dir

        num_of_atom = len(inp_file['material'])
        
#        average_atom_e = update_average_atom_from_file(continue_, original_path, output_path)

        average_atom_e = [0.0 for i in range(num_of_atom)]

        oldresults, pop_info, linelist = getting_info_from_files(continue_, tot_atom_num, original_path, src_path, num_of_atom)
        history = getting_history(src_path,original_path,continue_)
        
        global_best_results = []
        global_best_results, Emin   = update_best_structures(global_best_results,oldresults,inp_file, continue_, pop_info)

        rewrite_contcars_poscars(output_path, original_path, continue_, tot_atom_num)

        pop_num = len(oldresults['E'])

        num_of_best = getting_num_of_best(oldresults, Emin, tot_atom_num, inp_file)

        startnum = continue_ + 1
        

        ANTISEED = make_antiseed_continue(oldresults, global_best_results, inp_file, Emin)

        # write log
        with open(src_path+"/../"+output_dir+"/timelog","w") as ftimelog:
            ftimelog.write("-----------------------------------------------------------------------------\n")
            ftimelog.write('Start from folder '+str(original_dir)+' and generation '+str(continue_)+'\n')
            ftimelog.write("-----------------------------------------------------------------------------\n")
            ftimelog.write("\n")

        with open(output_path+"/specific_time","a") as f:
            f.write("                             Struct Lammps  RDF  Gets Commun Core#\n")

 
        write_logs(continue_, output_path, original_path)
        results={}
        results['E'] = [0.0 for i in range(pop_num)]
        results['V'] = [0.0 for i in range(pop_num)]
        results['coordinate'] = [[[0.0,0.0,0.0] for i in range(tot_atom_num)] for k in range(pop_num)]
        results['lattice']    = [[0.0 for i in range(6)] for k in range(pop_num)]
        results['rdf']        = [[0.0 for i in range(int(rdf_grid*num_of_atom*(num_of_atom+1)/2))] for k in range(pop_num)]
        results['atom_e']     = [[0.0 for i in range(tot_atom_num)] for k in range(pop_num)]

# etc preprocessing
if continue_ == 0:
    startnum = 1
else:
    startnum = continue_+1

if rank == 0:
    with open(output_path+"/structure","r") as fstr:
        atomnamelist = fstr.readline().split()
        atomnumlist  = list(map(int,fstr.readline().split()))

    tolerance_matrix = make_tolerance_matrix(atomnamelist,inp_file)

    Send = [atomnamelist, atomnumlist, tolerance_matrix]
        
    for i in range(1,corenum):
        comm.send(Send, dest = i, tag = 10)
else:
    [atomnamelist, atomnumlist, tolerance_matrix] = comm.recv(source = 0, tag=10)

if rank ==0 and continue_ == 0:
    average_atom_e = [0.0 for i in range(len(atomnumlist))]
    
# read structure info
with open(output_path+'/structure','r') as f:
        atom_names = f.readline().split()
        atom_nums_  = f.readline().split()
        atom_nums = []
        for atom in atom_nums_:
            atom_nums.append(int(atom))

# main loop
for gen in range(startnum,generation+1):

    if rank == 0:

    # preprocessing

        if gen == 1:
            ANTISEED = []

        t1 = time.time()
        with open(output_path+"/timelog","a") as ftimelog:
            ftimelog.write('Generation: %d / '%(gen))

        if gen == 1:
            pop_num = i_population
        else:
            pop_num = population + num_of_best

        with open(src_path+"/../"+output_dir+"/timelog","a") as ftimelog:
            ftimelog.write('Population: %d\n'%(pop_num))
            ftimelog.write('   Preprocessing + Local Optimization\n')

        Vmin, pop_info,random_num= get_pop_info(inp_file, gen, pop_num, oldresults, src_path, Emin, num_of_best, output_path, ANTISEED)    # pop_info: heredity, permutation, ancestor



        #main loop

        if gen == startnum:
          T_struct_random1 = 0.0
          initial_randoms = {}
          tr1 = time.time()
          for pop in range(1,pop_num+1):
            struct_type, parent = get_parent(pop_info,pop)
            if struct_type == 0:
              os.chdir(output_path+"/random_structure")
              coortemp, latttemp = random_gen(Vmin, inp_file, atomnumlist, atomnamelist, tolerance_matrix, str(pop), gen, pop)
              initial_randoms[pop] = [coortemp,latttemp]
              tr2 = time.time()
              T_struct_random1 += tr2-tr1
        
        # random structure generation in the last core
        Send = [gen, 0, 0, random_num, 0, 0, Vmin, 0, 0, 0, 0, 0]
        comm.send(Send, dest=corenum-1, tag=1)
        random_structure_done = 0
 

        pop = 1
        process = 1
        while process < corenum-1 and pop < pop_num + 1:

            struct_type, parent = get_parent(pop_info,pop)
            parent_data = {}
            t_pre = time.time()
              
            if gen == startnum and struct_type == 0:
                parent_data['coordinate'] = initial_randoms[pop][0]
                parent_data['lattice']    = initial_randoms[pop][1]

            else:
                if type(parent) == int:
                    if parent != 0:
                        parent_data['lattice']    = oldresults['lattice'][parent-1]
                        parent_data['coordinate'] = oldresults['coordinate'][parent-1]
                        parent_data['rdf']        = oldresults['rdf'][parent-1]
                        parent_data['E']          = oldresults['E'][parent-1]
                        parent_data['V']          = oldresults['V'][parent-1]
                        parent_data['atom_e']     = oldresults['atom_e'][parent-1]
                else:
                        parent_data['lattice']    = [oldresults['lattice'][parent[0]-1],   oldresults['lattice'][parent[1]-1]]
                        parent_data['coordinate'] = [oldresults['coordinate'][parent[0]-1],oldresults['coordinate'][parent[1]-1]]
                        parent_data['rdf']        = [oldresults['rdf'][parent[0]-1],       oldresults['rdf'][parent[1]-1]]
                        parent_data['E']          = [oldresults['E'][parent[0]-1],         oldresults['E'][parent[1]-1]]
                        parent_data['V']          = [oldresults['V'][parent[0]-1],         oldresults['V'][parent[1]-1]]
                        parent_data['atom_e']     = [oldresults['atom_e'][parent[0]-1],    oldresults['atom_e'][parent[1]-1]]

            t_start = time.time()
            Send                = [gen,Emin, t_start, pop, parent, struct_type, Vmin, parent_data,t_start-t_pre,average_atom_e, startnum, 1]
            comm.send(Send, dest=process, tag=1)
            pop += 1
            process += 1

   
        T_struct = 0.0
        T_lmp    = 0.0
        T_rdf    = 0.0
        T_get    = 0.0
        T_commun = 0.0
        Tmax     = [] 
        
        if gen != startnum:
            del results
            results = {}
            del poscars
            del contcars
        results['E'] = [0.0 for i in range(pop_num)]
        results['V'] = [0.0 for i in range(pop_num)]
        results['coordinate'] = [[[0.0,0.0,0.0] for i in range(tot_atom_num)] for k in range(pop_num)]
        results['lattice']    = [[0.0 for i in range(6)] for k in range(pop_num)]
        results['rdf']        = [[0.0 for i in range(int(rdf_grid*num_of_atom*(num_of_atom+1)/2))] for k in range(pop_num)]
        results['atom_e']     = [[0.0 for i in range(tot_atom_num)] for k in range(pop_num)]
        poscars  = ['' for i in range(pop_num)] 
        contcars = ['' for i in range(pop_num)] 

        received_processes = 1

        while received_processes < pop_num+1:
            [t_start, process, rpop, t_struct, t_lmp, t_rdf, t_get, rE, rV, poscar, contcar, rdfs, latt, coor, atom_e, status, job] = comm.recv(source=MPI.ANY_SOURCE, tag=2)
#            with open(output_path+"/received_process","a") as f:
#                f.write(str(process)+"\n")
 
            if job==2:
                random_structure_done = 1
                T_struct_random = t_start
                del t_start,rpop, t_struct, t_lmp, t_rdf, t_get, rE, rV, poscar, contcar, rdfs, latt, coor, atom_e, status
            else:
                t_end = time.time()
                t_tot = t_struct + t_lmp + t_rdf + t_get 

                T_lmp    += t_lmp + t_struct
                T_rdf    += t_rdf
                T_get    += t_get

                t_commun  = t_end - t_start - t_lmp - t_rdf - t_get - t_struct
                T_commun += t_commun

                with open(output_path+"/specific_time","a") as f:
                    f.write("generation %d / population %d : %.2f   %.2f  %.2f   %.2f   %.2f  %d   %s\n"%(gen,rpop,t_struct,t_lmp,t_rdf,t_get,t_commun,process,status)) 

                poscars[rpop-1]  = poscar
                contcars[rpop-1] = contcar

                Tmax.append(t_tot)
                results['E'][rpop-1]   = rE
                results['V'][rpop-1]   = rV
                results['coordinate'][rpop-1] = coor
                results['lattice'][rpop-1]    = latt
                results['rdf'][rpop-1]        = rdfs
                results['atom_e'][rpop-1]     = atom_e
                received_processes += 1
        
                del t_start,rpop, t_struct, t_lmp, t_rdf, t_get, rE, rV, poscar, contcar, rdfs, latt, coor, atom_e, status

            if pop < pop_num+1:
                struct_type, parent = get_parent(pop_info,pop)
                t_pre = time.time()
                parent_data = {}
            
                if gen == startnum and struct_type == 0:
                  parent_data['coordinate'] = initial_randoms[pop][0]
                  parent_data['lattice']    = initial_randoms[pop][1]

                else:
                  if type(parent) == int:
                      if parent != 0:
                          parent_data['lattice']    = oldresults['lattice'][parent-1]
                          parent_data['coordinate'] = oldresults['coordinate'][parent-1]
                          parent_data['rdf']        = oldresults['rdf'][parent-1]
                          parent_data['E']          = oldresults['E'][parent-1]
                          parent_data['V']          = oldresults['V'][parent-1]
                          parent_data['atom_e']     = oldresults['atom_e'][parent-1]
                  else:
                    parent_data['lattice']    = [oldresults['lattice'][parent[0]-1],   oldresults['lattice'][parent[1]-1]]
                    parent_data['coordinate'] = [oldresults['coordinate'][parent[0]-1],oldresults['coordinate'][parent[1]-1]]
                    parent_data['rdf']        = [oldresults['rdf'][parent[0]-1],       oldresults['rdf'][parent[1]-1]]
                    parent_data['E']          = [oldresults['E'][parent[0]-1],         oldresults['E'][parent[1]-1]]
                    parent_data['V']          = [oldresults['V'][parent[0]-1],         oldresults['V'][parent[1]-1]]
                    parent_data['atom_e']     = [oldresults['atom_e'][parent[0]-1],    oldresults['atom_e'][parent[1]-1]]

                t_start = time.time()
                Send                      = [gen, Emin, t_start, pop, parent, struct_type, Vmin, parent_data,  t_start-t_pre, average_atom_e, startnum, 1]
                comm.send(Send, dest=process, tag=1)
                pop += 1
                del parent_data, Send

        if random_structure_done == 0: 
            [T_struct_random, process, rpop, t_struct, t_lmp, t_rdf, t_get, rE, rV, poscar, contcar, rdfs, latt, coor, atom_e, status, job] = comm.recv(source=corenum-1, tag=2)

        t2 = time.time()
        with open(output_path+"/timelog","a") as ftimelog:
            if rerelax_best == False:
                ftimelog.write('    - TOTAL random struct_generation: %.2f s\n'%(float(T_struct_random)))
                ftimelog.write('    - average lammps time:            %.2f s\n'%(float(T_lmp/population)))
                ftimelog.write('    - average RDF time:               %.2f s\n'%(float(T_rdf/population)))
                ftimelog.write('    - getting info time:              %.2f s\n'%(float(T_get/population)))
                ftimelog.write('    - longest population              %.2f s\n'%(float(max(Tmax))))
                ftimelog.write('    - total time:                     %.2f s\n\n'%(t2-t1))
            if rerelax_best:
                ftimelog.write('    - TOTAL random struct_generation: %.2f s\n'%(float(T_struct_random)))
                ftimelog.write('    - average crossover+lammps time:  %.2f s\n'%(float(T_lmp/pop_num)))
                ftimelog.write('    - average RDF time:               %.2f s\n'%(float(T_rdf/pop_num)))
                ftimelog.write('    - getting info time:              %.2f s\n'%(float(T_get/pop_num)))
                ftimelog.write('    - longest population              %.2f s\n'%(float(max(Tmax))))
                ftimelog.write('    - total time:                     %.2f s\n\n'%(t2-t1))

    # Post processing
        with open(output_path+"/timelog","a") as ftimelog:
            ftimelog.write('   Postprocessing\n')

        # write structures
        write_structures(poscars,contcars,output_path)
        global_best_results, Emin   = update_best_structures(global_best_results,results,inp_file,gen, pop_info)

        with open(src_path+"/../Emin","w") as f:
            f.write("%f"%(Emin))

        t3 = time.time()
        with open(output_path+"/timelog","a") as ftimelog:
                ftimelog.write('    - write structures:               %.2f s\n'%(t3-t2))
        
    # Post processing: check if same structure + antiseed

        if inp_file['similarity_metric']['type'] == 'not-use':
            same_info        = [[0,0] for i in range(pop_num)]
        else:
            if num_of_best_input < 1:  
                num_of_best = update_num_of_best(results, Emin, inp_file)
            
            struct_pair = get_pairs(pop_num, results, num_of_best, inp_file, pop_info)
            results, same_info, same_sentence =  check_if_same_structure_from_distance(results, struct_pair,inp_file,pop_num,gen,history,pop_info,output_path)

            with open(output_path+"/distance_info","a") as fdd:
                fdd.write(same_sentence)

        t4 = time.time()
        with open(output_path+"/timelog","a") as ftimelog:
                ftimelog.write('    - check similarity:               %.2f s (%d pairs)\n'%(t4-t3,len(struct_pair)))

        # calculate antiseed energy
        ANTISEED = make_antiseed(results, global_best_results, inp_file, Emin)
        t5 = time.time()
        with open(output_path+"/timelog","a") as ftimelog: #################################################3
                ftimelog.write('    - Calculate antiseed:             %.2f s\n'%(t5-t4))
    
        # etc
        f_individual     = open(output_path+'/totalGen','a')  # current generation
        f_best           = open(output_path+'/totalbest','a') # current generation

        history          = write_info(gen, results, pop_info, f_best, f_individual, same_info, history)   # write Individuals, BESTIndividuals and also write history, and modify results (increase energy of the same structure)
        f_individual.close()
        f_best.close()
        del f_individual
        del f_best
        # results          = check_max_kept_best(results,history,gen,max_iter) # simple metadynamics option
        
        
    
        if inp_file['similarity_metric']['type'] != 'not-use':
            global_best_results = remove_replica_from_global_best_results(gen, global_best_results, same_info, num_of_best_input)

        if num_of_best_input < 1:  
            num_of_best = update_num_of_best(results, Emin, inp_file)

        if gen > 2:
            del(history[gen-2])
       
        if gen != startnum:
            del average_atom_e 
        average_atom_e = update_average_atom_e(global_best_results,atomnumlist)

        with open(output_path+"/ave_atomic_e","a") as f_atom_e:
            f_atom_e.write(str(gen))
            for atom in range(len(atomnumlist)):
                f_atom_e.write(" ")
                f_atom_e.write(str(average_atom_e[atom]))
            f_atom_e.write("\n")


        
        # log files
        with open(output_path+"/BestResults","w") as fbestresults:# writing best results
            fbestresults.write("  Generation  Population    E    V  antiseed \n")
            fbestresults.write("---------------------------------------------\n")
            for i in range(len(global_best_results)):
                fbestresults.write(str(global_best_results[i][0]))
                fbestresults.write("    ")
                fbestresults.write(str(global_best_results[i][1]))
                fbestresults.write("    ")
                fbestresults.write("%.3f"%(global_best_results[i][2]))
                fbestresults.write("    ")
                fbestresults.write("%.3f"%(global_best_results[i][3]))
                fbestresults.write("    ")
                structtype = struct_type_func(history[global_best_results[i][0]][global_best_results[i][1]-1])
                fbestresults.write(structtype)
                fbestresults.write("   %.4f\n"%(ANTISEED[global_best_results[i][1]-1]))
        
        with open(output_path+"/best_history","a") as fbestresults:# writing best results
            for i in range(len(global_best_results)):
                fbestresults.write(str(global_best_results[i][0]))
                fbestresults.write("    ")
                fbestresults.write(str(global_best_results[i][1]))
                fbestresults.write("    ")
                fbestresults.write("%.3f"%(global_best_results[i][2]))
                fbestresults.write("    ")
                fbestresults.write("%.3f"%(global_best_results[i][3]))
                fbestresults.write("    ")
                structtype = struct_type_func(history[global_best_results[i][0]][global_best_results[i][1]-1])
                fbestresults.write(structtype)
                fbestresults.write("   %.4f\n"%(ANTISEED[global_best_results[i][1]-1]))
        

        # save old results
        
        del oldresults
        oldresults = {} 
        oldresults['E'] = [0.0 for i in range(pop_num)]
        oldresults['V'] = [0.0 for i in range(pop_num)]
        oldresults['coordinate'] = [[[0.0,0.0,0.0] for i in range(tot_atom_num)] for k in range(pop_num)]
        oldresults['lattice']    = [[0.0 for i in range(6)] for k in range(pop_num)]
        oldresults['rdf']        = [[0.0 for i in range(int(rdf_grid*num_of_atom*(num_of_atom+1)/2))] for k in range(pop_num)]
        oldresults['atom_e']     = [[0.0 for i in range(tot_atom_num)] for k in range(pop_num)]
 
        oldresults['E'] = results['E']
        oldresults['V'] = results['V']
        oldresults['coordinate'] = results['coordinate']
        oldresults['lattice']    = results['lattice']
        oldresults['rdf']        = results['rdf']
        oldresults['atom_e']     = results['atom_e']



        gc.collect()

        for process in range(1,corenum):
            comm.send([-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1], dest=process, tag=1)
        

        os.chdir(output_path+"/random_structure")
        for pop in range(1,random_num+1):
            os.rename("coo_"+str(pop)+"_in_advance", "coo_"+str(pop))

        t6 = time.time()

        with open(output_path+"/timelog","a") as ftimelog: #################################################3
                ftimelog.write('    - etc:                            %.2f s\n'%(t6-t5))
                ftimelog.write('    - total time:                     %.2f s\n\n'%(t6-t2))
   

    elif rank < corenum-1: # slave node
        while True:
            [gen, Emin, rt_start, rpop, rparent, rstruct_type, rVmin, parent_data, t_struct0, average_atom_e, startnum, job] = comm.recv(source = 0, tag = 1)
            if rpop == -1:
                gc.collect()
                break
            if job == 1:  ## lammps calculation
                t_struct1, t_lmp, t_rdf, t_get, E, V, poscar, contcar, rdfs, latt, coor, atom_e, status, attempt = core_job(gen, rank, rpop, rparent, rstruct_type, rVmin, inp_file, src_path, output_path, atom_names, atom_nums,parent_data,continue_, Emin, tolerance_matrix, average_atom_e, startnum)
                comm.send([rt_start, rank, rpop, t_struct1, t_lmp, t_rdf, t_get, E, V, poscar, contcar, rdfs, latt, coor, atom_e, status, 0], dest = 0, tag = 2)

            if job == 2:   ## similarity calculation
                struct1 = Emin
                struct2 = rt_start
                rdf1    = rpop
                rdf2    = rparent
                pop_num = rstruct_type
                pop_info= rVmin
                history = parent_data
                Einfo   = t_struct0
                Vinfo   = average_atom_e
                same_info_specific, distance = check_if_same_structure_from_distance(rdf1,rdf2,struct1,struct2,inp_file,pop_num,gen,history,pop_info, Einfo, Vinfo) 

                sigma  = inp_file['antiseed']['gaussian_width']
                ANTISEED_specific  = np.exp(-distance**2.0/2.0/sigma**2.0)

                with open(src_path+"/../inheriting","a") as f:
                    f.write("%d %d %d %f\n"%(gen,struct1,struct2,distance))       

                comm.send([rank, struct1, struct2, same_info_specific, ANTISEED_specific], dest = 0, tag = 4) 


    else:    # slave node who make random structures
        while True:
            [gen, Emin, rt_start, random_num, rparent, rstruct_type, rVmin, parent_data, t_struct0, average_atom_e, startnum, job] = comm.recv(source = 0, tag = 1)

            if job == 0:   ## random struture generation
                tr1 = time.time()
                os.chdir(output_path+"/random_structure")
                for pop in range(1,random_num+1):
                   temp1, temp2 = random_gen(rVmin, inp_file, atomnumlist, atomnamelist, tolerance_matrix, str(pop)+"_in_advance", gen+1, pop) 
                   del temp1, temp2
                tr2 = time.time()
                comm.send([tr2-tr1, rank, 0  ,       0,       0 ,    0,    0,   0, 0,    0,      0,     0,     0,   0,     0,         0, 2], dest = 0, tag = 2) # after job finished, send end massage: job = 2 

            elif job == 1: ## lammps calculation
                rpop = random_num
                t_struct1, t_lmp, t_rdf, t_get, E, V, poscar, contcar, rdfs, latt, coor, atom_e, status, attempt = core_job(gen, rank, rpop, rparent, rstruct_type, rVmin, inp_file, src_path, output_path, atom_names, atom_nums,parent_data,continue_, Emin, tolerance_matrix, average_atom_e, startnum)
                comm.send([rt_start, rank, rpop, t_struct1, t_lmp, t_rdf, t_get, E, V, poscar, contcar, rdfs, latt, coor, atom_e, status, 0], dest = 0, tag = 2)
            
            elif job == 2:   ## similarity calculation
                struct1 = Emin
                struct2 = rt_start
                rdf1    = random_num
                rdf2    = rparent
                pop_num = rstruct_type
                pop_info= rVmin
                history = parent_data
                Einfo   = t_struct0
                Vinfo   = average_atom_e

                same_info_specific, distance = check_if_same_structure_from_distance(rdf1,rdf2,struct1,struct2,inp_file,pop_num,gen,history,pop_info,Einfo,Vinfo) 

                sigma  = inp_file['antiseed']['gaussian_width']
                ANTISEED_specific  = np.exp(-distance**2.0/2.0/sigma**2.0)
                
                with open(src_path+"/../inheriting","a") as f:
                    f.write("%d %d %d %f\n"%(gen, struct1,struct2,distance))       

                comm.send([rank, struct1, struct2, same_info_specific, ANTISEED_specific], dest = 0, tag = 4) 
            else:
                gc.collect()
                break
            
