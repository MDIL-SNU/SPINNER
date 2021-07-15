import os
import shutil
import numpy as np
import gc
def write_structures(poscars,contcars,output_path):

    with open(output_path+"/POSCARs","a") as fw:
        for i in range(len(poscars)):
            fw.write(poscars[i])
    
    with open(output_path+"/CONTCARs","a") as fw:
        for i in range(len(contcars)):
            fw.write(contcars[i])
            
    gc.collect()

def struct_type_func(history):
    if history[4] == 0:
        return "random           "
    if history[4] == 1:
        return "c_permutation    "
    if history[4] == 2:
        return "a_permutation    "
    if history[4] == 3:
        return "lattice_mutation "
    if history[4] == 5:
        return "crossover        "

def remove_replica_from_global_best_results(gen, global_best_results, same_info, num_of_best_input):

    remove_array = []

    for i in range(len(global_best_results)):

        if (global_best_results[i][0] == gen)  and  (same_info[global_best_results[i][1]-1] != [0,0]):

            remove_array.append([global_best_results[i][0], global_best_results[i][1], global_best_results[i][2], global_best_results[i][3], global_best_results[i][4]])

    for i in range(len(remove_array)):
            
        global_best_results.remove(remove_array[i])

    del remove_array
    return global_best_results
    gc.collect()

def update_best_structures(global_best_results, results, inp_file, gen, pop_info):
    
    num_of_best_input = inp_file['structure']['num_of_best']
    tot_atom_num      = inp_file['tot_atom_num']

    # update num_of_best
    
    if gen == 1:
        if num_of_best_input < 1:

            Emin_limit = inp_file['energy_criteria']['energy_cut_for_best_structures']*float(tot_atom_num)

            Emin = min(results['E'])

            num_of_best = 0

            for i in range(len(results['E'])):

                if results['E'][i] - Emin < Emin_limit:
                    global_best_results.append([gen,i+1,results['E'][i],results['V'][i],results['atom_e'][i]])
                    num_of_best = num_of_best + 1

        else:
            info = []
            for i in range(len(results['E'])):
                info.append([i,results['E'][i]])

            info_sorted = sorted(info, key = lambda x : x[1])

            Emin = info_sorted[0][1]

            del info,info_sorted

            for i in range(num_of_best_input):
                tempi = info_sorted[i][0]
                global_best_results.append([gen, tempi+1, results['E'][tempi], results['V'][tempi],results['atom_e'][i]])

            num_of_best = num_of_best_input

    else:

        tot_results = []
        for i in range(len(results['E'])):
            tot_results.append([gen,i+1,results['E'][i],results['V'][i],results['atom_e'][i]])
#
#        THIS PART IS NECESSARY WHEN ANTISEED OPTION IS RUNNING: NOT AVAILABLE FOR NOW (have to be rewritten)
#        for i in range(len(global_best_results)):
#            kept_best = 0
#
#            if global_best_results[i][0] == gen-1:
#
#                for j in range(len(pop_info['kept_best_num'])):
#
#                    if (pop_info['kept_best_num'][j] == global_best_results[i][1]) and (results['E'][j] < 9999.0):
#
#                        kept_best = kept_best + 1
#
#            if kept_best == 0:
#                tot_results.append([global_best_results[i][0],global_best_results[i][1],global_best_results[i][2],global_best_results[i][3]])

        if num_of_best_input < 1:

            Emin_limit = inp_file['energy_criteria']['energy_cut_for_best_structures']*float(tot_atom_num)

            Emin = sorted(tot_results, key = lambda x : x[2])[0][2]

            del global_best_results
            global_best_results = []

            num_of_best = 0

            for i in range(len(tot_results)):

                if tot_results[i][2] - Emin < Emin_limit:
                    global_best_results.append([tot_results[i][0], tot_results[i][1], tot_results[i][2], tot_results[i][3], tot_results[i][4]])
                    num_of_best = num_of_best + 1
            
        else:
            info = []
            info_sorted = sorted(tot_results, key = lambda x : x[2])
            
            del global_best_results
            global_best_results = []

            for i in range(num_of_best_input):
                tempi = info_sorted[i][0]
                global_best_results.append([info_sorted[i][0],info_sorted[i][1],info_sorted[i][2],info_sorted[i][3],info_sorted[i][4]])

            Emin = info_sorted[0][2]
            
            del info,info_sorted
            num_of_best = num_of_best_input

    # sort global_best_results

    global_best_results = sorted(global_best_results, key = lambda x : x[2])

    if gen != 1:
        del tot_results
    return global_best_results, Emin
    gc.collect()


def get_pairs(pop_num, results, num_of_best, inp_file, pop_info):

    rerelax_best      = inp_file['structure']['re-relax_best']

    tot_atom_num   = inp_file['tot_atom_num']
    atom_num       = len(inp_file['material'])
    E_cut          = inp_file['similarity_metric']['energy_cut']*float(tot_atom_num)
    V_cut          = inp_file['similarity_metric']['volume_cut']*float(tot_atom_num)
    
    Edat = []
    for i in range(len(results['E'])):
        Edat.append([i+1,results['E'][i],results['V'][i]])

    Edat = sorted(Edat, key = lambda x : x[1])

    count = 0

    struct_pair = []

    for i in range(num_of_best):

        for j in range(i+1, num_of_best):

            pop1 = Edat[i][0]
            E1   = Edat[i][1]
            V1   = Edat[i][2]
            pop2 = Edat[j][0]
            E2   = Edat[j][1]
            V2   = Edat[j][2]

            if rerelax_best == 0:
                if (abs(E1-E2) < E_cut) and (abs(V1-V2) < V_cut) and (pop_info['struct_type'][pop1-1] != 4 or pop_info['struct_type'][pop2-1] != 4):
                    struct_pair.append([pop1,pop2])
            else:
                if (abs(E1-E2) < E_cut) and (abs(V1-V2) < V_cut):
                    struct_pair.append([pop1,pop2])

    del Edat

    return struct_pair
    gc.collect()

def calculate_rdf_vector_distance(results,pop1,pop2):
    rdf1 = results['rdf'][pop1-1]
    rdf2 = results['rdf'][pop2-1]

    F12 = 0.0
    F11 = 0.0
    F22 = 0.0
    
    for i in range(len(rdf1)):
        F11 += rdf1[i]*rdf1[i]
        F22 += rdf2[i]*rdf2[i]
        F12 += rdf1[i]*rdf2[i]

    distance = 0.5*(1.0 - F12/np.sqrt(F11*F22))


    del rdf1, rdf2
    #print (distance)
    return distance
    gc.collect()

def check_if_same_structure_from_distance(results, struct_pair,inp_file,pop_num,gen,history,pop_info, output_path): 

    pair_num = len(struct_pair)
    distance_limit = inp_file['similarity_metric']['limit']

    same_info = [[0,0] for i in range(pop_num)]

    sentence = "" 

    for i in range(pair_num):
        pop1 = struct_pair[i][0]
        pop2 = struct_pair[i][1]

        distance = calculate_rdf_vector_distance(results,pop1,pop2) 

        sentence += str(gen)+" "+str(pop1)+" "+str(pop2)+" "+str(distance)+"\n"

        if distance < distance_limit:
                if (pop_info['kept_best_num'][pop1-1] != 0) and (pop_info['kept_best_num'][pop2-1] == 0):
                        same_info[pop2-1][0] = gen
                        same_info[pop2-1][1] = pop1
                        results['E'][pop2-1] = 10000.0
                elif (pop_info['kept_best_num'][pop1-1] != 0) and (pop_info['kept_best_num'][pop2-1] == 0):
                        same_info[pop1-1][0] = gen
                        same_info[pop1-1][1] = pop2
                        results['E'][pop1-1] = 10000.0
                elif (pop_info['kept_best_num'][pop1-1] != 0) and (pop_info['kept_best_num'][pop2-1] != 0):
                    
                    previous_num1 = pop_info['kept_best_num'][pop1-1]
                    origin_gen1   = history[gen-1][previous_num1-1][5]
                    previous_num2 = pop_info['kept_best_num'][pop2-1]
                    origin_gen2   = history[gen-1][previous_num2-1][5]

                    if origin_gen1 > origin_gen2:
                        same_info[pop1-1][0] = gen
                        same_info[pop1-1][1] = pop2
                        results['E'][pop1-1] = 10000.0
                    else:
                        same_info[pop2-1][0] = gen
                        same_info[pop2-1][1] = pop1
                        results['E'][pop2-1] = 10000.0 
                else: 
                    same_info[pop1-1][0] = gen
                    same_info[pop1-1][1] = pop2
                    results['E'][pop1-1] = 10000.0

                    
    return results, same_info, sentence
    gc.collect()

def write_info(gen, results, pop_info, f_best, f_individual, same_info, history):
    structure_type = ['random       ',
                      'c_permutation',
                      'a_permutation',
                      'lat_mutation ',
                      'keptbest     ',
                      'crossover    ']
    new_history = []

    if gen == 1:
        #individual

        info = []
        for j in range(len(results['E'])):
            if same_info[j] != [0,0]:
                f_individual.write("       %3d          %3d       %.3f         %.3f           same_as_%d_%d          -                 -\n"%(gen,j+1,results['E'][j],results['V'][j],same_info[j][0],same_info[j][1]))
                new_history.append([gen,j+1,results['E'][same_info[j][1]-1],results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   # 0,0: keepbest info
            else:
                f_individual.write("       %3d          %3d       %.3f         %.3f           %s      -                 -\n"%(gen,j+1,results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]]))
                new_history.append([gen,j+1,results['E'][j],results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   # 0,0: keepbest info
 
            info.append([j+1, results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]]])

        info = sorted(info, key = lambda x : x[1])

        #best individual
        f_best.write("       %3d          %3d       %.3f         %.3f           %s      -                 -\n"%(gen,info[0][0],info[0][1],info[0][2],info[0][3]))
            
    else:
     # if track info of keepbest
        info = []
        for j in range(len(results['E'])):
            if structure_type[pop_info['struct_type'][j]] == 'keptbest     ':
                previous_num = pop_info['kept_best_num'][j]
                origin_gen   = history[gen-1][previous_num-1][5]
                origin_pop   = history[gen-1][previous_num-1][6]
                origin_type  = history[gen-1][previous_num-1][4]
                keptbestnum  = history[gen-1][previous_num-1][8] + 1

                if same_info[j] != [0,0]:
                    if same_info[j][0] == gen:
                        Esame = results['E'][same_info[j][1]-1]
                    else:
                        Esame = history[same_info[j][0]][same_info[j][1]-1][2]
                    f_individual.write("       %3d          %3d       %.3f         %.3f           same_as_%d_%d      -               (%d,%d)  \n"%(gen,j+1,Esame,results['V'][j],same_info[j][0],same_info[j][1], origin_gen, origin_pop))
                    new_history.append([gen,j+1,Esame,results['V'][j],origin_type,origin_gen,origin_pop,same_info[j],keptbestnum])
                else:
                    f_individual.write("       %3d          %3d       %.3f         %.3f           %s      -                (%d,%d)      \n"%(gen,j+1,results['E'][j],results['V'][j],structure_type[origin_type], origin_gen, origin_pop))
                    new_history.append([gen,j+1,results['E'][j],results['V'][j],origin_type,origin_gen,origin_pop,same_info[j],keptbestnum])   

            elif structure_type[pop_info['struct_type'][j]] == 'c_permutation':
                if same_info[j] != [0,0]:
                    if same_info[j][0] == gen:
                        Esame = results['E'][same_info[j][1]-1]
                    else:
                        Esame = history[same_info[j][0]][same_info[j][1]-1][2]
                    f_individual.write("       %3d          %3d       %.3f         %.3f           same_as_%d_%d     (%d)               -\n"%(gen,j+1,Esame,results['V'][j],same_info[j][0],same_info[j][1],pop_info['con_permutation_num'][j]))
                    new_history.append([gen,j+1,Esame,results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   
                else:
                    f_individual.write("       %3d          %3d       %.3f         %.3f           %s     (%d)                -\n"%(gen,j+1,results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]],pop_info['con_permutation_num'][j]))
                    new_history.append([gen,j+1,results['E'][j],results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   


            elif structure_type[pop_info['struct_type'][j]] == 'a_permutation':
                if same_info[j] != [0,0]:
                    if same_info[j][0] == gen:
                        Esame = results['E'][same_info[j][1]-1]
                    else:
                        Esame = history[same_info[j][0]][same_info[j][1]-1][2]
                    f_individual.write("       %3d          %3d       %.3f         %.3f           same_as_%d_%d      (%d)               -\n"%(gen,j+1,Esame,results['V'][j],same_info[j][0],same_info[j][1],pop_info['all_permutation_num'][j]))
                    new_history.append([gen,j+1,Esame,results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   
                else:
                    f_individual.write("       %3d          %3d       %.3f         %.3f           %s     (%d)                -\n"%(gen,j+1,results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]],pop_info['all_permutation_num'][j]))
                    new_history.append([gen,j+1,results['E'][j],results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   

            elif structure_type[pop_info['struct_type'][j]] == 'lat_mutation ':
                if same_info[j] != [0,0]:
                    if same_info[j][0] == gen:
                        Esame = results['E'][same_info[j][1]-1]
                    else:
                        Esame = history[same_info[j][0]][same_info[j][1]-1][2]
                    f_individual.write("       %3d          %3d       %.3f         %.3f           same_as_%d_%d      (%d)               -\n"%(gen,j+1,Esame,results['V'][j],same_info[j][0],same_info[j][1],pop_info['latticemutation_num'][j]))
                    new_history.append([gen,j+1,Esame,results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   
                else:
                    f_individual.write("       %3d          %3d       %.3f         %.3f           %s     (%d)                -\n"%(gen,j+1,results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]],pop_info['latticemutation_num'][j]))
                    new_history.append([gen,j+1,results['E'][j],results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   

            elif structure_type[pop_info['struct_type'][j]] == 'softmutation ':
                if same_info[j] != [0,0]:
                    if same_info[j][0] == gen:
                        Esame = results['E'][same_info[j][1]-1]
                    else:
                        Esame = history[same_info[j][0]][same_info[j][1]-1][2]
                    f_individual.write("       %3d          %3d       %.3f         %.3f           same_as_%d_%d      (%d)               -\n"%(gen,j+1,Esame,results['V'][j],same_info[j][0],same_info[j][1],pop_info['softmutation_num'][j]))
                    new_history.append([gen,j+1,Esame,results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   
                else:
                    f_individual.write("       %3d          %3d       %.3f         %.3f           %s     (%d)                -\n"%(gen,j+1,results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]],pop_info['softmutation_num'][j]))
                    new_history.append([gen,j+1,results['E'][j],results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   

            elif structure_type[pop_info['struct_type'][j]] == 'random       ':
                if same_info[j] != [0,0]:
                    if same_info[j][0] == gen:
                        Esame = results['E'][same_info[j][1]-1]
                    else:
                        Esame = history[same_info[j][0]][same_info[j][1]-1][2]
                    f_individual.write("       %3d          %3d       %.3f         %.3f           same_as_%d_%d           -                 -\n"%(gen,j+1,Esame,results['V'][j],same_info[j][0],same_info[j][1]))
                    new_history.append([gen,j+1,Esame,results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   # 0,0: keepbest info
                else:
                    f_individual.write("       %3d          %3d       %.3f         %.3f           %s      -                 -\n"%(gen,j+1,results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]]))
                    new_history.append([gen,j+1,results['E'][j],results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   # 0,0: keepbest info
            
            elif structure_type[pop_info['struct_type'][j]] == 'crossover    ':
                if same_info[j] != [0,0]:
                    if same_info[j][0] == gen:
                        Esame = results['E'][same_info[j][1]-1]
                    else:
                        Esame = history[same_info[j][0]][same_info[j][1]-1][2]
                    f_individual.write("       %3d          %3d       %.3f         %.3f           same_as_%d_%d      (%d,%d)            -\n"%(gen,j+1,Esame,results['V'][j],same_info[j][0],same_info[j][1],pop_info['crossover_num'][j][0],pop_info['crossover_num'][j][1]))
                    new_history.append([gen,j+1,Esame,results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   
                else:
                    f_individual.write("       %3d          %3d       %.3f         %.3f           %s     (%d,%d)             -\n"%(gen,j+1,results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]],pop_info['crossover_num'][j][0],pop_info['crossover_num'][j][1]))
                    new_history.append([gen,j+1,results['E'][j],results['V'][j],pop_info['struct_type'][j],gen,j+1,same_info[j],1])   

            if same_info[j] == [0,0]:
                info.append([j+1, results['E'][j],results['V'][j],structure_type[pop_info['struct_type'][j]]])

        info = sorted(info, key = lambda x : x[1])

        #best individual
        if 'random' in info[0][3]:
            f_best.write("       %3d          %3d       %.3f         %.3f           %s      -                 -\n"%(gen,info[0][0],info[0][1],info[0][2],info[0][3]))
        if 'c_permutation' in info[0][3]:
            f_best.write("       %3d          %3d       %.3f         %.3f           %s     (%d)                -\n"%(gen,info[0][0],info[0][1],info[0][2],info[0][3],pop_info['con_permutation_num'][info[0][0]-1]))
        if 'a_permutation' in info[0][3]:
            f_best.write("       %3d          %3d       %.3f         %.3f           %s     (%d)                -\n"%(gen,info[0][0],info[0][1],info[0][2],info[0][3],pop_info['all_permutation_num'][info[0][0]-1]))
        if 'lat_mutation' in info[0][3]:
            f_best.write("       %3d          %3d       %.3f         %.3f           %s     (%d)                -\n"%(gen,info[0][0],info[0][1],info[0][2],info[0][3],pop_info['latticemutation_num'][info[0][0]-1]))
        if 'crossover'    in info[0][3]:
            f_best.write("       %3d          %3d       %.3f         %.3f           %s     (%d,%d)             -\n"%(gen,info[0][0],info[0][1],info[0][2],info[0][3],pop_info['crossover_num'][info[0][0]-1][0],pop_info['crossover_num'][info[0][0]-1][1]))
        if 'keptbest' in info[0][3]:
            previous_num = pop_info['kept_best_num'][info[0][0]-1]
            origin_gen   = history[gen-1][previous_num-1][5]
            origin_pop   = history[gen-1][previous_num-1][6]
            origin_type  = history[gen-1][previous_num-1][4]

            f_best.write("       %3d          %3d       %.3f         %.3f           %s      -                (%d,%d)      \n"%(gen,info[0][0],info[0][1],info[0][2],structure_type[origin_type], origin_gen, origin_pop))

        
    history[gen] = new_history
    del new_history

    return history
    gc.collect()

#print (write_info(gen, results, pop_info, f_best, f_individual, same_info, history))


def check_max_kept_best(results,history,gen,max_iter):
    for i in range(len(results['E'])):
        if history[gen][i][8] > max_iter:
            results['E'][i] = 10000.0

    return results
    gc.collect()


def update_num_of_best(results, Emin, inp_file):
    tot_atom_num = inp_file['tot_atom_num']
    Elimit       = float(tot_atom_num) * inp_file['energy_criteria']['energy_cut_for_best_structures']

    num_of_best = 0
    for i in range(len(results['E'])):
        if results['E'][i] - Emin < Elimit:
            num_of_best = num_of_best + 1

    return num_of_best

    gc.collect()


def update_average_atom_e(global_best_results, atomnumlist):

    average_atom_e = [0.0 for i in range(len(atomnumlist))]

    num_of_best_struct = len(global_best_results)
    num_of_elements    = len(atomnumlist)

    for i in range(num_of_best_struct):

        #calculate average atomic energy of structure i

        temp_i = -1
        for j in range(num_of_elements): 

            for k in range(atomnumlist[j]):

                temp_i += 1

                average_atom_e[j] += global_best_results[i][4][temp_i]

    for j in range(num_of_elements):

        average_atom_e[j] = average_atom_e[j] / atomnumlist[j] / num_of_best_struct 

    return average_atom_e

    gc.collect()

def make_antiseed(results, global_best_results, inp_file, Emin):
    
    Ecut  = inp_file['energy_criteria']['energy_cut_for_inheriting_structures']

    tot_atom_num = inp_file['tot_atom_num']
    
    ANTISEED = [10000.0 for i in range(len(results['E']))]

    nums = 0.0

    if inp_file['antiseed']['activate_antiseed'] == True:

        sigma  = inp_file['antiseed']['gaussian_width']

        for struct1 in range(1,len(results['E'])+1):

            if results['E'][struct1-1]-Emin < Ecut*float(tot_atom_num):

              nums += 1.0

              ANTISEED[struct1-1] = 0.0

              for struct2 in range(1,len(results['E'])+1):
                if (struct1 != struct2)  and results['E'][struct2-1]-Emin < Ecut*float(tot_atom_num):

                   distance = calculate_rdf_vector_distance(results,struct1,struct2) 
                   ANTISEED[struct1-1] += np.exp(-distance**2.0/2/sigma**2.0)

        for i in range(len(ANTISEED)):
            if ANTISEED[i] != 10000.0:
                ANTISEED[i] = float(ANTISEED[i] / nums)

    else:

        for struct1 in range(1,len(results['E'])+1):

            if results['E'][struct1-1]-Emin < Ecut*float(tot_atom_num):

              ANTISEED[struct1-1] = 0.0


    return ANTISEED
