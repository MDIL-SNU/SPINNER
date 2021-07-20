#from operators import *
import os,random
import subprocess
import sys
import bisect
import numpy as np
from module_structure import *
import shutil
from time import sleep
import gc

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

def make_potcar_cutoff_10(output_path):
    with open(output_path+"/potential",'r') as f:
        freadlines = f.readlines()
    with open(output_path+"/potential10","w") as f:
        for fstr in freadlines:
            if 'POT' in fstr:
                fstr_separate = fstr.split()
                f.write(fstr_separate[0])
                f.write(" ")
                f.write(fstr_separate[1])
                f.write(" 15.0\n")
            else:
                f.write(fstr)
    del freadlines
    gc.collect()


def get_parent(pop_info,pop):
    
    struct_type = pop_info['struct_type'][pop-1]

    if struct_type == 0:
        parent = int(0)
    if struct_type == 1:
        parent = int(pop_info['con_permutation_num'][pop-1])
    if struct_type == 2:
        parent = int(pop_info['all_permutation_num'][pop-1])
    if struct_type == 3:
        parent = int(pop_info['latticemutation_num'][pop-1])
    if struct_type == 4:
        parent = int(pop_info['kept_best_num'][pop-1])
    if struct_type == 5:
        parent = pop_info['crossover_num'][pop-1]

    return struct_type, parent
    gc.collect()

def get_pop_info(inp_file, gen, pop_num, results, src_path, Emin, num_of_best, output_path, ANTISEED):

    atom_num = len(inp_file['material'])
   
    tot_atom_num = inp_file['tot_atom_num']
 
    # determination of volume of random structure
    if gen == 1:
        Vmin = inp_file['initial_volume']
    else:
        info = []
        for j in range(len(results['V'])):
            info.append([results['E'][j],results['V'][j]])
        volumes_sorted = sorted(info, key = lambda x : x[0])
        Vmin = volumes_sorted[0][1]
        del info

    #read from inp_file
    output_dir = inp_file['output_dir']
        
    struct_type = [0 for i in range(pop_num)] # 0: random, 1: conditional permutation, 2: all permutation 3: lattice permutation 4: keepbest 5: crossover
    operators_num   = [0 for i in range(6)]   # [0]: random, [1]: random topology, ...

    random           = float(inp_file['operator']['random'])
    con_permutation  = float(inp_file['operator']['con_permutation'])
    all_permutation  = float(inp_file['operator']['all_permutation'])
    latticemutation  = float(inp_file['operator']['latticemutation'])
    crossover        = float(inp_file['operator']['crossover'])

    total_frac  = random +  con_permutation + all_permutation + latticemutation + crossover

    pop_operator  = pop_num - num_of_best

    # decide the number of operators
    if gen == 1:
        operators_num[0]      = pop_num
        con_permutation_nums  = [0]*pop_num
        all_permutation_nums  = [0]*pop_num
        latticemutation_nums  = [0]*pop_num
        kept_best_nums        = [0]*pop_num
        crossover_nums        = [[0,0]]*pop_num

        struct_types          = [0]*operators_num[0]

    else:
        operators_num[1] = int(round(con_permutation/total_frac * float(pop_operator)))
        operators_num[2] = int(round(all_permutation/total_frac * float(pop_operator)))
        operators_num[3] = int(round(latticemutation/total_frac * float(pop_operator)))
        operators_num[4] = num_of_best
        operators_num[5] = int(round(crossover/total_frac * float(pop_operator)))
        
        operators_num[0] = pop_operator - sum(operators_num[0:4]) - operators_num[5]

        # decide what structures to use (random / or consider energy)
        struct_types          = [0]*operators_num[0] + [1]*operators_num[1] + [2]*operators_num[2] + [3]*operators_num[3] + [4]*operators_num[4] + [5]*operators_num[5] 
        con_permutation_nums  = decide_con_permutation(struct_type, operators_num, results, inp_file, Emin, tot_atom_num, ANTISEED)
        all_permutation_nums  = decide_all_permutation(struct_type, operators_num, results, inp_file, Emin, tot_atom_num, ANTISEED)
        latticemutation_nums  = decide_latticemutation(struct_type, operators_num, results, inp_file, Emin, tot_atom_num, ANTISEED)
        kept_best_nums        = decide_kept_best(struct_type, operators_num, results, inp_file, num_of_best, Emin)
        crossover_nums        = decide_crossover(struct_type, operators_num, results, inp_file, Emin, tot_atom_num, ANTISEED)

    pop_info = {}
    pop_info['struct_type']          = struct_types
    pop_info['con_permutation_num']  = con_permutation_nums
    pop_info['all_permutation_num']  = all_permutation_nums
    pop_info['latticemutation_num']  = latticemutation_nums
    pop_info['kept_best_num']        = kept_best_nums
    pop_info['crossover_num']        = crossover_nums

    del operators_num
    operators_num   = [0 for i in range(6)]   # [0]: random, [1]: random topology, ...
    operators_num[1] = int(round(con_permutation/total_frac * float(pop_operator)))
    operators_num[2] = int(round(all_permutation/total_frac * float(pop_operator)))
    operators_num[3] = int(round(latticemutation/total_frac * float(pop_operator)))
    operators_num[4] = num_of_best
    operators_num[5] = int(round(crossover/total_frac * float(pop_operator)))
    
    operators_num[0] = pop_operator - sum(operators_num[0:4]) - operators_num[5]
    random_number            = operators_num[0]
    del operators_num, con_permutation_nums, all_permutation_nums, latticemutation_nums, kept_best_nums, crossover_nums, struct_types
    gc.collect()

    return Vmin, pop_info, random_number


################ decide structure number randomly ####################

def decide_con_permutation(struct_type, operators_num, results, inp_file, Emin, tot_atom_num, ANTISEED):

    Ecut  = inp_file['energy_criteria']['energy_cut_for_inheriting_structures']

    bestnum = 0
    for i in range(len(results['E'])):
        if results['E'][i]-Emin < Ecut*float(tot_atom_num):
            bestnum += 1

    sigma = inp_file['antiseed']['selection_gaussian']

    mutationnum = [0 for i in range(len(struct_type))]

    # place in order
    sumEN = []
    for i in range(len(results['E'])):
        sumEN.append([i+1,results['E'][i],np.exp(-ANTISEED[i]**2.0/2/sigma**2.0)])
    sumEN = sorted(sumEN, key = lambda x : x[1])

    population = []
    weights = []
    for i in range(len(sumEN)):
        population.append(sumEN[i][0])
        weights.append(sumEN[i][2])

    for i in range(sum(operators_num[0:1]),sum(operators_num[0:2])):
          if random.random() > inp_file['antiseed']['selection_fraction']:
            mutationnum[i] = sumEN[random.randint(0,bestnum-1)][0]
          else:
            mutationnum[i] = choice(population,weights) 
    del sumEN

    gc.collect()
    return mutationnum

def decide_all_permutation(struct_type, operators_num, results, inp_file,Emin,tot_atom_num, ANTISEED):

    Ecut  = inp_file['energy_criteria']['energy_cut_for_inheriting_structures']

    bestnum = 0
    for i in range(len(results['E'])):
        if results['E'][i]-Emin < Ecut*float(tot_atom_num):
            bestnum += 1

    sigma = inp_file['antiseed']['selection_gaussian']

    mutationnum = [0 for i in range(len(struct_type))]

    # place in order
    sumEN = []
    for i in range(len(results['E'])):
        sumEN.append([i+1,results['E'][i],np.exp(-ANTISEED[i]**2.0/2/sigma**2.0)])
    sumEN = sorted(sumEN, key = lambda x : x[1])

    population = []
    weights = []
    for i in range(len(sumEN)):
        population.append(sumEN[i][0])
        weights.append(sumEN[i][2])

    for i in range(sum(operators_num[0:2]),sum(operators_num[0:3])):
          if random.random() > inp_file['antiseed']['selection_fraction']:
            mutationnum[i] = sumEN[random.randint(0,bestnum-1)][0]
          else:
            mutationnum[i] = choice(population,weights) 
    del sumEN

    gc.collect()
    return mutationnum
    
def decide_latticemutation(struct_type, operators_num, results, inp_file,Emin,tot_atom_num, ANTISEED):

    Ecut  = inp_file['energy_criteria']['energy_cut_for_inheriting_structures']

    bestnum = 0
    for i in range(len(results['E'])):
        if results['E'][i]-Emin < Ecut*float(tot_atom_num):
            bestnum += 1

    sigma = inp_file['antiseed']['selection_gaussian']

    mutationnum = [0 for i in range(len(struct_type))]

    # place in order
    sumEN = []
    for i in range(len(results['E'])):
        sumEN.append([i+1,results['E'][i],np.exp(-ANTISEED[i]**2.0/2/sigma**2.0)])
    sumEN = sorted(sumEN, key = lambda x : x[1])

    population = []
    weights = []
    for i in range(len(sumEN)):
        population.append(sumEN[i][0])
        weights.append(sumEN[i][2])

    for i in range(sum(operators_num[0:3]),sum(operators_num[0:4])):
          if random.random() > inp_file['antiseed']['selection_fraction']:
            mutationnum[i] = sumEN[random.randint(0,bestnum-1)][0]
          else:
            mutationnum[i] = choice(population,weights) 
    del sumEN

    gc.collect()
    return mutationnum

def decide_kept_best(struct_type, operators_num, results, inp_file, num_of_best, Emin):
   
    num_of_best_input = inp_file['structure']['num_of_best'] 

    # place in order
    sumEN = []
    for i in range(len(results['E'])):
        sumEN.append([i+1,results['E'][i]])
    sumEN = sorted(sumEN, key = lambda x : x[1])


    mutationnum = [0 for i in range(len(struct_type))]

    for i in range(sum(operators_num[0:4]),sum(operators_num[0:5])):
        mutationnum[i] = sumEN[i-sum(operators_num[0:4])][0]

    del sumEN
    gc.collect()
    return mutationnum

def decide_crossover(struct_type, operators_num, results, inp_file,Emin,tot_atom_num, ANTISEED):

    Ecut  = inp_file['energy_criteria']['energy_cut_for_inheriting_structures']

    bestnum = 0
    for i in range(len(results['E'])):
        if results['E'][i]-Emin < Ecut*float(tot_atom_num):
            bestnum += 1

    sigma = inp_file['antiseed']['selection_gaussian']

    mutationnum = [0 for i in range(len(struct_type))]

    # place in order
    sumEN = []
    for i in range(len(results['E'])):
        sumEN.append([i+1,results['E'][i],np.exp(-ANTISEED[i]**2.0/2/sigma**2.0)])
    sumEN = sorted(sumEN, key = lambda x : x[1])

    population = []
    weights = []
    for i in range(len(sumEN)):
        population.append(sumEN[i][0])
        weights.append(sumEN[i][2])

    #listnum = [i for i in range(bestnum)]

    for i in range(sum(operators_num[0:5]),sum(operators_num[0:6])):
        if random.random() > inp_file['antiseed']['selection_fraction']:
           ran1 = sumEN[random.randint(0,bestnum-1)][0]

        else:
           ran1 = choice(population,weights) 

        if random.random() > 0.5:
           ran2 = ran1
        else:
          if random.random() > inp_file['antiseed']['selection_fraction']:
           ran2 = sumEN[random.randint(0,bestnum-1)][0]

          else:
           ran2 = choice(population,weights) 

        mutationnum[i] = [ran1, ran2]

    del sumEN
    gc.collect()
    return mutationnum
