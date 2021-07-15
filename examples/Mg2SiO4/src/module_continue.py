import os,sys,subprocess,shutil
import re

def update_average_atom_from_file(continue_, original_path, output_path):
    with open(original_path+"/ave_atomic_e","r") as f:
        f_ave = f.readlines()

    with open(output_path+"/ave_atomic_e","w") as f:
        for i in range(continue_):
            f.write(f_ave[i])

    aveE = list(map(float,f_ave[continue_-1].split()[1:len(f_ave[continue_-1])]))

    return aveE

def rewrite_contcars_poscars(output_path, original_path, continue_, tot_atom_num):

    with open(original_path+"/POSCARs", "r") as f:
        poscar_orig = f.readlines()

    with open(original_path+"/CONTCARs", "r") as f:
        contcar_orig = f.readlines()

    poscar_line = 0

    for i in range(len(poscar_orig)):
        if i%(tot_atom_num+8) ==0:
            if int(poscar_orig[i].split()[1]) <= continue_:
                poscar_line = i + tot_atom_num + 8

    with open(output_path+"/CONTCARs", "w") as f:
        for i in range(poscar_line):
            f.write(contcar_orig[i])
            
    with open(output_path+"/POSCARs", "w") as f:
        for i in range(poscar_line):
            f.write(poscar_orig[i])

def getting_global_best_results(continue_, inp_file, src_path, output_dir):

    num_of_best_input = inp_file['structure']['num_of_best']

    os.chdir(src_path+"/../"+output_dir)

    with open("totalGen","r") as f:
        freadlines = f.readlines()
    
    linenum = len(freadlines)

    best_structures = []

    num_of_best = 0

    Es =[]

    generation_info = [] 

    for i in range(2,linenum):
        gen         = int(freadlines[i].split()[0])
        pop         = int(freadlines[i].split()[1])
        E           = float(freadlines[i].split()[2])
        V           = float(freadlines[i].split()[3])
        struct_type = str(freadlines[i].split()[4])

        if gen == continue_:
            generation_info = [gen,pop,E,V] 



        best_structures.append([gen, pop, E, V])

        Es.append(E)

    Emin = min(Es)

    return best_structures, Emin

def getting_info_from_files(continue_, tot_atom_num, original_path, src_path, num_of_atom):
    os.chdir(original_path)
    with open("totalGen","r") as f:
        freadlines = f.readlines()
    linenum = len(freadlines)

    results = {}
    results['E'] = []
    results['V'] = []
    results['coordinate'] = []
    results['lattice'] = []
    
    linelist = []
    pop_info = []

    for i in range(2,linenum):
        gen         = int(freadlines[i].split()[0])
        pop         = int(freadlines[i].split()[1])
        E           = float(freadlines[i].split()[2])
        V           = float(freadlines[i].split()[3])
        struct_type = str(freadlines[i].split()[4])
        keptbest    = str(freadlines[i].split()[6])
 
        if 'same' in struct_type:
            E = 10000.0

        if gen == continue_:
            results['E'].append(E)
            results['V'].append(V)
            linelist.append((i-2)*(tot_atom_num+8))
            
        if struct_type == 'random':
            struct_num = 0
        elif struct_type == 'c_permutation':
            struct_num = 1
        elif struct_type == 'a_permutation':
            struct_num = 2
        elif struct_type == 'lat_mutation':
            struct_num = 3
        elif struct_type == 'crossover':
            struct_num = 5 
        elif 'same' in struct_type:
            struct_num = 0
        
        if keptbest != '-':
            struct_num = 4

        pop_info.append(struct_type)

    with open("CONTCARs","r") as f:
        CONTCARreadlines = f.readlines()
    
    for i in range(len(linelist)):
        line = linelist[i]
        latt = []
        latt.append(float(CONTCARreadlines[line+2].split()[0]))
        latt.append(float(CONTCARreadlines[line+3].split()[1]))
        latt.append(float(CONTCARreadlines[line+4].split()[2]))
        latt.append(float(CONTCARreadlines[line+3].split()[0]))
        latt.append(float(CONTCARreadlines[line+4].split()[0]))
        latt.append(float(CONTCARreadlines[line+4].split()[1]))

        results['lattice'].append(latt)

        coor = []
        for j in range(line+8,line+8+tot_atom_num):
            coor1 = float(CONTCARreadlines[j].split()[0])
            coor2 = float(CONTCARreadlines[j].split()[1])
            coor3 = float(CONTCARreadlines[j].split()[2])
            coor.append([coor1,coor2,coor3])
        results['coordinate'].append(coor)

    results['rdf']        = [[0.0 for I in range(int(250*num_of_atom*(num_of_atom+1)/2))] for K in range(len(linelist))]
    results['atom_e']     = [[0.0 for I in range(tot_atom_num)] for K in range(len(linelist))]
    
    return results, pop_info, linelist

def getting_history(src_path,original_path, continue_):
    os.chdir(original_path)
    with open("totalGen","r") as f:
        freadlines = f.readlines()

    linenum = len(freadlines)

    historyold = {}

    for i in range(2,linenum):
        gen         = int(freadlines[i].split()[0])
        pop         = int(freadlines[i].split()[1])
        E           = float(freadlines[i].split()[2])
        V           = float(freadlines[i].split()[3])
        struct_type = str(freadlines[i].split()[4])
        parents     = str(freadlines[i].split()[5])
        keptbest    = str(freadlines[i].split()[6])

        if gen <= continue_ and continue_-2 <= gen:

            if struct_type == 'random':
                struct_num = 0
            elif struct_type == 'c_permutation':
                struct_num = 1
            elif struct_type == 'a_permutation':
                struct_num = 2
            elif struct_type == 'lat_mutation':
                struct_num = 3
            elif struct_type == 'softmutation':
                struct_num = 4
            elif struct_type == 'crossover':
                struct_num = 5


            if keptbest == '-':
                origin_gen = gen
                origin_pop = pop
                keep_num = 1
            else:
                temp1 = keptbest.split(',')[0]
                temp2 = keptbest.split(',')[1]
                origin_gen = int(temp1[1:len(temp1)])
                origin_pop = int(temp2[0:len(temp2)-1])

                keep_num   = gen - origin_gen + 1

            if 'same' in struct_type:
                temp = struct_type.split('_')
                same_info = [int(temp[2]),int(temp[3])]

            else:
                same_info = [0,0]

            if gen in historyold:
                historyold[gen].append([gen,pop,E,V,struct_num,origin_gen,origin_pop,same_info,keep_num])
            else:
                historyold[gen] = []
                historyold[gen].append([gen,pop,E,V,struct_num,origin_gen,origin_pop,same_info,keep_num])

    return historyold


def penalty_if_same_structure_or_max_iter(continue_, results, history, inp_file):

    max_iter = inp_file['meta_dynamics']['max_iter']
    
    for i in range(len(results['E'])):

        if history[continue_][i][7] != [0,0]:

            results['E'][i] = 10000.0

        if history[continue_][i][0] - history[continue_][i][5] > max_iter - 1:

            results['E'][i] = 10000.0

    return results

def getting_pop_num(continue_, src_path, output_dir):
    os.chdir(src_path+"/../"+output_dir)
    f = open("totalGen","r")
    freadlines = f.readlines()
    f.close()

    linenum = len(freadlines)

    historyold = {}

    pop_num = 0

    for i in range(2,linenum):
        gen         = int(freadlines[i].split()[0])
        pop         = int(freadlines[i].split()[1])
        E           = float(freadlines[i].split()[2])
        V           = float(freadlines[i].split()[3])
        struct_type = str(freadlines[i].split()[4])
        parents     = str(freadlines[i].split()[5])
        keptbest    = str(freadlines[i].split()[6])

        if gen == continue_:
            pop_num = pop_num + 1

    return pop_num

def getting_num_of_best(results, Emin, tot_atom_num, inp_file):

    if inp_file['structure']['num_of_best'] > 0:
        num_of_best = inp_file['structure']['num_of_best']

    else:

        num_of_best = 0

        E_limit = inp_file['energy_criteria']['energy_cut_for_best_structures']*float(tot_atom_num)

        for i in range(len(results['E'])):
            if results['E'][i] - Emin < E_limit:
                num_of_best = num_of_best + 1


    return num_of_best



def write_logs(continue_, output_path, original_path):

    with open(original_path+"/totalGen","r") as f:
        totalgen = f.readlines()
    
    with open(original_path+"/totalbest","r") as f:
        totalbest = f.readlines()

    linenum = len(totalgen)
    with open(output_path+"/totalGen","a") as fw:
        for i in range(2,linenum):
            gen         = int(totalgen[i].split()[0])
            if gen <= continue_:
                fw.write(totalgen[i])

    linenum = len(totalbest)
    with open(output_path+"/totalbest","a") as fw:
        for i in range(2,linenum):
            gen         = int(totalbest[i].split()[0])
            if gen <= continue_:
                fw.write(totalbest[i])


def make_antiseed_continue(results, global_best_results, inp_file, Emin):
    
    Ecut  = inp_file['energy_criteria']['energy_cut_for_inheriting_structures']

    tot_atom_num = inp_file['tot_atom_num']
    
    ANTISEED = [10000.0 for i in range(len(results['E']))]

    nums = 0.0

    for struct1 in range(1,len(results['E'])+1):

        if results['E'][struct1-1]-Emin < Ecut*float(tot_atom_num):

              ANTISEED[struct1-1] = 0.0

    return ANTISEED
