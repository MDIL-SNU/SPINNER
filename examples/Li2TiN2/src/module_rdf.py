import numpy as np
import sys
import math

def read_rdf_from_lammps(atom_nums, inp_file, rdf0):
    num_of_atom = len(atom_nums)
    sigma       = inp_file['similarity_metric']['gaussian_dist']  
    rdf_grid    = inp_file['similarity_metric']['rdf_grid']  
    num_of_atom = len(atom_nums)

    rdfs = []
    weights = []
    total_weight = 0
    
    # calculate weights
    for i in range(0,num_of_atom):
        for j in range(i,num_of_atom):
            total_weight += float(atom_nums[i]*atom_nums[j])
            weights.append(atom_nums[i]*atom_nums[j])

    for i in range(len(weights)):
        weights[i] = float(weights[i])/total_weight

    # rdf
    weight_num = -1
    rdf_num = -1

    for i in range(1,(num_of_atom+1)):
        for j in range(i,num_of_atom+1):
            weight_num += 1
            rdf = rdf0[str(i)+str(j)]

            # gaussian of rdfs
            for k in range(rdf_grid):
                gaussian_rdf = 0.0
                for l in range(int(rdf_grid*2.4)):
                    gaussian_rdf += rdf[l][1] * 1.0/sigma/np.sqrt(2.0*math.pi) * np.exp(-1.0*(rdf[l][0]-rdf[k][0]*2)**2.0/2.0/sigma**2.0) 

                gaussian_rdf += -1.0
                gaussian_rdf = gaussian_rdf * np.sqrt(weights[weight_num])

                rdfs.append(gaussian_rdf)
    return rdfs
    gc.collect()
