import sys
import os
import shutil
import yaml
from module_input import *


# input
work_directory = sys.argv[1]
output_dir        = sys.argv[2]
iteration         = int(sys.argv[3])
input_name        = sys.argv[4]

inp_file, error_message = input_yaml(work_directory+"/"+input_name)

# read simplenn.yaml
simplennyaml = inp_file['NNP_setting']['input_file']
with open(work_directory+"/"+simplennyaml,'r') as inp:
    simplenn_file = yaml.safe_load(inp)

with open("../input.yaml",'r') as inp:
    original_file = yaml.safe_load(inp)

# write new input

simplenn_file['symmetry_function']['params'] = original_file['symmetry_function']['params']
simplenn_file['symmetry_function']['scale_rho'] = {}
simplenn_file['atom_types'] = []

for atom in inp_file['material']:
  simplenn_file['atom_types'].append(atom)
  simplenn_file['symmetry_function']['scale_rho'][atom] = 0.03

simplenn_file['generate_features'] = False
simplenn_file['preprocess']  = False
simplenn_file['train_model'] = True

simplenn_file['neural_network']['pca'] = True
simplenn_file['neural_network']['pca_whiten'] = True
simplenn_file['symmetry_function']['scale_type'] = False

simplenn_file['neural_network']['total_iteration'] = 2000

# make inputs
with open('input.yaml','w') as f:
    yaml.dump(simplenn_file,f)

# write train, valid list
with open("../train_list","r") as f:
    trainlist = f.readlines() 
with open("../valid_list","r") as f:
    validlist = f.readlines() 


lista = os.listdir('data')

with open("train_list","w") as f:
    for flist in lista:
        if 'training' in flist:
            f.write('data/'+flist)
            f.write("\n")
    for fstr in trainlist:
        f.write("../"+fstr)

with open("valid_list","w") as f:
    for flist in lista:
        if 'valid' in flist:
            f.write('data/'+flist)
            f.write("\n")
    for fstr in validlist:
        f.write("../"+fstr)
