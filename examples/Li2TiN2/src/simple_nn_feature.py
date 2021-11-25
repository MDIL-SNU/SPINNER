import sys
import os
import shutil
import yaml
from module_input import *


# input
work_directory    = sys.argv[1]
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

simplenn_file['generate_features'] = True
simplenn_file['preprocess']  = True
simplenn_file['train_model'] = False

simplenn_file['neural_network']['pca'] = False
simplenn_file['neural_network']['pca_whiten'] = False
simplenn_file['symmetry_function']['scale_type'] = False

# make inputs
with open('input.yaml','w') as f:
  yaml.dump(simplenn_file,f)
