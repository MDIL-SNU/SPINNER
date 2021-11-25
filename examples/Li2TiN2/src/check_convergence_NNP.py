import sys
import os
import shutil
import yaml
from module_input import *


# input
input_name        = sys.argv[1]
inp_file, error_message = input_yaml(input_name)

E_limit = inp_file['NNP_setting']['training_energy_limit']
F_limit = inp_file['NNP_setting']['training_force_limit']
S_limit = inp_file['NNP_setting']['training_stress_limit']

# investigate
with open("LOG","r") as f:
    freadlines = f.readlines()

linelength = len(freadlines)

for i in range(linelength-1,0,-1):
    fstr = freadlines[i]
    if 'RMSE' in fstr and 'iteration' in fstr:
        startline = i
        break

converged = 1

MQ_E_error = []
MQ_F_error = []
MQ_S_error = []

for i in range(startline+4,linelength-3):
    fstr = freadlines[i]
    
    if 'iteration' not in fstr:
        MQ_E_error.append(float(fstr.split()[2]))
        MQ_F_error.append(float(fstr.split()[4]))
        MQ_S_error.append(float(fstr.split()[6]))

    elif 'high' not in fstr and 'low' not in fstr:
        E_error = float(fstr.split()[2])
        F_error = float(fstr.split()[4])
        S_error = float(fstr.split()[6])

        if E_error > E_limit or F_error > F_limit or S_error > S_limit:
            converged = 0


if min(MQ_E_error) > E_limit or min(MQ_F_error) > F_limit or min(MQ_S_error) > S_limit:
    converged = 0

if converged == 1:
    with open("nnp_done","w") as f:
        f.write("nnp_done")
