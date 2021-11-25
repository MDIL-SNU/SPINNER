import os
import sys
import shutil
from module_input import *
import random
import subprocess

input_name = sys.argv[1]
iteration  = sys.argv[2]
inp_file, error_message = input_yaml(input_name)
folder= inp_file['output_dir']
num_sample   = int(inp_file['sampling_trainingset']['num_of_samples'])

os.chdir(folder)
absolute_path = os.getcwd()

os.mkdir('relax'+iteration)

os.chdir('dft')

a = os.listdir("Done")

os.chdir('Done')
        

fstr = open("../../str_list"+iteration,"a")

if len(a) != 0:
    fstr.write("\n[final"+iteration+"-relax: 20.0]\n")

for i in range(1,num_sample+1):
    for j in range(1,10000):
        if os.path.isfile(str(i)+'/relax_GGA/OUTCAR'+str(j)):
            shutil.copy(str(i)+'/relax_GGA/OUTCAR'+str(j),'../../relax'+iteration+'/OUTCAR'+str(i)+'-'+str(j))
            fstr.write(absolute_path+"/relax"+iteration+"/OUTCAR"+str(i)+'-'+str(j)+"  ::100000\n")
        else:
            break


fstr.close()
