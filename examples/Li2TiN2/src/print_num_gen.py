import yaml
import sys
from module_input import *

input_name = sys.argv[1]

inp_file, error_message = input_yaml(input_name)

generation = inp_file['structure']['generation']
period     = inp_file['quality_monitoring']['monitoring_period']

num = int(generation/period)

with open(input_name+"_num_gen","w") as f:
    f.write("number: %d\n"%(num))
    f.write("period: %d\n"%(period)) 
