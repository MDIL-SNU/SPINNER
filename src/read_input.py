import yaml
import sys
from module_input import *

input_name = sys.argv[1]


inp_file, error_message = input_yaml(input_name)

with open(input_name+"_read","w") as f:
    yaml.dump(inp_file,f)


