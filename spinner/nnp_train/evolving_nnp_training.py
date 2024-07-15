from simple_nn import run
import os, sys, shutil
import yaml
from check_convergence_NNP import check_converge


main_dir = os.getcwd()
if 'STATUS' in os.listdir('.'):
    inp_file = main_dir+'/'+sys.argv[1]
    with open(sys.argv[2], 'r') as f:
        total_yaml = yaml.safe_load(f)
    current_iteration = int(sys.argv[3])

    retraining_dir = total_yaml['output_dir']['nnp_retraining']+'_%s'%current_iteration
    os.makedirs(retraining_dir, exist_ok=True)
    os.chdir(retraining_dir)
else:
    inp_file = 'input.yaml'

run(inp_file)
#for i in range(1, 51):
#    run(inp_file)
#    check_converge(total_yaml)
#    os.rename('LOG', 'LOG%s'%i)
#    if 'nnp_done' in os.listdir():
#        break

######### Additional code ############
os.chdir(main_dir)
if 'STATUS' in os.listdir(main_dir):
    last_idx = 0
    for n in os.listdir(retraining_dir):
        if 'potential_saved_epoch' in n:
            idx = int(n.split('_')[-1])
            if idx > last_idx:
                last_idx = idx

    #pot_path = total_yaml['SPINNER_iterative']['input_dir']
    #shutil.copy(retraining_dir+'/potential_saved_epoch_%s'%last_idx, main_dir+'/'+pot_path+'/potential')
    if last_idx != 0:
        with open(main_dir+'/STATUS', 'w') as s:
            s.write('NNP_training%s Done'%(current_iteration+1))
####################################
