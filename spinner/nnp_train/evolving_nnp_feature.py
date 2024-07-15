from simple_nn import run
import os, sys, shutil
import yaml


main_dir = os.getcwd()
if 'STATUS' in os.listdir('.'):
    inp_file = main_dir+'/'+sys.argv[1]
    with open(sys.argv[2], 'r') as f:
        total_yaml = yaml.safe_load(f)
    current_iteration = int(sys.argv[3])

    initial_nnp_dir = total_yaml['output_dir']['initial_NNP_training']
    csp_dir = total_yaml['output_dir']['csp_iteration']
    retraining_dir = total_yaml['output_dir']['nnp_retraining']
    os.makedirs(retraining_dir+'_%s'%current_iteration, exist_ok=True)
    os.chdir(retraining_dir+'_%s'%current_iteration)

    shutil.copy(main_dir+'/'+csp_dir+'_%s/str_list0'%(current_iteration), 'structure_list')
else:
    inp_file = 'input.yaml'

run(inp_file)

######### Additional code ############
if 'STATUS' in os.listdir(main_dir):
    os.chdir(main_dir)
    shutil.copy(initial_nnp_dir+'/pca', retraining_dir+'_%s/pca'%current_iteration)
    shutil.copy(initial_nnp_dir+'/scale_factor', retraining_dir+'_%s/scale_factor'%current_iteration)
    os.rename(retraining_dir+'_%s/train_list'%current_iteration, retraining_dir+'_%s/train_list_tmp'%current_iteration)
    os.rename(retraining_dir+'_%s/valid_list'%current_iteration, retraining_dir+'_%s/valid_list_tmp'%current_iteration)
    os.rename(retraining_dir+'_%s/LOG'%current_iteration, retraining_dir+'_%s/LOG_feature'%current_iteration)


    with open(retraining_dir+'_%s/train_list'%current_iteration, 'w') as s:
        with open(initial_nnp_dir+'/train_list', 'r') as f:
            lines = f.readlines()
            for line in lines:
                s.write(line)
        for i in range(1, current_iteration+1):
            with open(retraining_dir+'_%s/train_list_tmp'%i, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    s.write(line)
    with open(retraining_dir+'_%s/valid_list'%current_iteration, 'w') as s:
        with open(initial_nnp_dir+'/valid_list', 'r') as f:
            lines = f.readlines()
            for line in lines:
                s.write(line)
        for i in range(1, current_iteration+1):
            with open(retraining_dir+'_%s/valid_list_tmp'%i, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    s.write(line)

    prev_pot = total_yaml['SPINNER_iterative']['input_dir']
    shutil.copy(main_dir+'/'+prev_pot+'/potential', retraining_dir+'_%s/potential_saved'%current_iteration)

####################################
