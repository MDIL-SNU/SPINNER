from spinner.simple_nn import run
import os, sys, re
import yaml


def check_previous_status(total_yaml):
    main_dir = os.getcwd()
    md_dir = total_yaml['output_dir']['ab_initio_mqa']
    if 'anneal' not in os.listdir(os.getcwd()+'/'+md_dir):
        print("Finish ab_initio_mqa first")
        exit()

    output_dir = total_yaml['output_dir']['initial_nnp_training']
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)

    key = 1
    if 'pca' not in os.listdir() or 'scale_factor' not in os.listdir() or 'train_list' not in os.listdir() or 'valid_list' not in os.listdir():
        key = 0

    if key == 0: # if previous run is not exist
        # Make structure_list
        with open('structure_list', 'w') as s:
            s.write('[Melting]\n')
            for n in os.listdir(main_dir+'/'+md_dir+'/melt/'):
                if 'OUTCAR' in n:
                    s.write(main_dir+'/'+md_dir+'/melt/%s ::5\n'%n)
            s.write('[Quenching]\n')
            for n in os.listdir(main_dir+'/'+md_dir+'/quench/'):
                if 'OUTCAR' in n:
                    s.write(main_dir+'/'+md_dir+'/quench/%s ::5\n'%n)
            s.write('[Annealing]\n')
            for n in os.listdir(main_dir+'/'+md_dir+'/anneal/'):
                if 'OUTCAR' in n:
                    s.write(main_dir+'/'+md_dir+'/anneal/%s ::5\n'%n)

        # Make config.yaml file
        src_path = os.path.dirname(os.path.abspath(__file__))
        if 'SIMPLE-NN' in total_yaml.keys():
            initial_NNP_training_config = total_yaml['SIMPLE-NN']
        else:
            initial_NNP_training_config = dict()
        initial_NNP_training_config['params'] = dict()
        p = re.compile('[A-Z][a-z]?')
        elements = p.findall(total_yaml['composition'])
        if len(elements) == 1:
            params_file = src_path+'/params_un'
        elif len(elements) == 2:
            params_file = src_path+'/params_bi'
        elif len(elements) == 3:
            params_file = src_path+'/params_ter'
        elif len(elements) == 4:
            params_file = src_path+'/params_quater'
        elif len(elements) == 5:
            params_file = src_path+'/params_quin'

        for element in elements:
            initial_NNP_training_config['params'][element] = params_file

        with open('config.yaml', 'w') as s:
            yaml.safe_dump(initial_NNP_training_config, s, default_flow_style=False)

    elif key == 1: # If previous run exist (config.yaml)
        with open('config.yaml', 'r') as f2:
            config_yaml = yaml.safe_load(f2)
            config_yaml['generate_features'] = False
            config_yaml['preprocess'] = False
        with open('config.yaml', 'w') as s:
            yaml.safe_dump(config_yaml, s, default_flow_style=False)

    inp_file = 'config.yaml'

    return inp_file


def main():
    inp_file = sys.argv[1]
    main_dir = os.getcwd()

    # Check if yaml is serial mode
    with open(inp_file, 'r') as f:
        tmp_yaml = yaml.safe_load(f)

    if 'total_actions' in tmp_yaml: # read setting from serial mode ('total.yaml')
        if tmp_yaml['total_actions']['initial_nnp_training'] != True:
            print('Set initial_nnp_training: True in total_actions of yaml file if you want to progress')
            exit()
        inp_file = check_previous_status(tmp_yaml)

    done = 0
    if 'LOG' in os.listdir():
        with open('LOG', 'r') as f:
            if 'NNP_training Done' in f.readlines()[-1]:
                done = 1

    # SIMPLE-NN main
    if done == 0:
        run(inp_file)

    with open('LOG', 'a') as s:
        s.write('\nNNP_training Done')

if __name__ == '__main__':
    main()
