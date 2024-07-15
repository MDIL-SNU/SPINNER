import sys, yaml, collections, os
from spinner.auto_md.module_util import create_and_move_to_directory

def initialize_auto_md(md_setting_yaml, num_tasks, working_dir):
    src_dir = os.path.dirname(os.path.abspath(__file__))
    with open(md_setting_yaml, 'r') as O:
        tmp_yaml = yaml.safe_load(O)

    if 'total_actions' in tmp_yaml: # read setting from sequential mode ('total.yaml')
        if tmp_yaml['total_actions']['ab_initio_mqa'] != True:
            print('Set ab_initio_mqa: True in total_actions of yaml file if you want to progress')
            exit()
        input_yaml = tmp_yaml['ab_initio_mqa']
        input_yaml['pot_dir'] = tmp_yaml['pot_path']
        input_yaml['composition'] = tmp_yaml['composition']
        input_yaml['vasp_config']['std'] = tmp_yaml['vasp_std']
        input_yaml['vasp_config']['gam'] = tmp_yaml['vasp_gam']
        input_yaml['vasp_config']['ncl'] = tmp_yaml['vasp_ncl']
        input_yaml['vasp_config']['mpicommand'] = tmp_yaml['mpi_command']
        input_yaml['mode'] = 'serial'
        input_yaml['lammps_simd'] = tmp_yaml['lammps_simd']
        input_yaml['Actions']['prdf'] = True
        input_yaml['working_dir'] = os.getcwd()+'/'+tmp_yaml['output_dir']['ab_initio_mqa']
    else:  # read setting from single mode ('auto_mqa.yaml')
        input_yaml = tmp_yaml
        input_yaml['mode'] = 'single'
        input_yaml['lammps_simd'] = False
        input_yaml['working_dir'] = os.getcwd()+'/'+working_dir

    input_yaml = load_input(src_dir, input_yaml)
    input_yaml['vasp_config']['num_tasks'] = str(num_tasks)

    # Setting directories
    create_and_move_to_directory(input_yaml['working_dir'])

    write_input_to_yaml(input_yaml)
    os.makedirs("Inputs", exist_ok=True)

    log = check_log(input_yaml)
    input_yaml['src_dir'] = src_dir
    input_yaml['log_path'] = log
    return input_yaml


def update_input(inp0, inp):
    for key in list(inp.keys()):
        if isinstance(inp0, collections.Mapping):
            if isinstance(inp[key], collections.Mapping) and inp[key]:
                returned = update_input(inp0.get(key, {}), inp[key])
                inp0[key] = returned
            else:
                inp0[key] = inp[key]
        else:
            inp0 = {key: inp[key]}

    return inp0

def convert(inp0, key):
    check = inp0[key].upper()
    if check[0] == '.':
        if check[1] == 'T':
            check = True
        elif check[1] == 'F':
            check = False
    elif check[0] == 'T':
        check = True
    elif check[0] == 'F':
        check = False

    if isinstance(check, bool):
        inp0[key] = check

def load_input(inp0_dir, inp):
    # Open default yaml and update yaml
    with open(inp0_dir+'/configure_default.yaml', 'r') as O:
        inp0=yaml.safe_load(O)
    #with open(inp_dir, 'r') as O:
    #    inp=yaml.safe_load(O)
    update_input(inp0, inp)

    # Convert boolean strings to boolean (ex. 'T' -> True, '.F.' -> False)
    #for key in inp0.keys():
        #if not isinstance(inp0[key], bool) and isinstance(inp0[key], str):
            #convert(inp0, key)

    return inp0

def write_input_to_yaml(inp_yaml):
    with open('./final_configure.yaml', 'w') as O:
        yaml.safe_dump(inp_yaml, O, default_flow_style=False)

def check_log(inp_yaml):
    if 'LOG' in os.listdir('./'):
        with open('LOG', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'Premelting is done' in line:
                    inp_yaml['Actions']['premelt'] = False
                elif 'Convergence test is done' in line:
                    inp_yaml['Actions']['convergence_test'] = False
                elif 'Relax-MD is done' in line:
                    inp_yaml['Actions']['relax_md'] = False
                elif 'Predicting melting temperature is done' in line:
                    inp_yaml['Actions']['find_Tm'] = False
                elif 'Melting is done' in line:
                    inp_yaml['Actions']['melt'] = False
                elif 'Quenching is done' in line:
                    inp_yaml['Actions']['quench'] = False
                elif 'Annealing is done' in line:
                    inp_yaml['Actions']['anneal'] = False
                elif 'PRDF is done' in line:
                    inp_yaml['Actions']['prdf'] = False
    else:
        with open('LOG', 'a') as s:
            s.write('Calculation start\n')

    log = os.getcwd()+'/LOG' 

    return log
