import yaml
import collections
import sys, os, re, shutil

def main():
    input_file = sys.argv[1]
    configs_path = sys.argv[2]
    if len(sys.argv) > 3:
        iteration = sys.argv[3]
    else:
        iteration = '0'

    with open(input_file, 'r') as f:
        inp_yaml = yaml.safe_load(f)

    os.makedirs(configs_path, exist_ok=True)

    make_final_spinner_config(inp_yaml, configs_path)
    nnp_dir = inp_yaml['output_dir']['initial_nnp_training']
    copy_potential(inp_yaml, nnp_dir, configs_path)
    """
    csp_iter_config = make_iterative_spinner_config(inp_yaml, configs_path, iteration)
    make_dft_oneshot_relax_config(inp_yaml, csp_iter_config, configs_path, iteration)

    if iteration == '1':
        make_final_spinner_config(inp_yaml, configs_path)
        make_evolving_nnp_feature_config(inp_yaml, configs_path)
        make_evolving_nnp_training_config(inp_yaml, configs_path)
        nnp_dir = inp_yaml['output_dir']['initial_nnp_training']
        copy_potential(inp_yaml, nnp_dir, configs_path)
    else:
        nnp_dir = inp_yaml['output_dir']['nnp_retraining']+'_%s'%(int(iteration)-1)
        copy_potential(inp_yaml, nnp_dir, configs_path)
    """

def make_iterative_spinner_config(inp_yaml, configs_path, iteration):
    p1 = re.compile('[A-Z][a-z]?')
    p2 = re.compile('[0-9][0-9]?')
    composition = inp_yaml['composition']
    elements = p1.findall(composition)
    comps = list(map(int, p2.findall(composition)))
    z, comps, npop = select_Z(comps, inp_yaml)
    md_dir = inp_yaml['output_dir']['ab_initio_mqa']
    volume = get_md_volume_per_formula(md_dir+'/melt/POSCAR')

    #if 'SPINNER_iterative' in inp_yaml.keys(): # overwirte additional key in 'input.yaml'
    #    iterative_spinner_config = inp_yaml['SPINNER_iterative']
    #else:
    iterative_spinner_config = dict()

    iterative_spinner_config['output_dir'] = inp_yaml['output_dir']['csp_iteration']+'_%s'%iteration
    iterative_spinner_config['lammps_simd'] = inp_yaml['lammps_simd']
    iterative_spinner_config['initial_volume'] = volume*z
    iterative_spinner_config['structure'] = dict()
    iterative_spinner_config['structure']['i_population'] = npop
    iterative_spinner_config['structure']['population'] = npop
    iterative_spinner_config['structure']['generation'] = 50
    iterative_spinner_config['material'] = dict()
    for i, element in enumerate(elements):
        iterative_spinner_config['material'][element] = comps[i]

    set_distance_constraint(elements, iterative_spinner_config, inp_yaml['output_dir']['ab_initio_mqa'])

    if 'SPINNER_iterative' in inp_yaml.keys():
        yaml_update(iterative_spinner_config, inp_yaml['SPINNER_iterative'])

    with open(configs_path+'/csp_iteration.yaml', 'w') as s:
        yaml.safe_dump(iterative_spinner_config, s, default_flow_style=False)

    return iterative_spinner_config

def get_encut_and_kspacing(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            en_key = 0
            lat_key = 0
            if 'ENCUT = ' in line:
                encut = float(line.split()[2])
                en_key = 1
            if 'reciprocal lattice vectors' in line:
                k_spacing = float(lines[i+6].split()[3])
                lat_key = 1
            if en_key == 1 and lat_key == 1:
                break

    return encut, k_spacing

def make_dft_oneshot_relax_config(inp_yaml, csp_iter_config, configs_path, iteration):
    dft_oneshot_relax_config = dict()
    dft_oneshot_relax_config['potential_directory'] = inp_yaml['pot_path']
    dft_oneshot_relax_config['vasp_gam'] = inp_yaml['vasp_gam']
    dft_oneshot_relax_config['vasp_std'] = inp_yaml['vasp_std']
    dft_oneshot_relax_config['vasp_ncl'] = inp_yaml['vasp_std']
    dft_oneshot_relax_config['mpi_command'] = inp_yaml['mpi_command']
    dft_oneshot_relax_config['input_dir'] = csp_iter_config['input_dir']
    dft_oneshot_relax_config['output_dir'] = csp_iter_config['output_dir']

    dft_oneshot_relax_config['num_of_samples'] = inp_yaml['DFT_oneshot_relax']['num_of_samples']
    dft_oneshot_relax_config['NPAR'] = inp_yaml['DFT_oneshot_relax']['NPAR']
    dft_oneshot_relax_config['KPAR'] = inp_yaml['DFT_oneshot_relax']['KPAR']

    md_dir = inp_yaml['output_dir']['ab_initio_mqa']
    encut, kspacing = get_encut_and_kspacing(md_dir+'/melt/OUTCAR')
    dft_oneshot_relax_config['ENCUT'] = encut
    dft_oneshot_relax_config['k-spacing'] = kspacing

    if 'DFT_oneshot_relax' in inp_yaml.keys():
        yaml_update(dft_oneshot_relax_config, inp_yaml['DFT_oneshot_relax'])

    with open(configs_path+'/dft_oneshot_relax.yaml', 'w') as s:
        yaml.safe_dump(dft_oneshot_relax_config, s, default_flow_style=False)

def make_final_spinner_config(inp_yaml, configs_path):
    p1 = re.compile('[A-Z][a-z]?')
    p2 = re.compile('[0-9][0-9]?')
    composition = inp_yaml['composition']
    elements = p1.findall(composition)
    comps = list(map(int, p2.findall(composition)))
    z, comps, npop = select_Z(comps, inp_yaml)
    md_dir = inp_yaml['output_dir']['ab_initio_mqa']
    volume = get_md_volume_per_formula(md_dir+'/melt/POSCAR')

    #if 'SPINNER_final' in inp_yaml.keys(): # overwirte additional key in 'input.yaml'
    #    final_spinner_config = inp_yaml['SPINNER_final']
    #else:
    final_spinner_config = dict()

    final_spinner_config['output_dir'] = inp_yaml['output_dir']['final_csp']
    final_spinner_config['lammps_simd'] = inp_yaml['lammps_simd']
    final_spinner_config['initial_volume'] = volume*z
    final_spinner_config['structure'] = dict()
    final_spinner_config['structure']['i_population'] = 300
    final_spinner_config['structure']['population'] = 300
    final_spinner_config['structure']['generation'] = 200
    final_spinner_config['material'] = dict()
    for i, element in enumerate(elements):
        final_spinner_config['material'][element] = comps[i]

    set_distance_constraint(elements, final_spinner_config, inp_yaml['output_dir']['ab_initio_mqa'])

    if 'SPINNER_final' in inp_yaml.keys():
        yaml_update(final_spinner_config, inp_yaml['SPINNER_final'])

    with open(configs_path+'/final_spinner.yaml', 'w') as s:
        yaml.safe_dump(final_spinner_config, s, default_flow_style=False)

def select_Z(comps, inp_yaml):
    # devide by GDC
    gdc = 1
    for ii in range(max(comps), 0, -1):
        key = 1
        for jj in range(len(comps)):
            if comps[jj] % ii != 0:
                key = 0
                break
        if key == 1:
            gdc = ii
            break
    for i in range(len(comps)):
        comps[i] = int(comps[i]/gdc)


    ntot = sum(comps)
    if ntot*4 <= inp_yaml['maximum_crystal_atoms']:
        z = 4
    elif ntot*2 <= inp_yaml['maximum_crystal_atoms']:
        z = 2
    else:
        z = 1
    comps = [i*z for i in comps]
    ntot = sum(comps)
    if ntot*2 > inp_yaml['maximum_pop']:
        npop = inp_yaml['maximum_pop']
    else:
        npop = ntot*2
    
    return z, comps, npop

def get_md_volume_per_formula(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        a=list(map(float,lines[2].split()))
        b=list(map(float,lines[3].split()))
        c=list(map(float,lines[4].split()))
        cross=[]
        cross.append(a[1]*b[2]-a[2]*b[1])
        cross.append(a[2]*b[0]-a[0]*b[2])
        cross.append(a[0]*b[1]-a[1]*b[0])
        v=cross[0]*c[0]+cross[1]*c[1]+cross[2]*c[2]

        comps = list(map(int, lines[6].split()))
        tot_atom = sum(comps)
        for i in range(1, tot_atom):
            key = 1
            for comp in comps:
                if comp % i != 0:
                    key = 0
            if key == 1:
                multiple = i

    return v/multiple

def set_distance_constraint(elements, iterative_spinner_config, md_dir):
    iterative_spinner_config['distance_constraint'] = dict()
    index = {0:'a', 1:'b', 2:'c', 3:'d', 4:'e'}
    for i, element in enumerate(elements):
        for j in range(i, len(elements)):
            with open(md_dir+'/prdf/dft_%s%s.dat'%(index[i], index[j]), 'r') as f:
                lines = f.readlines()
                max = 0
                max_d = 0
                for l, line in enumerate(lines):
                    if l < 4:
                        continue
                    d = float(line.strip().split()[1])
                    k = float(line.strip().split()[2])
                    if d > 6:
                        break
                    if k > max:
                        max = k
                        max_d = d
                min_val = max/10
                for l, line in enumerate(lines):
                    if l < 4:
                        continue
                    d = float(line.strip().split()[1])
                    k = float(line.strip().split()[2])
                    if k > min_val:
                        break
                iterative_spinner_config['distance_constraint'][elements[i]+'-'+elements[j]] = d

def yaml_update(inp0, inp):    
    for key in list(inp.keys()):
        if isinstance(inp0, collections.Mapping):
            if isinstance(inp[key], collections.Mapping) and inp[key]:
                returned = yaml_update(inp0.get(key, {}), inp[key])
                inp0[key] = returned
            else:
                inp0[key] = inp[key]
        else:
            inp0 = {key: inp[key]}

    return inp0

def make_evolving_nnp_feature_config(inp_yaml, configs_path):
    src_path = os.path.dirname(os.path.abspath(__file__))
    if 'SIMPLE-NN' in inp_yaml.keys():
        evolving_nnp_config = inp_yaml['SIMPLE-NN']

    evolving_nnp_config['generate_features'] = True
    evolving_nnp_config['preprocess'] = True
    evolving_nnp_config['train_model'] = False

    evolving_nnp_config['params'] = dict()
    p = re.compile('[A-Z][a-z]?')
    elements = p.findall(inp_yaml['composition'])
    if len(elements) == 1:
        params_file = src_path+'/nnp_train/params_un'
    elif len(elements) == 2:
        params_file = src_path+'/nnp_train/params_bi'
    elif len(elements) == 3:
        params_file = src_path+'/nnp_train/params_ter'
    elif len(elements) == 4:
        params_file = src_path+'/nnp_train/params_quater'
    elif len(elements) == 5:
        params_file = src_path+'/nnp_train/params_quin'

    for element in elements:
        evolving_nnp_config['params'][element] = params_file

    evolving_nnp_config['preprocessing']['calc_scale'] = False
    evolving_nnp_config['preprocessing']['calc_pca'] = False

    with open(configs_path+'/simplenn_feature.yaml', 'w') as s:
        yaml.safe_dump(evolving_nnp_config, s, default_flow_style=False)

def make_evolving_nnp_training_config(inp_yaml, configs_path):
    src_path = os.path.dirname(os.path.abspath(__file__))
    if 'SIMPLE-NN' in inp_yaml.keys():
        evolving_nnp_config = inp_yaml['SIMPLE-NN']

    evolving_nnp_config['generate_features'] = False
    evolving_nnp_config['preprocess'] = False
    evolving_nnp_config['train_model'] = True

    evolving_nnp_config['params'] = dict()
    p = re.compile('[A-Z][a-z]?')
    elements = p.findall(inp_yaml['composition'])
    if len(elements) == 1:
        params_file = src_path+'/nnp_train/params_un'
    elif len(elements) == 2:
        params_file = src_path+'/nnp_train/params_bi'
    elif len(elements) == 3:
        params_file = src_path+'/nnp_train/params_ter'
    elif len(elements) == 4:
        params_file = src_path+'/nnp_train/params_quater'
    elif len(elements) == 5:
        params_file = src_path+'/nnp_train/params_quin'

    for element in elements:
        evolving_nnp_config['params'][element] = params_file

    evolving_nnp_config['neural_network']['continue'] = 'weights'
    evolving_nnp_config['neural_network']['total_epoch'] = 500
    evolving_nnp_config['neural_network']['learning_rate'] = 1e-05

    with open(configs_path+'/simplenn_training.yaml', 'w') as s:
        yaml.safe_dump(evolving_nnp_config, s, default_flow_style=False)

def copy_potential(inp_yaml, nnp_dir, config_paths):
    #with open(config_paths+'/csp_iteration.yaml', 'r') as f:
    with open(config_paths+'/final_spinner.yaml', 'r') as f:
        tmp_yaml = yaml.safe_load(f)
        pot_path = tmp_yaml['input_dir']
    os.makedirs(pot_path, exist_ok=True)

    last_idx = 0
    for n in os.listdir(nnp_dir):
        if 'potential_saved_epoch' in n:
            idx = int(n.split('_')[-1])
            if idx > last_idx:
                last_idx = idx
    shutil.copy(nnp_dir+'/potential_saved_epoch_%s'%last_idx, pot_path+'/potential')



if __name__ == "__main__":
    main()
