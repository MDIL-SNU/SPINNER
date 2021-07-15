import yaml

def input_yaml(input_file):
    # read input
    with open(input_file,'r') as inp:
        inp_yaml = yaml.safe_load(inp)
    
    # check if there is missing input
    error_message = check_if_there_is_missing_input(inp_yaml)

    if error_message == "":
    
      # calculate total_atom_num
      tot_atom_num = 0
      for atom in inp_yaml['material']:
        tot_atom_num = tot_atom_num + inp_yaml['material'][atom]
      inp_yaml['tot_atom_num'] = int(tot_atom_num) 

      # default input
      inp_yaml = default_inputs(inp_yaml)

    return inp_yaml, error_message


def check_if_there_is_missing_input(inp_yaml):

    # requirements are input_dir, output_dir, structure, material, initial volume

    if   'output_dir' not in inp_yaml:
        return "output_dir: is missing in input.yaml"

    elif 'input_dir'  not in inp_yaml:
        return "input_dir: is missing in input.yaml"

    elif 'structure' not in inp_yaml:
        return "structure: generation: is missing in input.yaml"

    elif 'generation' not in inp_yaml['structure']:
        return "structure: generation: is missing in input.yaml"

    elif 'material' not in inp_yaml:
        return "material: generation: is missing in input.yaml"

    elif 'initial_volume' not in inp_yaml:
        return "initial_volume: is missing in input.yaml"

    else:
        return ""


def default_inputs(inp_yaml):
    
    # operator
    if 'operator' not in inp_yaml:
        inp_yaml['operator'] = {}
        inp_yaml['operator']['random'] = 0.7
        inp_yaml['operator']['crossover'] = 0.0
        inp_yaml['operator']['con_permutation'] = 0.0
        inp_yaml['operator']['all_permutation'] = 0.2
        inp_yaml['operator']['latticemutation'] = 0.1

    else:
        if 'random' not in inp_yaml['operator']:
            inp_yaml['operator']['random'] = 0.0
        if 'crossover' not in inp_yaml['operator']:
            inp_yaml['operator']['crossover'] = 0.0
        if 'con_permutation' not in inp_yaml['operator']:
            inp_yaml['operator']['con_permutation'] = 0.0
        if 'all_permutation' not in inp_yaml['operator']:
            inp_yaml['operator']['all_permutation'] = 0.0
        if 'latticemutation' not in inp_yaml['operator']:
            inp_yaml['operator']['latticemutation'] = 0.0

    # structure
    if 'i_population' not in inp_yaml['structure']:
        inp_yaml['structure']['i_population'] = inp_yaml['tot_atom_num']*2

    if 'population' not in inp_yaml['structure']:
        inp_yaml['structure']['population'] = inp_yaml['tot_atom_num']*2

    if 'num_of_best' not in inp_yaml['structure']:
        inp_yaml['structure']['num_of_best'] = 0

    if 're-relax_best' not in inp_yaml['structure']:
        inp_yaml['structure']['re-relax_best'] = True

    # energy_critetria:
    if 'energy_criteria' not in inp_yaml:
        inp_yaml['energy_criteria'] = {}

    if 'energy_cut_for_inheriting_structures' not in inp_yaml['energy_criteria']:
        inp_yaml['energy_criteria']['energy_cut_for_inheriting_structures'] = 0.10

    if 'energy_cut_for_best_structures' not in inp_yaml['energy_criteria']:
        inp_yaml['energy_criteria']['energy_cut_for_best_structures'] = 0.05

    if 'energy_cut_for_further_relax' not in inp_yaml['energy_criteria']:
        inp_yaml['energy_criteria']['energy_cut_for_further_relax'] = 100.0

    # crossover condition
    if 'crossover_condition' not in inp_yaml:
        inp_yaml['crossover_condition'] = {}

    if 'num_of_grid_for_cut' not in inp_yaml['crossover_condition']:
        inp_yaml['crossover_condition']['num_of_grid_for_cut'] = 10

    if 'energy_range_for_cut_select' not in inp_yaml['crossover_condition']:
        inp_yaml['crossover_condition']['energy_range_for_cut_select'] = 10

    if 'grid_for_shift' not in inp_yaml['crossover_condition']:
        inp_yaml['crossover_condition']['grid_for_shift'] = 3

    if 'iteration_for_add_atoms' not in inp_yaml['crossover_condition']:
        inp_yaml['crossover_condition']['iteration_for_add_atoms'] = 50

    # random_condition
    if 'random_condition' not in inp_yaml:
        inp_yaml['random_condition'] = {}

    if 'force_general_Wyckoff_site' not in inp_yaml['random_condition']:
        inp_yaml['random_condition']['force_general_Wyckoff_site'] = False
        
    if 'maximum_attempts_for_one_space_group_and_volume' not in inp_yaml['random_condition']:
        inp_yaml['random_condition']['maximum_attempts_for_one_space_group_and_volume'] = 100

    if 'scale_factor' not in inp_yaml['random_condition']:
        inp_yaml['random_condition']['scale_factor'] = 1.0

    if 'sublattice_generation' not in inp_yaml['random_condition']:
        inp_yaml['random_condition']['sublattice_generation'] = 0.0

    if 'max_sub_lattice' not in inp_yaml['random_condition']:
        inp_yaml['random_condition']['max_sub_lattice'] = 2

    # relax condition
    if 'relax_condition' not in inp_yaml:
        inp_yaml['relax_condition'] = {}

    if 'relax_iteration' not in inp_yaml['relax_condition']:
        inp_yaml['relax_condition']['relax_iteration'] = 5
    if 'method_of_first_relax' not in inp_yaml['relax_condition']:
        inp_yaml['relax_condition']['method_of_first_relax'] = 'cg'
    if 'further_calculate_with_accurate_potential' not in inp_yaml['relax_condition']:
        inp_yaml['relax_condition']['further_calculate_with_accurate_potential'] = False

    # distance constraint
    if 'distance_constraint' not in inp_yaml:
        inp_yaml['distance_constraint'] = {}

    for atom1 in inp_yaml['material']:
        for atom2 in inp_yaml['material']:

            if (atom1+'-'+atom2 not in inp_yaml['distance_constraint']) and (atom2+'-'+atom1 not in inp_yaml['distance_constraint']):
                inp_yaml['material'][atom1+'-'+atom2] = 0.7

    # vacuum constraint
    if 'vacuum_constraint' not in inp_yaml:
        inp_yaml['vacuum_constraint'] = {}

    if 'apply_vacuum_constraint' not in inp_yaml['vacuum_constraint']:
        inp_yaml['vacuum_constraint']['apply_vacuum_constraint'] = True

    if 'maximum_vacuum_length' not in inp_yaml['vacuum_constraint']:
        inp_yaml['vacuum_constraint']['maximum_vacuum_length'] = 10.0

    if 'grid' not in inp_yaml['vacuum_constraint']:
        inp_yaml['vacuum_constraint']['grid'] = 1.0

    # similarity metric
    if 'similarity_metric' not in inp_yaml:
        inp_yaml['similarity_metric'] = {}

    if 'type' not in inp_yaml['similarity_metric']:
        inp_yaml['similarity_metric']['type'] = 'pRDF'

    if 'limit' not in inp_yaml['similarity_metric']:
        inp_yaml['similarity_metric']['limit'] = min(0.08*40/float(tot_num_atom))

    if 'volume_cut' not in inp_yaml['similarity_metric']:
        inp_yaml['similarity_metric']['volume_cut'] = 0.1

    if 'energy_cut' not in inp_yaml['similarity_metric']:
        inp_yaml['similarity_metric']['energy_cut'] = 0.005

    if 'gaussian_dist' not in inp_yaml['similarity_metric']:
        inp_yaml['similarity_metric']['gaussian_dist'] = 0.1

    if 'rdf_grid' not in inp_yaml['similarity_metric']:
        inp_yaml['similarity_metric']['rdf_grid'] = 250

    # antiseed
    if 'antiseed' not in inp_yaml:
        inp_yaml['antiseed'] = {}

    if 'activate_antiseed' not in inp_yaml['antiseed']:
        inp_yaml['antiseed']['activate_antiseed'] = False

    if 'gaussian_width' not in inp_yaml['antiseed']:
        inp_yaml['antiseed']['gaussian_width'] = 0.2

    if 'selection_gaussian' not in inp_yaml['antiseed']:
        inp_yaml['antiseed']['selection_gaussian'] = 0.1

    if 'selection_fraction' not in inp_yaml['antiseed']:
        inp_yaml['antiseed']['selection_fraction'] = 0.5

    # continue
    if 'continue' not in inp_yaml:
        inp_yaml['continue'] = {}
        inp_yaml['continue']['continue_num'] = 0

    return inp_yaml


def make_tolerance_matrix(atomnamelist, inp_file):
    
    tolerance_matrix = [[0.0 for i in range(len(atomnamelist))] for j in range(len(atomnamelist))]

    for atom1num in range(len(atomnamelist)):
        for atom2num in range(len(atomnamelist)): 
            atom1 = atomnamelist[atom1num]
            atom2 = atomnamelist[atom2num]
            pair1 = atom1 + "-" + atom2
            pair2 = atom2 + "-" + atom1

            if pair1 in inp_file['distance_constraint']:
                tolerance_matrix[atom1num][atom2num] = inp_file['distance_constraint'][pair1]
                tolerance_matrix[atom2num][atom1num] = inp_file['distance_constraint'][pair1]
            elif pair2 in inp_file['distance_constraint']: 
                tolerance_matrix[atom2num][atom1num] = inp_file['distance_constraint'][pair2]
                tolerance_matrix[atom1num][atom2num] = inp_file['distance_constraint'][pair2]
                
    return tolerance_matrix    
