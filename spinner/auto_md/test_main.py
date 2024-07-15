from module_input import initialize_auto_md
from module_vasp import check_vasp_version
from module_premelt import run_premelting
from module_conv_test import run_convergence_test
from module_relaxation import run_initial_relaxation
from module_Tm_prediction import run_melting_temperature_prediction
from module_MQA import run_melting, run_quenching, run_annealing
from module_util import move_to_directory, write_log, calculate_total_elapsed_time

# initialize
input_yaml = initialize_auto_md()
is_action_needed = input_yaml['Actions']
log = input_yaml['log_path']
TEST_MODE = True

# MD start
if is_action_needed['premelt']:
    run_premelting(input_yaml)

if is_action_needed['convergence_test']:
    run_convergence_test(input_yaml)

# Check if KP is gamma or not
vasp_ver = check_vasp_version(log)
input_yaml['vasp_version'] = vasp_ver

if is_action_needed['relax_md']:
    run_initial_relaxation(input_yaml, TEST_MODE)
    
if is_action_needed['find_Tm']:
    run_melting_temperature_prediction(input_yaml, TEST_MODE)

if is_action_needed['melt']:
    run_melting(input_yaml)

if is_action_needed['quench']:
    run_quenching(input_yaml, TEST_MODE)

if is_action_needed['anneal']:
    run_annealing(input_yaml)

write_log(log, 'Calculation done')
write_log(log, f'Total time: {calculate_total_elapsed_time(log)}')
move_to_directory(input_yaml['working_dir'])
