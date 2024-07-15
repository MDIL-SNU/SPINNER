import argparse
from spinner.auto_md.module_input import initialize_auto_md
from spinner.auto_md.module_vasp import check_vasp_version
from spinner.auto_md.module_premelt import run_premelting
from spinner.auto_md.module_conv_test import run_convergence_test
from spinner.auto_md.module_relaxation import run_initial_relaxation
from spinner.auto_md.module_Tm_prediction import run_melting_temperature_prediction
from spinner.auto_md.module_MQA import run_melting, run_quenching, run_annealing
from spinner.auto_md.module_util import move_to_directory, write_log, calculate_total_elapsed_time
from spinner.auto_md.module_prdf import calculate_prdf
import os

description = (
    #f'spinner version={SPINNER_VERSION}, automated molecular dynamics are conducted based on input.yaml'
)

input_yaml_help = 'Input.yaml for running MD'
working_dir_help = 'Path to write output. Default is Auto_MD'
num_process_help = 'Number of process cores'
mode_help = 'Select mode. Set [auto_md] if run only auto_md. Set [full] if run full spinner process.'

def main(args=None):
    md_setting_yaml, num_tasks, working_dir, mode = cmd_parse_main(args)

    # initialize
    input_yaml = initialize_auto_md(md_setting_yaml, num_tasks, working_dir)
    is_action_needed = input_yaml['Actions']
    log = input_yaml['log_path']

    # MD start
    if is_action_needed['premelt']:
        run_premelting(input_yaml)

    if is_action_needed['convergence_test']:
        run_convergence_test(input_yaml)

    # Check if KP is gamma or not
    vasp_ver = check_vasp_version(log)
    input_yaml['vasp_version'] = vasp_ver

    if is_action_needed['relax_md']:
        run_initial_relaxation(input_yaml)

    if is_action_needed['find_Tm']:
        run_melting_temperature_prediction(input_yaml)

    if is_action_needed['melt']:
        run_melting(input_yaml)

    if is_action_needed['quench']:
        run_quenching(input_yaml)

    if is_action_needed['anneal']:
        run_annealing(input_yaml)

    if is_action_needed['prdf']:
        calculate_prdf(input_yaml, target_dir='melt/')

    write_log(log, 'Calculation done')
    write_log(log, f'Total time: {calculate_total_elapsed_time(log)}')
    move_to_directory(input_yaml['working_dir'])

def cmd_parse_main(args=None):
    ag = argparse.ArgumentParser(description=description)
    ag.add_argument('input_yaml', help=input_yaml_help, type=str)
    ag.add_argument(
        '-np',
        '--num_process',
        help=num_process_help,
        type=int,
    )
    ag.add_argument(
        '-w',
        '--working_dir',
        help=working_dir_help,
        type=str,
        default='Auto_MD',
    )
    ag.add_argument(
        '-m',
        '--mode',
        help=mode_help,
        type=str,
        default='auto_md',
    )

    args = ag.parse_args()
    input_yaml = args.input_yaml
    num_tasks = args.num_process
    wd = args.working_dir
    mode = args.mode

    return input_yaml, num_tasks, wd, mode

if __name__ == '__main__':
    main()
