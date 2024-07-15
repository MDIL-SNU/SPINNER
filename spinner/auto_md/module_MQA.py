import os, re
from spinner.auto_md.module_util import write_log_with_timestamp, write_log, create_and_move_to_directory, calculate_elapsed_time, move_to_directory, get_abs_path
from spinner.auto_md.module_vasp import copy_inputs, edit_INCAR, run_vasp, grep_nth_item

def run_melting(input_yaml):
    working_dir = input_yaml['working_dir']
    log = input_yaml['log_path']
    start_time = write_log_with_timestamp(log, "Melting start")
    create_and_move_to_directory('melt')
    do_melting(input_yaml)
    move_to_directory(working_dir)
    end_time = write_log_with_timestamp(log, "Melting is done")
    write_log(log, f"Melting time: {calculate_elapsed_time(start_time, end_time)}\n")

def run_quenching(input_yaml, test_mode=False):
    working_dir = input_yaml['working_dir']
    log = input_yaml['log_path']
    start_time = write_log_with_timestamp(log, "Quenching start")
    create_and_move_to_directory('quench')
    do_quenching(input_yaml, test_mode)
    move_to_directory(working_dir)
    end_time = write_log_with_timestamp(log, "Quenching is done")
    write_log(log, f"Quenching time: {calculate_elapsed_time(start_time, end_time)}\n")

def run_annealing(input_yaml):
    working_dir = input_yaml['working_dir']
    log = input_yaml['log_path']
    start_time = write_log_with_timestamp(log, "Annealing start")
    create_and_move_to_directory('anneal')
    do_annealing(input_yaml)
    move_to_directory(working_dir)
    end_time = write_log_with_timestamp(log, "Annealing is done")
    write_log(log, f"Annealing time: {calculate_elapsed_time(start_time, end_time)}\n")

def check_last():
    p = re.compile('\d+')
    step = 0
    last_idx = 0
    for n in os.listdir('.'):
        if 'OUTCAR' in n:
            if n != 'OUTCAR':
                idx = int(p.findall(n)[0])
                if idx > last_idx:
                    last_idx = idx
            with open(n, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'free  ' in line:
                        step += 1
    if 'OUTCAR' in os.listdir('.'):
        last_idx += 1
        os.rename('OUTCAR', 'OUTCAR%s'%last_idx)
        os.rename('XDATCAR', 'XDATCAR%s'%last_idx)
        os.rename('POSCAR', 'POSCAR%s'%last_idx)
        os.rename('CONTCAR', 'POSCAR')

    return step, last_idx

def do_melting(input_yaml):
    working_dir = input_yaml['working_dir']
    vasp_ver = input_yaml['vasp_version']
    step, _ = check_last()
    if step == 0:
        copy_inputs(get_abs_path(working_dir, 'Inputs'), './', ['POSCAR_to_melt','POTCAR','KPOINTS','INCAR'])
        os.rename('POSCAR_to_melt', 'POSCAR')
        edit_INCAR('./INCAR', {'NSW': input_yaml['melt_config']['steps']})
    else:
        new_nsw = int(input_yaml['melt_config']['steps']) - step
        edit_INCAR('./INCAR', {'NSW': f'{new_nsw}'})

    run_vasp(input_yaml['vasp_config']['mpicommand'], input_yaml['vasp_config']['num_tasks'], input_yaml['vasp_config'][vasp_ver])

def do_quenching(input_yaml, test_mode=False):
    working_dir = input_yaml['working_dir']
    vasp_ver = input_yaml['vasp_version']
    step, _ = check_last()
    if step == 0:
        copy_inputs(get_abs_path(working_dir, 'melt'), './', ['CONTCAR','POTCAR','KPOINTS','INCAR'])
        os.rename('CONTCAR', 'POSCAR')
        nsw = int((float(grep_nth_item('TEBEG', './INCAR',2)[0]) - input_yaml['quench_config']['Temp_end'])/input_yaml['quench_config']['quenching_rate']*1000/2)
        edit_INCAR('./INCAR', {'NSW':f'{nsw}', 'TEEND': input_yaml['quench_config']['Temp_end']})
    else:
        nsw=int((float(grep_nth_item('TEBEG', './INCAR',2)[0]) - input_yaml['quench_config']['Temp_end'])/input_yaml['quench_config']['quenching_rate']*1000/2)
        new_nsw = nsw - step
        new_tebeg = 300 + new_nsw/nsw*int((float(grep_nth_item('TEBEG', get_abs_path(working_dir, 'melt/INCAR'), 2)[0]) - input_yaml['quench_config']['Temp_end']))
        edit_INCAR('./INCAR', {'NSW': f'{new_nsw}', 'TEBEG': f'{new_tebeg}'})

    if test_mode:
        edit_INCAR('./INCAR', {'NSW': '30'})
    run_vasp(input_yaml['vasp_config']['mpicommand'], input_yaml['vasp_config']['num_tasks'], input_yaml['vasp_config'][vasp_ver])

def do_annealing(input_yaml):
    working_dir = input_yaml['working_dir']
    vasp_ver = input_yaml['vasp_version']
    step, _ = check_last()
    if step == 0:
        copy_inputs(get_abs_path(working_dir, 'quench'), './', ['CONTCAR','POTCAR','KPOINTS','INCAR'])
        os.rename('CONTCAR', 'POSCAR')
        edit_INCAR('./INCAR', {'NSW': input_yaml['annealing_config']['steps'], 'TEBEG': input_yaml['annealing_config']['Temp_start'],
                'TEEND': input_yaml['annealing_config']['Temp_end']})
    else:
        new_nsw = int(input_yaml['annealing_config']['steps']) - step
        edit_INCAR('./INCAR', {'NSW': f'{new_nsw}'})

    run_vasp(input_yaml['vasp_config']['mpicommand'], input_yaml['vasp_config']['num_tasks'], input_yaml['vasp_config'][vasp_ver])
