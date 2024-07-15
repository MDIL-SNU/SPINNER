import os
import shutil as shu
import numpy as np
from spinner.auto_md.module_vasp import edit_INCAR, run_vasp, get_velocity_and_write, grep_nth_item, copy_inputs
from spinner.auto_md.module_util import create_and_move_to_directory, write_log, write_log_with_timestamp, calculate_elapsed_time, move_to_directory, get_abs_path

def run_initial_relaxation(input_yaml, test_mode=False):
    log = input_yaml['log_path']
    start_time = write_log_with_timestamp(log, "Relax-MD start")
    create_and_move_to_directory('find_Tm')
    do_initial_relaxation(input_yaml, test_mode)
    end_time = write_log_with_timestamp(log, "Relax-MD is done")
    write_log(log, f"Relax&Md time: {calculate_elapsed_time(start_time, end_time)}\n")


def do_initial_relaxation(input_yaml, test_mode):
    find_tm_dir = os.getcwd()
    working_dir = input_yaml['working_dir']
    vasp_config = input_yaml['vasp_config']
    vasp_version = input_yaml['vasp_version']
    log = input_yaml['log_path']

    iter_to_start, done = check_relaxation_checkpoint(log)
    if done:
        return
    
    create_and_move_to_directory('Step0_relax')
    relaxation_dir = os.getcwd()
    initial_config = set_initial_config(working_dir, iter_to_start)

    # Main relaxation loop
    file_for_velocity = initial_config['POS_dir']
    for i in range(iter_to_start, 5):
        # relax
        create_and_move_to_directory(str(i))
        i_th_dir = os.getcwd()
        shu.copy(initial_config['INCAR_dir'], './INCAR')
        shu.copy(initial_config['KP_dir'], './KPOINTS')
        shu.copy(initial_config['POT_dir'], './POTCAR')
        shu.copy(initial_config['POS_dir'], './POSCAR')
        edit_INCAR('./INCAR', {'IBRION': '2', 'NSW': '10', 'ISIF': '7', 'SMASS': '#', 'MDALGO': '#', 
                               'POTIM': '#', 'TEBEG': '#', 'TEEND': '#', 'PSTRESS': '30'})
        run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config[vasp_version])

        # md
        shu.move('./CONTCAR', './POSCAR')
        shu.move('./OUTCAR', './OUTCAR_rlx')
        get_velocity_and_write(file_for_velocity, './POSCAR')
        edit_INCAR('./INCAR', {'IBRION': '0', 'NSW': '500', 'ISIF': '2', 'SMASS': '0', 'MDALGO': '2', 
                               'POTIM': '2', 'TEBEG': '4500', 'TEEND': '4500', 'PSTRESS': ''})
        if test_mode:
            edit_INCAR('./INCAR', {'NSW': '50'})
        run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config[vasp_version])
        
        ##get pressure after 0.5 ps to equilibriate
        pressures = np.array([float(p_line) for p_line in grep_nth_item('external', './OUTCAR', 3)[250:]]) 
        if np.average(pressures) <= 50 and np.average(pressures) >= -10:
            write_log(log, f"{i} iter relaxation - average pressure: {np.average(pressures)} kbar\nPressure is equilibrated")
            break
        else:
            file_for_velocity = get_abs_path(i_th_dir, 'CONTCAR')
            initial_config['INCAR_dir'] = get_abs_path(i_th_dir, 'INCAR')
            initial_config['KP_dir'] = get_abs_path(i_th_dir, 'KPOINTS')
            initial_config['POT_dir'] = get_abs_path(i_th_dir, 'POTCAR')
            initial_config['POS_dir'] = get_abs_path(i_th_dir, 'CONTCAR')
            move_to_directory(relaxation_dir)
            write_log(log, f"{i} iter relaxation - average pressure: {np.average(pressures)} kbar")
            if test_mode:
                break
            
            ### TODO: if not converged?
    move_to_directory(find_tm_dir)


def do_final_relaxation_and_md(input_yaml, Tm, test_mode):
    input_Tm = input_yaml['find_Tm_config']
    vasp_config = input_yaml['vasp_config']
    vasp_version = input_yaml['vasp_version']
    log = input_yaml['log_path']
    working_dir = input_yaml['working_dir']
    find_tm_dir = os.getcwd()

    # Final volume relax
    continue_dir = get_abs_path(find_tm_dir, f'Step1_T_{Tm}')
    pressures = np.array([float(p_line) for p_line in grep_nth_item('external', get_abs_path(continue_dir, 'OUTCAR'), 3)[:]])
    if  np.abs(np.mean(pressures)) >= 30:
        ### post-relax after finding Tm if needed
        velocity_from_poscar = get_abs_path(continue_dir, 'CONTCAR')
        create_and_move_to_directory(f'Step2_relax')
        relaxation_dir = os.getcwd()
        for i in range(input_Tm['maximum_post_rlx']):
            create_and_move_to_directory(str(i))
            # relax
            copy_inputs(continue_dir, './', ['INCAR','KPOINTS','POTCAR','CONTCAR'])
            shu.move('./CONTCAR', './POSCAR')
            edit_INCAR('./INCAR', {'IBRION': '2', 'NSW': '10', 'ISIF': '7', 'SMASS': '#',
                                   'MDALGO': '#', 'POTIM': '#', 'TEBEG': '#', 'TEEND': '#'})
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config[vasp_version])
            
            # md
            shu.move('./CONTCAR', './POSCAR')
            shu.move('./OUTCAR', './OUTCAR_rlx')
            get_velocity_and_write(velocity_from_poscar, './POSCAR')
            edit_INCAR('./INCAR', {'IBRION': '0', 'NSW': '500', 'ISIF': '2', 'MDALGO': '2',
                                   'SMASS': '0', 'TEBEG': f'{Tm}', 'TEEND': f'{Tm}', 'POTIM': '2'})
            if test_mode:
                edit_INCAR('./INCAR', {'NSW': '50'})
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config[vasp_version])

            continue_dir = os.getcwd()
            pressures = np.array([float(p_line) for p_line in grep_nth_item('external', './OUTCAR', 3)[250:]])
            if  np.abs(np.mean(pressures)) <= 30:
                shu.copy('./CONTCAR', get_abs_path(working_dir, 'Inputs/POSCAR_to_melt'))
                write_log(log, f"{i} iter final relaxation is done\nPressure is equilibrated")
                break
            else:
                velocity_from_poscar = get_abs_path(continue_dir, 'CONTCAR')
                shu.copy('./CONTCAR', get_abs_path(working_dir, 'Inputs/POSCAR_to_melt'))
                move_to_directory(relaxation_dir)
                write_log(log, f"{i} iter final relaxation is done")
            #TODO if converged pressure, break
        move_to_directory(find_tm_dir)
    else:
        shu.copy(f'Step1_T_{Tm}/CONTCAR', get_abs_path(working_dir, 'Inputs/POSCAR_to_melt'))
        write_log(log, "Pressure is equilibrated")


def set_initial_config(working_dir, iter_to_start):
    relaxation_dir = os.getcwd()
    initial_config={'POS_dir': get_abs_path(working_dir, 'premelt/CONTCAR'), 
                'INCAR_dir': get_abs_path(working_dir, 'Inputs/INCAR'),
                'POT_dir': get_abs_path(working_dir, 'premelt/POTCAR'), 
                'KP_dir': get_abs_path(working_dir, 'Inputs/KPOINTS')}
    if iter_to_start != 0:
        move_to_directory(str(iter_to_start-1))
        current_dir = os.getcwd()
        initial_config['INCAR_dir'] = get_abs_path(current_dir, 'INCAR')
        initial_config['KP_dir'] = get_abs_path(current_dir, 'KPOINTS')
        initial_config['POT_dir'] = get_abs_path(current_dir, 'POTCAR')
        initial_config['POS_dir'] = get_abs_path(current_dir, 'CONTCAR')
        move_to_directory(relaxation_dir)
    return initial_config


def check_relaxation_checkpoint(log):
    done = False
    iter_to_start = 0
    with open(log, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "Pressure is equilibrated" in line:
                done = True
            elif "iter initial relaxation is done" in line:
                iter_to_start = int(line.split()[0])+1
    return iter_to_start, done
