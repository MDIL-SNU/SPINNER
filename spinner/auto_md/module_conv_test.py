import os, time
import shutil as shu
import numpy as np
from spinner.auto_md.module_util import write_log, write_log_with_timestamp, calculate_elapsed_time, create_and_move_to_directory, move_to_directory, get_abs_path, get_cpu_core_number, get_divisor
from spinner.auto_md.module_vasp import edit_INCAR, edit_KP, run_vasp, copy_inputs, get_energy_force_stress, get_maximum_enmax_value, get_loop_plus_real_time

def run_convergence_test(input_yaml):
    working_dir = input_yaml['working_dir']
    log = input_yaml['log_path']
    start_time = write_log_with_timestamp(log, "Convergence test start")
    create_and_move_to_directory('conv_test')
    do_convergence_test(input_yaml)
    move_to_directory(working_dir)
    end_time = write_log_with_timestamp(log, "Convergence test is done")
    write_log(log, f"Convergence test time: {calculate_elapsed_time(start_time, end_time)}\n")


def do_convergence_test(input_yaml):
    is_kptest = input_yaml['conv_test_config']['kp_use'].upper() == 'AUTO'
    is_cutoff = input_yaml['conv_test_config']['cutoff_use'].upper() == 'AUTO'
    working_dir = input_yaml['working_dir']
    input_yaml['conv_test_config']['reference_cutoff'] = set_reference_cutoff_energy(get_abs_path(working_dir, 'Inputs/POTCAR'))

    reference_energy_force_stress = calculate_reference(input_yaml)
    chosen_kpoint = do_kpoint_test(input_yaml, reference_energy_force_stress)
    do_cutoff_test(input_yaml, reference_energy_force_stress, chosen_kpoint)
    #TODO NPAR TEST
    do_npar_test(input_yaml, chosen_kpoint)

def calculate_reference(input_yaml):
    working_dir = input_yaml['working_dir']
    conv_test_dir = os.getcwd()
    conv_test_config = input_yaml['conv_test_config']
    vasp_config = input_yaml['vasp_config']
    log = input_yaml['log_path']

    create_and_move_to_directory('reference')
    # Calculation of reference
    copy_inputs(from_dir=get_abs_path(working_dir, 'premelt'), to_dir='./', copy_list=['INCAR', 'KPOINTS', 'POTCAR', 'CONTCAR'])
    shu.move('CONTCAR', 'POSCAR')
    edit_INCAR(target='./INCAR', to_change={'ENCUT': conv_test_config['reference_cutoff'], 'PREC': 'normal', 'NSW': '0', 
                                            'POTIM': '#', 'TEBEG': '#', 'TEEND': '#', 'SMASS': '#', 'MDALGO': '#'})
    if vasp_config['npar']//2 >= 1 and vasp_config['npar']%2 == 0:
        edit_INCAR(target='./INCAR', to_change={'KPAR': '2', 'NPAR': str(vasp_config['npar']//2)})
    edit_KP(target='./KPOINTS', num=conv_test_config['reference_KP'])
    write_log(log, f"Reference calculation \t ENCUT: {conv_test_config['reference_cutoff']} | KP: {conv_test_config['reference_KP']}")
    run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['std'])
    num_of_ions = get_num_of_ions('./POSCAR')
    reference_energy_force_stress = get_energy_force_stress('./OUTCAR', num_of_ions)
    energy, force, stress = reference_energy_force_stress
    write_log(log, f"Reference calculation is done \t E: {energy/num_of_ions} meV/atom | S: {stress} kbar")
    move_to_directory(conv_test_dir)
    return reference_energy_force_stress

def get_num_of_ions(poscar_path):
    with open(poscar_path, 'r') as O:
        ions = O.readlines()[6].split()
    num_of_ions = sum([int(n_tmp) for n_tmp in ions])
    return num_of_ions


def do_kpoint_test(input_yaml, reference_energy_force_stress):
    conv_test_dir = os.getcwd()
    working_dir = input_yaml['working_dir']
    conv_test_config = input_yaml['conv_test_config']
    vasp_config = input_yaml['vasp_config']
    log = input_yaml['log_path']

    create_and_move_to_directory('kptest')
    kptest_dir = os.getcwd()
    for current_KP in ['G', 'B']:
        create_and_move_to_directory(current_KP)
        copy_inputs(from_dir=get_abs_path(conv_test_dir, 'reference'), to_dir='./', copy_list=['INCAR','KPOINTS','POTCAR','POSCAR'])
        edit_INCAR(target='./INCAR', to_change={'KPAR': '#', 'NPAR': str(vasp_config['npar'])})
        edit_KP(target='./KPOINTS', special_kp=current_KP)
        if current_KP=='G':
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['gam'])
        else:
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['std'])
        num_of_ions = get_num_of_ions('./POSCAR')
        calculated_energy_force_stress = get_energy_force_stress('./OUTCAR', num_of_ions)
        is_converged, delta_E, delta_F, delta_S = check_convergence(reference_energy_force_stress, calculated_energy_force_stress, num_of_ions, conv_test_config)
        if is_converged:
            shu.copy('KPOINTS', get_abs_path(working_dir, 'Inputs'))
            shu.copy('INCAR', get_abs_path(working_dir, 'Inputs'))
            write_log(log, f"Converged KP: {current_KP} \t delta E: {delta_E} meV/atom | max delta F: {delta_F} eV/A | delta S: {delta_S} kbar")
            break
        else:
            write_log(log, f"{current_KP} not converged \t delta E: {delta_E} meV/atom | max delta F: {delta_F} eV/A | delta S: {delta_S} kbar")
        move_to_directory(kptest_dir)
    else:
        shu.copy('B/KPOINTS', get_abs_path(working_dir, 'Inputs'))
        shu.copy('B/INCAR', get_abs_path(working_dir, 'Inputs'))
        write_log(log, "KP is not converged. Use B")
            
    move_to_directory(conv_test_dir)
    return current_KP


def do_cutoff_test(input_yaml, reference_energy_force_stress, chosen_kpoint):
    conv_test_dir = os.getcwd()
    working_dir = input_yaml['working_dir']
    conv_test_config = input_yaml['conv_test_config']
    vasp_config = input_yaml['vasp_config']
    log = input_yaml['log_path']

    create_and_move_to_directory('cutoff')
    cutoff_dir = os.getcwd()
    min_encut = 150
    encut_start = conv_test_config['reference_cutoff']
    encut_stepsize = conv_test_config['cutoff_stepsize']
    e_range =  range(encut_start-encut_stepsize, min_encut-1, -encut_stepsize)
    for current_cutoff in e_range:
        create_and_move_to_directory(str(current_cutoff))
        copy_inputs(from_dir=get_abs_path(conv_test_dir, 'reference'), to_dir='./', copy_list=['POTCAR','POSCAR'])
        copy_inputs(from_dir=get_abs_path(working_dir, 'Inputs'), to_dir='./', copy_list=['INCAR','KPOINTS'])
        edit_INCAR('./INCAR', {'ENCUT':str(current_cutoff)})
        if chosen_kpoint == 'G':
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['gam'])
        else:
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['std'])
        num_of_ions = get_num_of_ions('./POSCAR')
        calculated_energy_force_stress = get_energy_force_stress('./OUTCAR', num_of_ions)
        is_converged, delta_E, delta_F, delta_S = check_convergence(reference_energy_force_stress, calculated_energy_force_stress, num_of_ions, conv_test_config)
        if is_converged:
            shu.copy('./INCAR', get_abs_path(working_dir, 'Inputs/INCAR'))
            write_log(log, f"ENCUT {current_cutoff} done \t delta E: {delta_E} meV/atom | max delta F: {delta_F} eV/A | delta S: {delta_S} kbar")
            move_to_directory(cutoff_dir)
        else:
            write_log(log, f"ENCUT {current_cutoff} done \t delta E: {delta_E} meV/atom | max delta F: {delta_F} eV/A | delta S: {delta_S} kbar")
            write_log(log, f"Converged ENCUT: {current_cutoff+conv_test_config['cutoff_stepsize']}")
            break
    else:
        write_log(log, f"Converged ENCUT: {current_cutoff}")
    
    move_to_directory(conv_test_dir)

def do_npar_test(input_yaml, chosen_kpoint):
    working_dir = input_yaml['working_dir']
    conv_test_dir = os.getcwd()
    vasp_config = input_yaml['vasp_config']
    log = input_yaml['log_path']

    create_and_move_to_directory('npar')
    npar_dir = os.getcwd()

    # TODO: NPAR TEST
    npar_list = get_divisor(get_cpu_core_number())
    npar_times = []
    for npar in npar_list:
        create_and_move_to_directory(str(npar))
        # copy files
        copy_inputs(from_dir=get_abs_path(working_dir, 'Inputs'), to_dir='./', copy_list=['INCAR', 'KPOINTS', 'POTCAR'])
        copy_inputs(from_dir=get_abs_path(conv_test_dir, 'reference'), to_dir='./', copy_list=['POSCAR'])
        # edit INCAR
        edit_INCAR(target='./INCAR', to_change={'NPAR': str(npar)})
        # calculation
        if chosen_kpoint == 'G':
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['gam'])
        else:
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['std'])
        calculation_time = get_loop_plus_real_time('./OUTCAR')
        npar_times.append((npar, calculation_time))
        write_log(log, f"NPAR {npar} done \t time: {calculation_time} s")
        move_to_directory(npar_dir)

    # find NPAR with minimum time & copy INCAR
    min_npar, _ = min(npar_times, key=lambda x: x[1])
    write_log(log, f"Optimal NPAR: {min_npar}")
    shu.copy(get_abs_path(npar_dir, f'{str(min_npar)}/INCAR'), get_abs_path(working_dir, 'Inputs/INCAR'))
    move_to_directory(conv_test_dir)
        
def check_convergence(EFS_ref, EFS_check, NION, criteria):
    delta_E = abs(EFS_ref[0]-EFS_check[0])/NION
    delta_F = np.sqrt(np.max(np.sum(np.power(EFS_ref[1]-EFS_check[1], 2), axis=1)))
    delta_S = abs(EFS_ref[2]-EFS_check[2])
    if (delta_E <= criteria['E_tol'] and 
       delta_F <= criteria['F_tol'] and 
       delta_S <= criteria['S_tol']):
        return 1, delta_E, delta_F, delta_S
    else:
        return 0, delta_E, delta_F, delta_S
    
def set_reference_cutoff_energy(potcar_path):
    maximum_enmax = get_maximum_enmax_value(potcar_path)
    return maximum_enmax + 100




