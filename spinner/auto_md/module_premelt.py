import os
import shutil as shu
from spinner.auto_md.rand_structure import create_rand_structure
from spinner.auto_md.module_util import write_log, write_log_with_timestamp, calculate_elapsed_time, create_and_move_to_directory, move_to_directory, get_abs_path
from spinner.auto_md.module_vasp import edit_INCAR, run_vasp

def run_premelting(input_yaml):
    working_dir = input_yaml['working_dir']
    log = input_yaml['log_path']
    start_time_premelt = write_log_with_timestamp(log, "Premelting start")
    create_and_move_to_directory('premelt')
    do_premelting(input_yaml)
    move_to_directory(working_dir)
    end_time = write_log_with_timestamp(log, "Premelting is done")
    write_log(log, f"Premelting time: {calculate_elapsed_time(start_time_premelt, end_time)}\n")

def do_premelting(input_yaml):
    material = input_yaml['composition']
    working_dir = input_yaml['working_dir']
    potential_dir = input_yaml['pot_dir']
    target_num = input_yaml['target_num']
    src_dir = input_yaml['src_dir']
    nsw = input_yaml['premelt_config']['steps']
    vasp_config = input_yaml['vasp_config']
    
    # random spray
    create_rand_structure(material, potential_dir, target_num) 
    shu.copy('POTCAR', get_abs_path(working_dir, 'Inputs/POTCAR'))
    shu.copy(get_abs_path(src_dir, 'INCAR_premelt'), './INCAR')
    edit_INCAR('./INCAR', {'SYSTEM':material, 'NSW': nsw})
    shu.copy(get_abs_path(src_dir, 'KPOINTS'), './KPOINTS')
    run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config['gam'])
    

