import shutil as shu
import os
from spinner.auto_md.module_util import write_log, write_log_with_timestamp, calculate_elapsed_time, create_and_move_to_directory, move_to_directory, get_abs_path
from spinner.auto_md.module_vasp import edit_INCAR, copy_inputs, run_vasp, scale_velo, grep_nth_item
from spinner.auto_md.module_relaxation import do_final_relaxation_and_md
import numpy as np

def run_melting_temperature_prediction(input_yaml, test_mode=False):
    working_dir = input_yaml['working_dir']
    log = input_yaml['log_path']
    start_time = write_log_with_timestamp(log, "Predicting melting temperature start")
    try:
        move_to_directory('find_Tm')
    except:
        pass
    continue_dir = get_abs_path(os.getcwd(), 'Step0_relax', str(len(os.listdir('Step0_relax'))-1))
    input_yaml['continue_dir'] = continue_dir
    melting_temperature = find_melting_temperature(input_yaml)
    do_final_relaxation_and_md(input_yaml, melting_temperature, test_mode)
    move_to_directory(working_dir)
    end_time = write_log_with_timestamp(log, "Predicting melting temperature is done")
    write_log(log, f"Predicting melting temperature time: {calculate_elapsed_time(start_time, end_time)}\n")


def find_melting_temperature(input_yaml):
    working_dir = input_yaml['working_dir']
    continue_dir = input_yaml['continue_dir']
    input_Tm = input_yaml['find_Tm_config']
    vasp_config = input_yaml['vasp_config'] 
    vasp_version = input_yaml['vasp_version']
    log = input_yaml['log_path']
    find_tm_dir = os.getcwd()

    Tm, is_Tm_found = check_Tm_prediction_checkpoint(log)
    
    if Tm != 0:
        continue_dir = get_abs_path(find_tm_dir, f'/Step1_T_{Tm}')
        start_temp = Tm - input_Tm['StepSize_T']
    else:
        start_temp = input_Tm['Start_T']

    ###find MSD at each T
    if not is_Tm_found:
        for temperature in range(start_temp, input_Tm['End_T']-input_Tm['StepSize_T'], -input_Tm['StepSize_T']):
            create_and_move_to_directory(f'Step1_T_{temperature}')
            copy_inputs(continue_dir, './', ['INCAR', 'KPOINTS', 'POTCAR', 'CONTCAR'])
            shu.move('./CONTCAR', './POSCAR')
            scale_velo('./POSCAR', temperature)
            edit_INCAR('./INCAR', {'NSW': str(input_Tm['duration']),
                        'IBRION': '0', 'MDALGO': '2', 'SMASS': '0', 'TEBEG': str(temperature),
                        'TEEND': str(temperature), 'ISIF': '2', 'POTIM': '2'})
            run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config[vasp_version])
            if os.path.isfile('./XDATCAR_unwrapped.xyz'):
                sim_info = read_xdat_unwrap('./XDATCAR_unwrapped.xyz')
            else:
                sim_info = save_xdat_unwrap('./XDATCAR')

            msd = calculate_MSD(sim_info)

            if msd[-1,-1]/msd[-1,0]/6/2*1e-8 >= 4e-9:
                continue_dir = os.getcwd()
                shu.copy('./INCAR', get_abs_path(working_dir, 'Inputs/INCAR'))
                Tm = temperature
                move_to_directory(find_tm_dir)
                write_log(log, f"{temperature} K melting is done")
            else:
                move_to_directory(find_tm_dir)
                write_log(log, f"{temperature} K melting is done\nPredicted Tm: {Tm}")
                is_Tm_found = True
                break
        else:
            write_log(log, f"Tm is less than {Tm} K\nPredicted Tm: {Tm}")
            
    if is_Tm_found and Tm == 0:
        Tm = 4000
        create_and_move_to_directory(f'Step1_T_{Tm}')
        copy_inputs(continue_dir, './', ['INCAR', 'KPOINTS', 'POTCAR', 'CONTCAR'])
        shu.move('./CONTCAR', './POSCAR')
        scale_velo('./POSCAR', Tm)
        edit_INCAR('./INCAR', {'NSW': str(input_Tm['duration']),
                    'IBRION': '0', 'MDALGO': '2', 'SMASS': '0', 'TEBEG': '4000',
                    'TEEND': '4000', 'ISIF': '2', 'POTIM': '2'})
        run_vasp(vasp_config['mpicommand'], vasp_config['num_tasks'], vasp_config[vasp_version])
        shu.copy('./INCAR', get_abs_path(working_dir, 'Inputs/INCAR'))
        move_to_directory(find_tm_dir)
        write_log(log, f"{Tm} K melting is done\nPredicted Tm: {Tm}")

    return Tm


def check_Tm_prediction_checkpoint(log):
    Tm = 0
    is_Tm_found = False
    with open(log, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'K melting is done' in line:
                Tm = int(line.split()[0])
            elif 'Predicted Tm:' in line:
                Tm = int(line.split()[2])
                is_Tm_found = True
            
    return Tm, is_Tm_found


def save_xdat_unwrap(xdat_dir):
    vec=np.zeros((3,3))
    with open(xdat_dir,'r') as O:
        xdat=O.readlines()
    Atoms=[int(i) for i in xdat[6].split()]
    num_of_ions=sum(Atoms)
    for i in range(3):
        vec[i,:]=np.array([float(j) for j in xdat[i+2].split()])
    tot_len=(len(xdat)-7)//(num_of_ions+1)

    tot_coord=np.zeros((tot_len,num_of_ions,3))
    rel_to_move=np.zeros((tot_len,num_of_ions,3))

    for step in range(0,tot_len):
        for atom_coord in range(0,num_of_ions):
            tot_coord[step,atom_coord,:]= \
            [float(coord) for coord in xdat[step*(num_of_ions+1)+8+atom_coord].split()]
    for i in range(1,tot_len-1):
        rel_to_move[i:,:,:]+=np.round(tot_coord[i,:,:]-tot_coord[i-1,:,:])
    tot_coord+=rel_to_move
    tot_coord=np.matmul(tot_coord.reshape(-1,3),vec).reshape(tot_len,-1,3)
    with open(xdat_dir+'_unwrapped.xyz','w') as O:
        for i in range(7):
            O.write(xdat[i])
        for i in range(tot_len):
            O.write(f'Cartesian configuration= {i+1}\n')
            for j in range(num_of_ions):
                O.write(f'{tot_coord[i,j,0]} {tot_coord[i,j,1]} {tot_coord[i,j,2]}\n')

    return vec,Atoms,tot_coord


def read_xdat_unwrap(xdat_dir):
    vec=np.zeros((3,3))
    with open(xdat_dir,'r') as O:
        xdat=O.readlines()
    Atoms=[int(i) for i in xdat[6].split()]
    NION=sum(Atoms)
    for i in range(3):
        vec[i,:]=np.array([float(j) for j in xdat[i+2].split()])
    tot_len=(len(xdat)-7)//(NION+1)

    tot_coord=np.zeros((tot_len,NION,3))

    for step in range(0,tot_len):
        for atom_coord in range(0,NION):
            tot_coord[step,atom_coord,:]= \
            [float(coord) for coord in xdat[step*(NION+1)+8+atom_coord].split()]

    return vec,Atoms,tot_coord


def calculate_MSD(sim_info):
    vec=sim_info[0]
    Atoms=sim_info[1]
    Nspecies=len(Atoms)
    NION=sum(Atoms)
    tot_coord=sim_info[2]
    tot_len=tot_coord.shape[0]
    ### TODO: sim_time should be input & msd_len should be configurable
    sim_time=2e-3
    msd_len=1500
    msd_interval=5

    if tot_len>=msd_len+msd_interval:
        MSD=np.zeros((msd_len,Nspecies+2))
        MSD[:,0]=np.arange(0,msd_len)*sim_time

        ensemble_N=(tot_len-msd_len)//msd_interval
        for ensemble in range(ensemble_N):
            initial=tot_coord[ensemble*msd_interval,:,:].copy()
            for step in range(1,msd_len):
                start_count=0
                for e_at,at_num in enumerate(Atoms):
                    MSD[step,e_at+1]+=np.sum(np.power(tot_coord[ensemble*msd_interval+step,start_count:start_count+at_num,:]-initial[start_count:start_count+at_num,:],2))/at_num
                    start_count+=at_num
        for e_at,at_num in enumerate(Atoms):
            MSD[:,-1]+=MSD[:,e_at+1]*at_num
        MSD[:,-1]/=NION
        MSD[:,1:]/=ensemble_N

    else:
        MSD=np.zeros((tot_len,Nspecies+2))
        MSD[:,0]=np.arange(0,tot_len)*sim_time

        initial=tot_coord[0,:,:].copy()
        for step in range(1,tot_len):
            start_count=0
            for e_at,at_num in enumerate(Atoms):
                MSD[step,e_at+1]+=np.sum(np.power(tot_coord[step,start_count:start_count+at_num,:]-initial[start_count:start_count+at_num,:],2))/at_num
                start_count+=at_num
                for e_at,at_num in enumerate(Atoms):
                    MSD[step,-1]+=MSD[step,e_at+1]*at_num
                MSD[step,-1]/=NION
    writer=''
    for j in MSD:
        for k in j:
            writer+=f'{k} '
        writer+='\n'
    with open('XDATCAR.msd','w') as O:
        O.write('#Time ')
        for i in range(Nspecies):
            O.write(f'Atom{i+1} ')
        O.write('Average\n')
        O.write(writer)
    return MSD

def pearson(MSD):
    msd=MSD[:,-1]

    xy_dat=np.zeros((len(msd)+1,2),dtype=float)
    xy_dat[:,0]=np.arange(0,len(msd)+1).reshape(1,-1)
    for e,i in enumerate(msd):
        xy_dat[e+1,1]=i
    #print(xy_dat)
    return  np.sum((xy_dat[:,0]-np.mean(xy_dat[:,0]))*(xy_dat[:,1]-np.mean(xy_dat[:,1]))) /np.sqrt(np.sum(  np.power(xy_dat[:,0]-np.mean(xy_dat[:,0]),2.)  )) /np.sqrt(np.sum(  np.power(xy_dat[:,1]-np.mean(xy_dat[:,1]),2.)  ))


if __name__ == '__main__':
    #continue_dir = '/data/haekwan98/Auto-MD/md_test/Output/Na1Cl1/find_Tm/Step1_T_1000'
    #pressures = np.array([float(p_line) for p_line in grep_nth_item('external', continue_dir+'/OUTCAR',3)[:]])
    #print(pressures)
    pass
