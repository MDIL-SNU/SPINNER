import ase.io
import numpy as np
import os, re
from lammps import lammps
from spinner.auto_md.module_util import write_log, write_log_with_timestamp, calculate_elapsed_time


def calculate_prdf(input_yaml, target_dir='melt/'):
    working_dir = input_yaml['working_dir']
    log = input_yaml['log_path']
    melt_steps = input_yaml['melt_config']['steps']
    composition = input_yaml['composition']
    p = re.compile('[A-Z][a-z]?')
    elements = p.findall(composition)
    results_dir = 'prdf/'
    start_time = write_log_with_timestamp(log, "PRDF start")

    os.makedirs(results_dir, exist_ok=True)
    convert_poscar_to_lammps(target_dir, results_dir)
    tot_steps = convert_xdatcar_to_lammps(target_dir, results_dir, melt_steps)
    lammps_prdf_calculation(input_yaml, elements, results_dir, tot_steps)

    end_time = write_log_with_timestamp(log, "PRDF is done")
    write_log(log, f"PRDF time: {calculate_elapsed_time(start_time, end_time)}\n")

def convert_poscar_to_lammps(target_dir, results_dir):
    data = ase.io.read(target_dir+'/POSCAR', format='vasp')
    ase.io.write(results_dir+'/coo', data, format='lammps-data', force_skew=True)

def convert_xdatcar_to_lammps(target_dir, results_dir, melt_steps):
    formatter = {'float_kind':lambda x: "%.9f" %x}
    if melt_steps >= 200:
        interval_step = int(melt_steps/200)
        start_step = melt_steps - 200*interval_step
        tot_steps = 200
    else:
        interval_step = 1
        start_step = 0
        tot_steps = melt_steps

    with open(results_dir+'/xdatcar.dump', 'w') as fp:
        snap = ase.io.read(target_dir+'/XDATCAR', '%s::%s'%(start_step, interval_step), format='vasp-xdatcar')

        elems = list()
        for i in snap[0].numbers:
            if i not in elems:
                elems.append(i)

        timestep = 0
        for atoms in snap:
            #for i in range(len(elems)):
                #atoms.numbers[atoms.numbers == elems[i]] = i+1
            for i in range(len(atoms)):
                atoms.numbers[i] = elems.index(atoms.numbers[i])+1

            cell = atoms.get_cell_lengths_and_angles()
            lx = cell[0]
            xy = cell[1]*np.cos(cell[5]*np.pi/180)
            xz = cell[2]*np.cos(cell[4]*np.pi/180)
            ly = np.sqrt(cell[1]**2 - xy**2)
            yz = (cell[1]*cell[2]*np.cos(cell[3]*np.pi/180) - xy*xz)/ly
            lz = np.sqrt(cell[2]**2 - xz**2 - yz**2)
        
            atoms.set_cell(np.array([[lx, 0.0, 0.0],[xy, ly, 0.0], [xz, yz, lz]]), scale_atoms=True)
            fp.write("ITEM: TIMESTEP\n") 
            fp.write("{}\n".format(timestep)) 
            fp.write("ITEM: NUMBER OF ATOMS\n")
            fp.write("{}\n".format(len(atoms.get_chemical_symbols())))
            fp.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
            fp.write("{} {} {}\n".format(min([0.0, xy, xz, xy+xz]),lx + max([0.0, xy, xz, xy+xz]), xy))
            fp.write("{} {} {}\n".format(min([0.0, yz]), ly + max(0.0, yz), xz))
            fp.write("0.0 {} {}\n".format(lz, yz))
            fp.write("ITEM: ATOMS id type x y z\n")
            
            for idx, (types, pos) in enumerate(zip(atoms.numbers, atoms.get_positions())):
                pos = np.array2string(pos, formatter=formatter).strip('[]')
                fp.write("{} {} {}\n".format((idx+1), types, pos))
            
            timestep += 1

    return tot_steps

def lammps_prdf_calculation(input_yaml, elements, results_dir, tot_steps):
    os.chdir(results_dir)

    if input_yaml['lammps_simd'] == True:
        lmp = lammps('simd_serial')
    else:
        lmp = lammps()
    lmp.command("boundary p p p")
    lmp.command("processors * * * grid numa")
    lmp.command("units metal")
    lmp.command("read_data coo")
    lmp.command("pair_style zero 10.0")
    lmp.command("pair_coeff * *")
    for i, element in enumerate(elements):
        lmp.command("mass %s 1"%(i+1))
    lmp.command("thermo 1")
    lmp.command("thermo_style custom step")
    lmp.command("compute rdf_tot all rdf 1000 * * cutoff 10.0")
    index = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e'}
    for i, element in enumerate(elements):
        for j in range(i, len(elements)):
            lmp.command("compute rdf_%s%s all rdf 1000 %s %s cutoff 10.0"%(index[i], index[j], i+1, j+1))
    lmp.command("fix print1 all ave/time 1 %s %s c_rdf_tot[*] file dft_tot.dat mode vector"%(tot_steps-1, tot_steps-1))
    idx = 2
    for i, element in enumerate(elements):
        for j in range(i, len(elements)):
            lmp.command("fix print%s all ave/time 1 %s %s c_rdf_%s%s[*] file dft_%s%s.dat mode vector"%(idx, tot_steps-1, tot_steps-1, index[i], index[j], index[i], index[j]))
            idx += 1
    lmp.command("rerun ./xdatcar.dump dump x y z purge yes add yes box yes replace no")

    os.chdir('..')

