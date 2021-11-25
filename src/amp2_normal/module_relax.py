###########################################
### Date: 2018-12-05                        ###
### yybbyb@snu.ac.kr                        ###
###########################################
# This is a package of modules for structure optimization.
from module_log import *
from module_vasprun import *
from module_amp2_input import read_poscar
from module_vector import *

# set nsw for relax
def set_nsw(poscar,incar):
        pos = open(poscar,'r').readlines()[6].split()
        natom = 0
        for i in range(len(pos)) :
                natom = natom + int(pos[i])
#        nsw = natom*3
#        if nsw < 20:
#                nsw = 20
        nsw = 100
        nsw = str(nsw)
        wincar(incar,incar,[['NSW',nsw]],[])
        return int(nsw)

# Lattice mismatch check between relaxed poscar and original poscar
def set_pos_compare(poscar_rlx,poscar_ori,target,err_percent):
        axis_rlx = poscar_to_axis(poscar_rlx)
        axis_ori = poscar_to_axis(poscar_ori)

        min_diff = 0.01

        warning = 0
        for i in range(3):
                for j in range(3):
                        if abs(axis_rlx[i][j]-axis_ori[i][j]) > min_diff and abs((axis_rlx[i][j]-axis_ori[i][j])/(axis_ori[i][j]+0.001)) > float(err_percent/100):
                                warning = 1
        if warning == 1:
                make_amp2_log(target,'Warning!! More than One of lattice parameter changes more than 1 %.')
        return warning

def spacing_check(outcar_initial,kp_initial):
        import math
        initial_tmp = []
        tagtag = False
        tt = 0
        with open(outcar_initial,'r') as out:
                for line in out :
                        if "length of vectors" in line:
                                tagtag = True
                        if tagtag and tt ==0:
                                tt += 1
                        elif tagtag and tt !=0:
                                ll = line.split()
                                ll = [float(k) for k in ll]
                                initial_tmp.append(ll)
                                tt += 1
                                tagtag = False
                                break
        with open(kp_initial,'r') as kp:
                kp.readline()        
                kp.readline()        
                kp.readline()
                ini_kp = kp.readline().split()
        ini_kp = [float(k) for k in ini_kp]
        ini_spacing1 = initial_tmp[0][3]/ini_kp[0]
        ini_spacing2 = initial_tmp[0][4]/ini_kp[1]
        ini_spacing3 = initial_tmp[0][5]/ini_kp[2]
        return [ini_spacing1,ini_spacing2,ini_spacing3]


def reciprocal_compare(outcar_now,kp_now,a_spacing,b_spacing,c_spacing):
        import math
        recalc = False
        with open(kp_now,'r') as kp:
                kp.readline()        
                kp.readline()        
                kp.readline()
                now_kp = kp.readline().split()

        now_kp = [int(k) for k in now_kp]
        ref_kp = [int(k) for k in now_kp]
        tagtag = False
        tt = 0
        with open(outcar_now,'r') as out:
                for line in out :
                        if "length of vectors" in line:
                                tagtag = True
                        if tagtag and tt ==0:
                                tt += 1
                        elif tagtag and tt !=0:
                                ll = line.split()
                                ll = [float(k) for k in ll]
                                tt += 1
                                if tt == 2:
                                        tagtag = False
                                        tt = 0
                                        ref1 = int(math.ceil(max(1.,ll[3]/a_spacing)))
                                        ref2 = int(math.ceil(max(1.,ll[4]/b_spacing)))
                                        ref3 = int(math.ceil(max(1.,ll[5]/c_spacing)))
                                        now_kp[0] = ref1
                                        now_kp[1] = ref2
                                        now_kp[2] = ref3
                                        if ref1 != ref_kp[0] or  ref2 != ref_kp[1] or  ref3 != ref_kp[2] : 
                                                recalc = True
                                        else: 
                                                recalc = False




        return [recalc,now_kp[0],now_kp[1],now_kp[2]]

# Calculate the change of lattice parameters
def calc_lattice_change(pos_ini,pos_fin):
        import math
        axis_ini = poscar_to_axis(pos_ini)
        axis_fin = poscar_to_axis(pos_fin)
        length_ini = [dist_point(axis_ini[x],[0,0,0]) for x in range(3)]
        length_fin = [dist_point(axis_fin[x],[0,0,0]) for x in range(3)]
        angle_ini = [calc_angle(axis_ini[0],axis_ini[1]),calc_angle(axis_ini[1],axis_ini[2]),calc_angle(axis_ini[2],axis_ini[0])]
        angle_fin = [calc_angle(axis_fin[0],axis_fin[1]),calc_angle(axis_fin[1],axis_fin[2]),calc_angle(axis_fin[2],axis_fin[0])]
        length_diff = max([abs(length_fin[x]-length_ini[x])/length_ini[x] for x in range(3)])
        angle_diff = math.degrees(max([abs(angle_fin[x]-angle_ini[x]) for x in range(3)]))
        return [length_diff,angle_diff]

# Set EDIFFG to meet the condition for pressure and force
def set_ediffg(poscar,force,pressure,energy):
        if pressure < 0 and force < 0 and energy < 0:
                return 0
        elif pressure < 0 and force < 0 and enegy > 0:
                return energy
        else:
                [axis,atom_pos] = read_poscar(poscar)
                volume = calc_volume(axis)
                nion = len(atom_pos)
                conv_press = pressure*volume/nion/1602.17733 # it is identical to VASP precision.
                if conv_press < force:
                        return -1.0*conv_press
                else:
                        return -1.0*force
        
