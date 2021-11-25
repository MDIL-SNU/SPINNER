###########################################
### Date: 2019-12-17			###
### yybbyb@snu.ac.kr			###
###########################################
import os, sys, subprocess, yaml, shutil, platform
from input_conf import *
from module_amp2_input import *
from module_log import *
# Check python version
py_ver = platform.sys.version.split()[0].split('.')
if '3' in py_ver[0]:
	print ('Python3 is used.')
elif '7' in py_ver[1]:
	print ('Python2.7 is used.')
else:
	print ('You must use python3 or python2.7. Please check your system.')
	sys.exit()


# input from shell
conf_file = sys.argv[1]
node = os.path.abspath(sys.argv[2])
nproc = sys.argv[3]
encut_input = sys.argv[4]
pypath = sys.executable

home = os.getcwd()

inp_file = input_conf(conf_file)
with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)
cal_dic = inp_yaml['calculation']
src_path = inp_yaml['directory']['src_path']
ERROR_path = inp_yaml['directory']['error']
Done_path = inp_yaml['directory']['done']
large_off = inp_yaml['calculation']['large_off']

calc_list = make_list(inp_file)
calc_out = 0

while len(make_list(inp_file)) > 0:
	target_mat = make_list(inp_file)[0]
	try:
		target = subprocess.check_output([pypath,src_path+'/amp2_input.py',inp_file,node,target_mat[0],target_mat[1],encut_input],universal_newlines=True).splitlines()[0]
	except:
		target = '0'
	if target == '0':
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		continue
	if isinstance(large_off,int) and large_off > 1:
		num_pos = len(read_poscar(target+'/INPUT0/POSCAR')[1])
		if num_pos > large_off:
			make_amp2_log(target,'The number of atoms is larger than cutoff value. We pass this material due to computational cost.')
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue

	if set_on_off(cal_dic['kp_test']) == 1:
		try:
			notice = subprocess.check_output([pypath,src_path+'/kpoint.py',target,inp_file,node,nproc],universal_newlines=True)
		except:
			notice = '0'
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue
	# check existance of follow calculation and K-pts file
	if 1 in list(cal_dic.values()) and not os.path.isfile(target+'/INPUT0/KPOINTS'):
		make_amp2_log(target,'Warning!!! KPOINTS file should be located in INPUT0.')
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		continue
	if set_on_off(cal_dic['encut_test']) == 1:
		try:
			notice = subprocess.check_output([pypath,src_path+'/cutoff.py',target,inp_file,node,nproc],universal_newlines=True)
		except:
			notice = '0'
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue
	if set_on_off(cal_dic['relaxation']) == 1:
		pot_type = inp_yaml['magnetic_ordering']['potential_type']
		try:
				notice = subprocess.check_output([pypath,src_path+'/relax.py',target,inp_file,node,nproc,pot_type],universal_newlines=True)
		except:
				notice = '0'
		if not notice.splitlines()[-1][0] == '1':
				shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
				continue
	# Extract calculation results in Results directory
	subprocess.call([pypath,src_path+'/get_result.py',target])
	shutil.move(target,Done_path+'/'+target.split('/')[-1])
