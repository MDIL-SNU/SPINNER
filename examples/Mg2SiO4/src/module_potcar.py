import os
import subprocess
import shutil

def make_potcar(pot_dir):

    with open("POSCAR","r") as f:
        for i in range(5):
            f.readline()
        atoms = f.readline().split()

    fw = open("POTCAR","a")

    for atom in atoms:
        if atom == 'Cs':
            with open(pot_dir+'/'+atom+"_sv_GW/POTCAR","r") as f:
                for fstr in f.readlines():
                    fw.write(fstr)
        elif atom == 'Ca':
            with open(pot_dir+'/'+atom+"_sv/POTCAR","r") as f:
                for fstr in f.readlines():
                    fw.write(fstr)
        elif os.path.isdir(pot_dir+"/"+atom):
            with open(pot_dir+'/'+atom+"/POTCAR","r") as f:
                for fstr in f.readlines():
                    fw.write(fstr)
        elif os.path.isdir(pot_dir+"/"+atom+'_pv'):
            with open(pot_dir+'/'+atom+"_pv/POTCAR","r") as f:
                for fstr in f.readlines():
                    fw.write(fstr)
        elif os.path.isdir(pot_dir+"/"+atom+'_sv'):
            with open(pot_dir+'/'+atom+"_sv/POTCAR","r") as f:
                for fstr in f.readlines():
                    fw.write(fstr)

    fw.close()
