import os
from os.path import isfile, join

def listfile(mypath):
        onlyfiles = [f for f in os.listdir(mypath) if isfile(join(mypath, f))]
        return onlyfiles

def listDir(mypath):
        onlyfiles = [f for f in os.listdir(mypath) if not isfile(join(mypath, f))]
        return onlyfiles

def clean_files(gen, pop_num, src_path, output_dir):

    for i in range(1,pop_num+1):
#        print (i)

        path=src_path+"/../"+output_dir+"/generation"+str(gen)+"/pop"+str(i)

        os.remove(path+"/log.lammps")
        os.remove(path+"/coo")
        os.remove(path+"/lammps.in")
        os.remove(path+"/potential")
        os.remove(path+"/do.sh")


#        print (t2-t1)

