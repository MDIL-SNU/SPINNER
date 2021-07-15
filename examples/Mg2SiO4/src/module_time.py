import os,sys,subprocess

def get_time(gen, pop_num, src_path, output_dir):

    T = []

    for pop in range(1,pop_num+1):

        os.chdir(src_path+"/../"+output_dir+"/generation"+str(gen)+"/pop"+str(pop))

        f = open("lammps.out","r")
        time_ = f.readlines()[-1].split()[-1].split(':')
        f.close()

        time  = float(time_[0])*3600.0 + float(time_[1])*60.0 + float(time_[2])
        T.append(time)

    return max(T), float(sum(T)/float(len(T)))

#src_path  = os.getcwd()
#print (get_cpu_time(1, 30, src_path, "output"))
