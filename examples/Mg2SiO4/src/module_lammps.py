import lammps
from lammps import lammps
import os
import math
import module_spacegroup
from module_spacegroup import *
import gc

def start_lammps(atomnamelist,lmp,number): 
    lmp.command("units           metal")
    lmp.command("newton          on")
    lmp.command("dimension       3")
    lmp.command("boundary        p p p")
    lmp.command("atom_style 	atomic ")
    lmp.command("atom_modify     map yes")
    lmp.command("box tilt large")
    if number == 0:
        lmp.command("read_data coo  ")
    else:
        lmp.command("read_data coo_result{}  ".format(number))
    lmp.command("pair_style nn   ")
    sentence = "pair_coeff * *  potential "
    for atom in atomnamelist:
        sentence +=  atom + " "
    lmp.command(sentence)
    #lmp.command("pair_coeff * *  potential Li Nb S")
    for i in range(len(atomnamelist)):
        sentence = "mass {} 15.0".format(i+1)
        lmp.command(sentence)
    lmp.command("neighbor 2 bin")
    lmp.command("neigh_modify 	every 1 delay 0 check yes")
    lmp.command("compute 	_rg all gyration ")               # All atoms will be involved in the simulation http://lammps.sandia.gov/doc/compute.html 
    gc.collect()


def fix_holdem(atomnumlist, tolerance_matrix,lmp):
    atomnum_array = []
    for i in range(len(atomnumlist)):
        atomnum_array += [i]*atomnumlist[i]

    number = 0
    for i in range(len(atomnum_array)):
        for j in range(i,len(atomnum_array)):
            number += 1
            sentence = "fix holdem"+ str(number) + " all restrain lbound " + str(i+1) + "  " + str(j+1) + " 100.0 100.0 " + str(tolerance_matrix[int(atomnum_array[i])][int(atomnum_array[j])])
            lmp.command(sentence) 
            sentence = "fix_modify holdem"+str(number)+" energy yes"
            lmp.command(sentence) 
    gc.collect()

def relax_with_fixed_lattice(tot_atom_num,lmp,relax_method, relax_iter):
    if relax_method == 'fire':
        lmp.command("min_style fire")
        lmp.command("min_modify dmax 0.1")
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
    elif relax_method == 'quickmin':
        lmp.command("min_style quickmin")
        lmp.command("min_modify dmax 0.1")
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
    elif relax_method == 'cg':
        lmp.command('min_style cg')
        lmp.command("min_modify line quadratic dmax 0.1")
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num)))
    elif relax_method == 'mix':
        lmp.command('min_style cg')
        lmp.command("min_modify line quadratic dmax 0.1")
        lmp.command("minimize 0 0.1  %d 10000"%(int(relax_iter*tot_atom_num/2)))
        lmp.command("min_style fire")
        lmp.command("min_modify dmax 0.05")
        lmp.command("minimize 0 0.02 %d 10000"%(int(relax_iter*tot_atom_num/2)))
    gc.collect()


def relax_with_unfixed_lattice(tot_atom_num,lmp, relax_iter):
    lmp.command("min_style cg")
    lmp.command("min_modify      line quadratic dmax 0.1 ")
    lmp.command("fix 2 all box/relax tri 0.0 vmax 0.001")
    lmp.command("minimize 1.0e-6 0 %d 10000"%(relax_iter*tot_atom_num))
    lmp.command("unfix 2")
    
    lmp.command("fix 3 all box/relax tri 0.0 vmax 0.0001")
    lmp.command("min_modify      line quadratic dmax 0.1 ")
    lmp.command("minimize 1.0e-6 0 %d 10000"%(relax_iter*tot_atom_num))
    lmp.command("unfix 3")
    lmp.command('write_data coo_result2')
    gc.collect()

def relax_with_accurate_potential(tot_atom_num,lmp,atomnamelist):
    lmp.command("units           metal")
    lmp.command("newton          on")
    lmp.command("dimension       3")
    lmp.command("boundary        p p p")
    lmp.command("atom_style 	atomic ")
    lmp.command("atom_modify     map yes")
    lmp.command("box tilt large")
    lmp.command("read_data coo_result2  ")
    lmp.command("pair_style nn   ")
    sentence = "pair_coeff * *  potential_accurate "
    for atom in atomnamelist:
        sentence +=  atom + " "
    lmp.command(sentence)
    lmp.command("neighbor 2 bin")
    lmp.command("neigh_modify 	every 1 delay 0 check yes")
    lmp.command("compute 	_rg all gyration ")               # All atoms will be involved in the simulation http://lammps.sandia.gov/doc/compute.html 

    lmp.command("min_style cg")
    lmp.command("fix 4 all box/relax tri 0.0 vmax 0.0001")
    lmp.command("min_modify      line quadratic dmax 0.1 ")
    lmp.command("minimize 1.0e-6 0 1 10000")
    lmp.command("unfix 4")
    lmp.command('write_data coo_result3')
    gc.collect()
    
def unfixholdem(atomnumlist,lmp):
    atomnum_array = []
    for i in range(len(atomnumlist)):
        atomnum_array += [i]*atomnumlist[i]
    number = 0
    for i in range(len(atomnum_array)):
        for j in range(i,len(atomnum_array)):
            number += 1
            lmp.command("unfix holdem"+str(number))
    del atomnum_array
    gc.collect()

def evaluation_E(lmp):
    lmp.command("run 0")
    lmp.command("variable e equal etotal")
    lmp.command("variable v equal vol")
    e1 = lmp.extract_variable("e",0,0)
    v1 = lmp.extract_variable("v",0,0)
    return e1, v1
    gc.collect()

def get_latt_coor(tot_atom_num,lmp):
    lmp.command("variable v equal vol")
    lmp.command("variable e equal etotal")
    lmp.command("variable a equal cellalpha")
    lmp.command("variable b equal cellbeta")
    lmp.command("variable o equal cellgamma")
    lmp.command("variable g equal cella")
    lmp.command("variable h equal cellb")
    lmp.command("variable y equal cellc")
    lmp.command("variable a1 equal $g")
    lmp.command("variable b1 equal $h*cos($o*PI/180)")
    lmp.command("variable b2 equal $h*sin($o*PI/180)")
    lmp.command("variable c1 equal $y*cos($b*PI/180)")
    lmp.command("variable c2 equal $y*(cos($a*PI/180)-cos($b*PI/180)*cos($o*PI/180))/sin($o*PI/180)")
    lmp.command("variable c3 equal $y*sqrt(1+2*cos($a*PI/180)*cos($b*PI/180)*cos($o*PI/180)-cos($a*PI/180)^2-cos($b*PI/180)^2-cos($o*PI/180)^2)/sin($o*PI/180)")
    a1 = lmp.extract_variable("a1",0,0)
    b1 = lmp.extract_variable("b1",0,0)
    b2 = lmp.extract_variable("b2",0,0)
    c1 = lmp.extract_variable("c1",0,0)
    c2 = lmp.extract_variable("c2",0,0)
    c3 = lmp.extract_variable("c3",0,0)

    latt = [[a1,0.0,0.0],[b1,b2,0.0],[c1,c2,c3]]
    x = lmp.extract_atom('x',3)

    coor = []
    for i in range(tot_atom_num):
        coor.append([x[i][0],x[i][1],x[i][2]])

    return latt,coor
    gc.collect()

def calculate_rdf(atomnamelist, rdf_grid,lmp,number):
    lmp.command("units           metal")
    lmp.command("newton          on")
    lmp.command("dimension       3")
    lmp.command("boundary        p p p")
    lmp.command("atom_style 	atomic ")
    lmp.command("atom_modify     map yes")
    lmp.command("box tilt large")
    lmp.command("read_data coo_result{}  ".format(number))
    lmp.command("pair_style nn   ")
    sentence = "pair_coeff * *  potential10 "
    for atom in atomnamelist:
        sentence +=  atom + " "
    lmp.command(sentence)
    lmp.command("neighbor 2 bin")
    lmp.command("neigh_modify 	every 1 delay 0 check yes one 50000 page 500000")
    lmp.command("compute 	_rg all gyration ")               # All atoms will be involved in the simulation http://lammps.sandia.gov/doc/compute.html 
    lmp.command("variable g equal cella")
    lmp.command("variable h equal cellb")
    lmp.command("variable y equal cellc")
    g = lmp.extract_variable("g",0,0)
    h = lmp.extract_variable("h",0,0)
    y = lmp.extract_variable("y",0,0)
    
    if g<15.0:
        r1 = math.ceil(15.0/g)
    else:
        r1 = 1

    if h < 15.0:
        r2 = math.ceil(15.0/h)
    else:
        r2 = 1
    
    if y < 15.0:
        r3 = math.ceil(15.0/y)
    else:
        r3 = 1
    
    lmp.command("replicate "+str(r1)+" "+str(r2)+" "+str(r3))

    for i in range(1,len(atomnamelist)+1):
        for j in range(i,len(atomnamelist)+1):
            sentence = 'compute myRDF{}{} all rdf {} {} {} cutoff 15.0'.format(i,j,int(3.0*rdf_grid),i,j)
            lmp.command(sentence)
            sentence = 'fix myRDF{}{}_ all ave/time 1 1 1 c_myRDF{}{}[*] mode vector'.format(i,j,i,j)
            lmp.command(sentence)
    lmp.command('run 1')

    rdf = {}
    
    for i in range(1,len(atomnamelist)+1):
        for j in range(i,len(atomnamelist)+1):
            rdf[str(i)+str(j)] = []
            for k in range(int(rdf_grid*3.0)):
                rdf[str(i)+str(j)].append([lmp.extract_fix("myRDF"+str(i)+str(j)+"_",0,2,k,1),lmp.extract_fix("myRDF"+str(i)+str(j)+"_",0,2,k,1)])

    return rdf
    gc.collect()

def atomic_energy(tot_atom_num,lmp):
    lmp.command("compute atompe all pe/atom")
    #lmp.command("dump mydump     all custom 1 atomic_energy id c_atompe")
    #lmp.command("dump_modify     mydump sort id")

    lmp.command("minimize 0 0 0 0")
    atom_e = lmp.extract_compute("atompe",1,1)    
    return atom_e[0:tot_atom_num]
    gc.collect()

def run_lammps(Emin, tolerance_matrix, atomnamelist, atomnumlist, inp_file, accurate_potential, Vmin):

    tot_atom_num = inp_file['tot_atom_num']
    Ecut = inp_file['energy_criteria']['energy_cut_for_further_relax']
    rdf_grid = inp_file['similarity_metric']['rdf_grid']
    relax_method = inp_file['relax_condition']['method_of_first_relax']
    relax_iter = inp_file['relax_condition']['relax_iteration']

    lmp = lammps()
    start_lammps(atomnamelist,lmp,0)
    fix_holdem(atomnumlist, tolerance_matrix,lmp)
    relax_with_fixed_lattice(tot_atom_num,lmp,relax_method, relax_iter)
    unfixholdem(atomnumlist,lmp)
    e,v = evaluation_E(lmp) # include unfix holdem
    latt, coor = get_latt_coor(tot_atom_num,lmp)    
    atom_e = atomic_energy(tot_atom_num, lmp)

    lmp.command('write_data coo_result1')
    lmp.close()
    del lmp
    success, latt, coor = check_tolerance(latt, coor, tolerance_matrix, atomnumlist, 0, 1.0, inp_file,0,v,Vmin) 

    if success == 0:
        status = "distance cut in the first relax"
        e = 10000.0
        lmp = lammps()
        rdf0 = calculate_rdf(atomnamelist, rdf_grid,lmp,1)
        lmp.close()
        del lmp
    elif e-Emin < Ecut*tot_atom_num:
        lmp = lammps()
        start_lammps(atomnamelist,lmp,1)
        fix_holdem(atomnumlist, tolerance_matrix,lmp)
        relax_with_unfixed_lattice(tot_atom_num,lmp, relax_iter)
        atom_e = atomic_energy(tot_atom_num, lmp)
        unfixholdem(atomnumlist,lmp)
        e,v = evaluation_E(lmp) # include unfix holdem
        latt, coor = get_latt_coor(tot_atom_num,lmp)
        lmp.close()

        del lmp
        if accurate_potential  == True:
            lmp = lammps()
            relax_with_accurate_potential(tot_atom_num,lmp,atomnamelist)
            atom_e = atomic_energy(tot_atom_num, lmp)
            e,v = evaluation_E(lmp) # include unfix holdem
            latt, coor = get_latt_coor(tot_atom_num,lmp)
            lmp.close()
            del lmp
            lmp = lammps()
            rdf0 = calculate_rdf(atomnamelist, rdf_grid,lmp,3)
            lmp.close()
            del lmp
        else:
            lmp = lammps()
            rdf0 = calculate_rdf(atomnamelist, rdf_grid,lmp,2)
            lmp.close() 
            del lmp
        success, latt, coor = check_tolerance(latt, coor, tolerance_matrix, atomnumlist, 0, 1.0, inp_file,1,v,Vmin)  
        status = ""
        if success == 0:
            status = "distance cut in the second relax"
            e = 10000.0
        if success == -1:
            status = "too large vacuum region"
            e = 10000.0
    else:
        status = "not in the energy range in the first relax"
        lmp = lammps()
        rdf0 = calculate_rdf(atomnamelist, rdf_grid,lmp,1)
        lmp.close()
        del lmp

    return e,v,latt,coor,rdf0,atom_e,status
    gc.collect()

