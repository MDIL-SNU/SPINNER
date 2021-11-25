import numpy as np
import math
import random
pi=3.14159265358979
rad2deg=180.0/pi
 
def asin(a):
        return math.atan2(a,math.sqrt(1.0-a*a))
def acos(a):
        return pi/2.0-asin(a)
def norm(x):
        return (math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))
def dotprod (x,y):
        return ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] )
def angle (v1,v2):
        myacos = dotprod(v1,v2)/norm(v1)/norm(v2)
        if myacos>1.0 :
            myacos = 1.0
        if myacos<-1.0:
            myacos = -1.0
        return(acos(myacos)*180.0/3.14159265358979)

def choose_celltype(space_group):
    if space_group < 3:
        return 'triclinic'
    elif space_group < 16:
        return 'monoclinic'
    elif space_group < 75:
        return 'triclinic'
    elif space_group < 143:
        return 'tetragonal'
    elif space_group < 168:
        return 'trigonal'
    elif space_group < 195:
        return 'hexagonal'
    else:
        return 'cubic'

def check_tolerance(latt0, coor0, tolerance_matrix, atomnumlist, change_coor, scale_factor, inp_file, longrange_yes, v, Vmin, number_output):
    
    tot_atom_num = len(coor0)

    if change_coor == 1:
        # restate latt and coor in lammps way
        latt = [[0.0 for i in range(3)] for j in range(3)]
        for i in range(3):
            for j in range(3):
                latt[i][j] = latt0[2-i][2-j]
        coor = [[0.0 for i in range(3)] for j in range(tot_atom_num)]
        for i in range(tot_atom_num):
            coor[i][0] = coor0[i][2]
            coor[i][1] = coor0[i][1]
            coor[i][2] = coor0[i][0]

        # change to cartesian coordinate
        latt = np.array(latt)
        coor = np.array(coor)
        coor = coor.dot(latt)

    else:
        # restate latt and coor in lammps way
        latt = [[0.0 for i in range(3)] for j in range(3)]
        for i in range(3):
            for j in range(3):
                latt[i][j] = latt0[i][j]
        coor = [[0.0 for i in range(3)] for j in range(tot_atom_num)]
        for i in range(tot_atom_num):
            coor[i][0] = coor0[i][0]
            coor[i][1] = coor0[i][1]
            coor[i][2] = coor0[i][2]


    # make type array
    atom_type_num = len(tolerance_matrix)

    atom_type_array = []
    member = -1
    for i in atomnumlist:
        member = member + 1
        for j in range(i):
            atom_type_array.append(int(member))

    #print (1)
    success = 1

    mindis = 10000.0 

   # check the distance tolerance
    if scale_factor != 0.0:
      for atom1 in range(tot_atom_num):
        for atom2 in range(atom1,tot_atom_num):
            R = []
            for i in range(-1,2):
                for j in range(-1,2):
                    for k in range(-1,2):
                        if (atom1 == atom2) and (i==0 and j ==0 and k==0):
                            R.append(10000.0)
                        else:
                            I = float(i)
                            J = float(j)
                            K = float(k)
                            coorx=coor[atom2][0]+I*latt[0][0]+J*latt[1][0]+K*latt[2][0]
                            coory=coor[atom2][1]+I*latt[0][1]+J*latt[1][1]+K*latt[2][1]
                            coorz=coor[atom2][2]+I*latt[0][2]+J*latt[1][2]+K*latt[2][2]
                            coor2 = np.array([coorx,coory,coorz])
                            R.append(np.linalg.norm(coor2-coor[atom1]))
                        
            type1 = atom_type_array[atom1]
            type2 = atom_type_array[atom2]
            if min(R) < tolerance_matrix[type1][type2]*scale_factor:
                success = 0

            if min(R) < mindis:
                mindis = min(R) 

    # check the long range tolerance
    if inp_file['vacuum_constraint']['apply_vacuum_constraint'] == True and success !=0 and longrange_yes == 1:
        # make supercell coordinate
        coor_sup = []

        for atom2 in range(0,tot_atom_num):
            coor_sup.append([coor[atom2][0],coor[atom2][1],coor[atom2][2]])

        for atom2 in range(0,tot_atom_num):
            for i in range(-1,2):
                for j in range(-1,2):
                    for k in range(-1,2):
                        if not (i ==0 and j==0 and k ==0):
                          I = float(i)
                          J = float(j)
                          K = float(k)
                          coorx=coor[atom2][0]+I*latt[0][0]+J*latt[1][0]+K*latt[2][0]
                          coory=coor[atom2][1]+I*latt[0][1]+J*latt[1][1]+K*latt[2][1]
                          coorz=coor[atom2][2]+I*latt[0][2]+J*latt[1][2]+K*latt[2][2]
                          coor_sup.append([coorx,coory,coorz])

        coor_sup = np.array(coor_sup)

        # main loop
        coorx=coor[atom2][0]
        coory=coor[atom2][1]
        coorz=coor[atom2][2]
        lmax = inp_file['vacuum_constraint']['maximum_vacuum_length']
        inspection_grid = inp_file['vacuum_constraint']['grid']
        gridx = int(np.linalg.norm(latt[0])/inspection_grid) 
        gridy = int(np.linalg.norm(latt[1])/inspection_grid)
        gridz = int(np.linalg.norm(latt[2])/inspection_grid)

        shiftx = float(random.random()/gridx)
        shifty = float(random.random()/gridy)
        shiftz = float(random.random()/gridz)

        for T in range(gridx*gridy*gridz):
              i = T%gridx
              T1 = int((T-i)/gridx)
              j = int(T1%gridy)
              k = int((T1-j)/gridy)
              # Cartesian coordinate of the void point
              void_coorx = float(i/gridx+shiftx)*latt[0][0] + float(j/gridy+shifty)*latt[1][0] + float(k/gridz+shiftz)*latt[2][0]
              void_coory = float(i/gridx+shiftx)*latt[0][1] + float(j/gridy+shifty)*latt[1][1] + float(k/gridz+shiftz)*latt[2][1]
              void_coorz = float(i/gridx+shiftx)*latt[0][2] + float(j/gridy+shifty)*latt[1][2] + float(k/gridz+shiftz)*latt[2][2]

              void_coor = np.array([void_coorx,void_coory,void_coorz])

              Rmin = 100.0
              nn1 = -1
              nn2 = -1
              fail = 1
              for l in range(len(coor_sup)):
                Rtemp = np.linalg.norm(void_coor-coor_sup[l])
                if Rtemp < lmax/2.0:
                  fail = 0
                  break
 
              if fail == 1:
                  success = -1
                  break
 
    # convert into lammps format
    p_a= math.sqrt(latt[0][0]**2.0 + latt[0][1]**2.0 + latt[0][2]**2.0);
    p_b= math.sqrt(latt[1][0]**2.0 + latt[1][1]**2.0 + latt[1][2]**2.0);
    p_c= math.sqrt(latt[2][0]**2.0 + latt[2][1]**2.0 + latt[2][2]**2.0);
    alpha=  angle(latt[1],latt[2]);  beta=  angle(latt[0],latt[2]); gamma=  angle(latt[0],latt[1]);  # Angles in degree
    alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians

    lx=   p_a
    p_xy= p_b * math.cos(gammar)
    p_xz= p_c * math.cos(betar)
    ly=   math.sqrt(p_b**2.0 - p_xy**2.0)
    p_yz= (p_b*p_c*math.cos(alphar)-p_xy*p_xz)/(ly)
    lz=   math.sqrt(p_c**2.0 - p_xz**2.0 - p_yz**2.0)

    lattr = [lx,ly,lz,p_xy,p_xz,p_yz]

    # print the mindis coor
    atomtype = []
    for i in range(len(atomnumlist)):
        atomtype += [i+1]*atomnumlist[i]

    with open("coo_result"+str(number_output),"w") as f:
        f.write(" RDF lammps\n")
        f.write("\n")
        f.write(str(sum(atomnumlist))+" atoms\n")
        f.write(str(len(atomnumlist))+" atom types\n")
        f.write("\n")
        f.write("0.000000 %10.6f xlo xhi\n"%(lx/mindis))
        f.write("0.000000 %10.6f ylo yhi\n"%(ly/mindis))
        f.write("0.000000 %10.6f zlo zhi\n"%(lz/mindis))
        f.write("\n")
        f.write("%10.6f %10.6f %10.6f  xy xz yz\n"%(p_xy/mindis, p_xz/mindis, p_yz/mindis))
        f.write("\n")
        f.write("Atoms\n")
        f.write("\n")
        for i in range(sum(atomnumlist)):
            f.write(" %d %d  %10.6f %10.6f %10.6f\n"%(i+1, atomtype[i], coor[i][0]/mindis, coor[i][1]/mindis, coor[i][2]/mindis))

    return success, lattr, coor, mindis
