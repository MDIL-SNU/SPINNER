import sys
import re
import math
import os, sys, random, datetime
import subprocess
# import shutil

ion_radii_table = {"Zn2":0.88,"Pb4":0.915,"Ba2":1.49,"Pb2":1.33,"Ca2":1.14,"S-2":1.7,\
                    "F-1":1.19,"Ta4":0.82,"Ta5":0.78,"Cl7":0.41,"Ta3":0.86,"O-2":1.26,\
                    "Ge2":0.87,"Mo3":0.83,"F7":0.22,"Ga3":0.76,"Mo6":0.73,"Mo5":0.75,\
                    "Mo4":0.79,"Mg2":0.86,"W5":0.76,"Cl-1":1.67,"B3":0.41,"Bi5":0.9,\
                    "Bi3":1.17,"Tl3":1.025,"Tl1":1.64,"Pd4":0.755,"Pd3":0.9,"Pd2":1.0,\
                    "Pd1":0.73,"Sr2":1.32,"Ag3":0.89,"Ag2":1.08,"Ag1":1.29,"Os8":0.53,\
                    "W6":0.74,"Tc7":0.7,"Tc4":0.785,"Tc5":0.74,"Os5":0.715,"Os4":0.77,\
                    "Os7":0.665,"Os6":0.685,"W4":0.8,"Pt5":0.71,"Pt4":0.765,"Pt2":0.94,\
                    "Sc3":0.885,"P3":0.58,"Lu3":1.001,"Te-2":2.07,"P5":0.52,"Hg2":1.16,\
                    "Re6":0.69,"Re7":0.67,"Re4":0.77,"Re5":0.72,"Hg1":1.33,"Cs1":1.81,\
                    "Sb3":0.9,"Sb5":0.74,"H1":-0.04,"Ru8":0.5,"Ru3":0.82,"Ru7":0.52,\
                    "Ru4":0.76,"Ru5":0.705,"S6":0.43,"Rb1":1.66,"S4":0.51,"K1":1.52,\
                    "Be2":0.59,"Nb4":0.82,"Se6":0.56,"Se4":0.64,"Nb3":0.86,"Nb5":0.78,\
                    "In3":0.94,"Te4":1.11,"Te6":0.7,"C4":0.3,"Au1":1.51,"Au3":0.99,\
                    "Au5":0.71,"Cl5":0.26,"N3":0.3,"N5":0.27,"I-1":2.06,"Br-1":1.82,\
                    "Ge4":0.67,"Hf4":0.85,"I5":1.09,"N-3":1.32,"Br7":0.53,"Br5":0.45,\
                    "Br3":0.73,"I7":0.67,"Cd2":1.09,"Xe8":0.62,"Al3":0.675,"Zr4":0.86,\
                    "Si4":0.54,"Ti4":0.745,"Ir3":0.82,"Ir4":0.765,"Ir5":0.71,"Ti2":1.0,\
                    "Ti3":0.81,"Na1":1.16,"Li1":0.9,"Se-2":1.84,"As3":0.72,"As5":0.6,\
                    "Rh3":0.805,"Rh5":0.69,"Rh4":0.74,"Sn4":0.83,"Y3":1.04,'V5':0.68,'Cu2':0.87,\
                    "Ce3":1.01,"Pr3":0.99,"Nd3":0.983,"Pm3":0.97,"Sm3":0.958,"Eu3":0.947,\
                    "Gd3":0.935,"Tb3":0.923,"Dy3":0.912,"Ho3":0.901,"Er3":0.89,"Tm3":0.88,"Tb3":0.868,\
                    "Cr3":0.615,"Cr6":0.44,"Mn2":0.67,"Mn4":0.53,"Fe2":0.61,"Fe3":0.55,"Co2":0.65,"Co3":0.545,\
                    "Ni2":0.69,"Ac3":1.065,"Th4":0.94,"Pa5":0.78,"U4":1.025,"U6":0.76,"Np5":0.87,"Yb3":1.008}

non_metal = ['H','O','F','S','Cl','Se','Br','Te','I','N','C','Si','P','Ge','As','Sn','Sb']

zero = ['He','Ne','Ar','Xe']
one = ['Li','Na','K','Rb','Cs']
two = ['Be','Mg','Ca','Sr','Ba']
three = ['B','Al','Ga','In','Sc','Y','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb']

oxidation = {}

oxidation['H'] = [1,-1]
oxidation['D'] = [1,-1]

oxidation['C'] = [4,3,2,1,0,-1,-2,-3,-4]
oxidation['N'] = [3, -3, 5]
oxidation['O'] = [-2]
oxidation['F'] = [-1]

oxidation['Si'] = [4, -4]
oxidation['P'] = [5, 3, -3]
oxidation['S'] = [6, -2, 2, 4]
oxidation['Cl'] = [-1, 1, 3, 5, 7]

oxidation['Ti'] = [4]
oxidation['Ge'] = [2, 4, -4]
oxidation['As'] = [3, -3, 5]
oxidation['Se'] = [4, 6, 2, -2]
oxidation['Br'] = [-1, 1, 3, 5]
oxidation['Kr'] = [0, 2]

oxidation['Zr'] = [4]
oxidation['Nb'] = [5]
oxidation['Mo'] = [6, 4]
oxidation['Tc'] = [7, 4]
oxidation['Ru'] = [4, 3]
oxidation['Rh'] = [3]
oxidation['Pd'] = [2, 4, 0]
oxidation['Ag'] = [1]
oxidation['Cd'] = [2]
oxidation['Sn'] = [4, 2, -4]
oxidation['Sb'] = [3, 5, -3]
oxidation['Te'] = [4, 6, 2, -2]
oxidation['I'] = [-1, 1,3, 5, 7]
oxidation['Zn'] = [2]

oxidation['Hf'] = [4]
oxidation['Ta'] = [5]
oxidation['W'] = [6, 4]
oxidation['Re'] = [4]
oxidation['Os'] = [4]
oxidation['Ir'] = [3, 4]
oxidation['Pt'] = [4, 2]
oxidation['Au'] = [3, 1]
oxidation['Hg'] = [2, 1]
oxidation['Tl'] = [1, 3]
oxidation['Pb'] = [2, 4]
oxidation['Bi'] = [3]

oxidation['Lu'] = [3]
oxidation['La'] = [3]

oxidation['V'] = [5]
oxidation['Cr'] = [3,6]
oxidation['Mn'] = [2,4]
oxidation['Fe'] = [2,3]
oxidation['Co'] = [2,3]
oxidation['Ni'] = [2]
oxidation['Cu'] = [2]

oxidation['Ac'] = [3]
oxidation['Th'] = [4]
oxidation['Pa'] = [5]
oxidation['U'] = [4,6]
oxidation['Np'] = [5]

for i in zero:
    oxidation[i] = [0]
for i in one:
    oxidation[i] = [1]
for i in two:
    oxidation[i] = [2]

for i in three:
    oxidation[i] = [3]


ion_radii_table['Li0'] = 145.*0.01
ion_radii_table['Be0'] = 105.*0.01
ion_radii_table['B0'] =  85.*0.01
ion_radii_table['Na0'] = 180.*0.01
ion_radii_table['Mg0'] = 150.*0.01
ion_radii_table['Al0'] = 125.*0.01
ion_radii_table['K0'] = 220.*0.01
ion_radii_table['Ca0'] = 180.*0.01
ion_radii_table['Sc0'] = 160.*0.01
ion_radii_table['Ti0'] = 140.*0.01
ion_radii_table['Ga0'] = 130.*0.01
ion_radii_table['Ge0'] = 125.*0.01
ion_radii_table['Rb0'] = 235.*0.01
ion_radii_table['Sr0'] = 200.*0.01
ion_radii_table['Y0'] = 180.*0.01
ion_radii_table['Zr0'] = 155.*0.01
ion_radii_table['Nb0'] = 145.*0.01
ion_radii_table['Mo0'] = 145.*0.01
ion_radii_table['Ru0'] = 130.*0.01
ion_radii_table['Rh0'] = 135.*0.01
ion_radii_table['Pd0'] = 140.*0.01
ion_radii_table['Ag0'] = 160.*0.01
ion_radii_table['Cd0'] = 155.*0.01
ion_radii_table['In0'] = 155.*0.01
ion_radii_table['Sn0'] = 145.*0.01
ion_radii_table['Sb0'] = 145.*0.01
ion_radii_table['Cs0'] = 260.*0.01
ion_radii_table['Ba0'] = 215.*0.01
ion_radii_table['Lu0'] = 227.*0.01
ion_radii_table['Hf0'] = 155.*0.01
ion_radii_table['Ta0'] = 145.*0.01
ion_radii_table['W0'] = 135.*0.01
ion_radii_table['Re0'] = 135.*0.01
ion_radii_table['Os0'] = 130.*0.01
ion_radii_table['Ir0'] = 135.*0.01
ion_radii_table['Pt0'] = 135.*0.01
ion_radii_table['Au0'] = 135.*0.01
ion_radii_table['Hg0'] = 150.*0.01
ion_radii_table['Tl0'] = 190.*0.01
ion_radii_table['Pb0'] = 180.*0.01
ion_radii_table['Bi0'] = 160.*0.01
ion_radii_table['La0'] = 195.*0.01
ion_radii_table['Cu0'] = 135.*0.01
ion_radii_table['Zn0'] = 135.*0.01

ion_radii_table['V0'] = 135.*0.01
ion_radii_table['Cr0'] = 140.*0.01
ion_radii_table['Mn0'] = 140.*0.01
ion_radii_table['Fe0'] = 140.*0.01
ion_radii_table['Co0'] = 135.*0.01
ion_radii_table['Ni0'] = 135.*0.01
ion_radii_table['Ce0'] = 185.*0.01
ion_radii_table['Pr0'] = 185.*0.01
ion_radii_table['Nd0'] = 185.*0.01
ion_radii_table['Pm0'] = 185.*0.01
ion_radii_table['Sm0'] = 185.*0.01
ion_radii_table['Eu0'] = 185.*0.01
ion_radii_table['Gd0'] = 180.*0.01
ion_radii_table['Tb0'] = 175.*0.01
ion_radii_table['Dy0'] = 175.*0.01
ion_radii_table['Ho0'] = 175.*0.01
ion_radii_table['Er0'] = 175.*0.01
ion_radii_table['Tm0'] = 175.*0.01
ion_radii_table['Yb0'] = 175.*0.01
ion_radii_table['Ac0'] = 195.*0.01
ion_radii_table['Th0'] = 180.*0.01
ion_radii_table['Pa0'] = 180.*0.01
ion_radii_table['U0'] = 175.*0.01
ion_radii_table['Np0'] = 175.*0.01


atom_vol = {'H':5.08,'Li':22.6,'Be':36,'B':13.24,'C':13.87,'N':11.8,'O':11.39,'F':11.17,\
            'Na':26,'Mg':36,'Al':39.6,'Si':37.3,'P':29.5,'S':25.2,'Cl':25.8,'K':36,'Ca':45,\
            'Sc':42,'Ti':27.3,'V':24,'Cr':28.1,'Mn':31.9,'Fe':30.4,'Co':29.4,'Ni':26,\
            'Cu':26.9,'Zn':39,'Ga':37.8,'Ge':41.6,'As':36.4,'Se':30.3,'Br':32.7,'Rb':42,\
            'Sr':47,'Y':44,'Zr':27,'Nb':37,'Mo':38,'Tc':38,'Ru':37.3,'Rh':31.2,'Pd':35,\
            'Ag':35,'Cd':51,'In':55,'Sn':52.8,'Sb':48,'Te':46.7,'I':46.2,'Xe':45,'Cs':46,\
            'Ba':66,'La':58,'Ce':54,'Pr':57,'Nd':50,'Sm':50,'Eu':53,'Gd':56,'Tb':45,'Dy':50,\
            'Ho':42,'Er':54,'Tm':49,'Yb':59,'Lu':35,'Hf':40,'Ta':43,'W':38.8,'Re':42.7,\
            'Os':41.9,'Ir':34.3,'Pt':38,'Au':43,'Hg':38,'Tl':54,'Pb':52,'Bi':60,'Ac':74,\
            'Th':56,'Pa':60,'U':58,'Np':45,'Am':17}

def create_rand_structure(composition, pot_dir, target_num):
    at = re.compile('[A-Z][a-z]?')
    elements = at.findall(composition)
    number = list(map(int, re.findall("\d+", composition)))

    # 1.Calculate stoichiometry
    GCD = number[0]
    if len(number) >= 2:
        for i in number[1:]:
            GCD = math.gcd(GCD, i)
    smallest = [i//GCD for i in number]
    multiple = target_num//sum(smallest)
    if multiple*sum(smallest) < target_num:
            multiple += 1
    int_number = [i*multiple for i in smallest]
    float_number = [float(i) for i in int_number]
    print(int_number)
    print(float_number)

    # 2.Find neutral oxidation number in given stoichiometry
    is_contain_non_metal = False
    for j, jtem in enumerate(elements):
        if jtem in non_metal :
            is_contain_non_metal = True

    index = [0]*len(elements)
    reference = [len(oxidation[i]) for i in elements]

    oxi_state = []
    if is_contain_non_metal:
        tag = True
        while tag:
            oxi_tmp = [oxidation[elements[i]][index[i]] for i in range(len(elements))]
            charge = 0
            for i in range(len(elements)):
                charge += oxi_tmp[i]*int_number[i]

            if charge == 0: # satisfy neutral condition
                tag = False
                for i, item in enumerate(elements):
                    oxi_state.append(oxi_tmp[i])
            else:
                index[-1] += 1
                index_valid = [index[i] < reference[i] for i in range(len(index))]
                while False in index_valid:
                    point = index_valid.index(False)
                    index[point-1] += 1
                    for i in range(len(index)-point):
                        index[point+i] = 0
                    index_valid = [index[i] < reference[i] for i in range(len(index))]
                    if False in index_valid:
                        if index_valid.index(False) == 0:
                            tag = False
                            break
    else: # No non metal
        for i, item in enumerate(elements):
            oxi_state.append(0)

    str_oxi_state = [str(i) for i in oxi_state]
    print(f'oxi state : {oxi_state}')
    # 3.Calculate volume of cell using oxidation number info
    if len(oxi_state) == 0 :
        vol = 0
        atom_num = 0
        for i, item in enumerate(elements):
            vol += atom_vol[item]*float_number[i]
            atom_num += float_number[i]

        atom_vol0 = (vol/atom_num) / 1.391
        print("Net charge cannot be 0 with common oxidation states.")
        print("Volume predicted by atomic volume. vol/atom=", round(atom_vol0, 5))
        print("Total vol=", round(atom_vol0*atom_num, 1))
        print()
    else:
        if all(element+state in ion_radii_table for element, state in zip(elements, str_oxi_state)):
            # Code block to execute if all elements and states are present in ion_radii_table
            ion_vol = []
            vol = 0
            atom_num = 0
            for i,item in enumerate(elements):
                vol += (4./3.)*math.pi*float_number[i]*(ion_radii_table[item+str_oxi_state[i]]**3)
                atom_num += float_number[i]
            vol = vol/sum(float_number)
            atom_vol0 = vol/0.6654

            print("Predicted oxidation states:",[elements[i] for i in range(len(elements))],[oxi_state[i] for i in range(len(elements))])
            print("Volume predicted by ionic radius. vol/atom=", round(atom_vol0,5))
            print("Total vol=", round(atom_vol0*atom_num,1))
            print()
        else:
            vol = 0
            atom_num = 0

            for i,item in enumerate(elements):
                vol += atom_vol[item]*float_number[i]
                atom_num += float_number[i]

            atom_vol0 = (vol/atom_num) / 1.391
            print("Predicted oxidation states:",[elements[i] for i in range(len(elements))],[oxi_state[i] for i in range(len(elements))])
            print("No ionic radius available.")
            print("Volume predicted by atomic volume. vol/atom=", round(atom_vol0,5))
            print("Total vol=", round(atom_vol0*atom_num,1))
            print()

    print("--------------------")
    print("structure generation" )
    print("--------------------")

    V = atom_vol0*atom_num
    density = V

    atomname = []
    for atomi in elements:
        atomname.append(atomi)
    natom = int_number

    # 4.Make POTCAR
    print("-------------")
    print("making potcar")
    print("-------------")
    for i in range(len(atomname)):
        atomi = atomname[i]
        if atomi=="Cs":
            atomname[i] +="_sv_GW"
        elif  atomi=="Ca":
            atomname[i] +="_sv"
        elif atomi in ['Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']:
            atomname[i] += '_3'
        elif os.path.isdir(pot_dir+'/'+atomi) == False:
            if os.path.isdir(pot_dir+'/'+atomi+"_pv"):
                atomname[i] = atomname[i]+"_pv"
            elif os.path.isdir(pot_dir+'/'+atomi+"_sv"):
                atomname[i] = atomname[i]+"_sv"
            else:
                print("ERROR! no atom name in potcar")

    if len(atomname)==1:
        subprocess.call('cat '+pot_dir+"/"+atomname[0]+'/POTCAR > POTCAR',shell=True)
    else:
        subprocess.call('cat '+pot_dir+"/"+atomname[0]+'/POTCAR > POTCAR',shell=True)
        for atomi in atomname[1:]:
            subprocess.call('cat '+pot_dir+"/"+atomi+'/POTCAR >> POTCAR',shell=True)

    # 5.Create random positions
    def distance(a,b) :
        return sum([(x-y)**2.0 for x,y in zip(a,b)])**0.5 ;

    def dist_pbc(a,b,cell) :
        imagecell = [-1,0,1]
        pbc = [[i,j,k] for i in imagecell for j in imagecell for k in imagecell]
        b_pbc = [[b[i] + cell*pbc[j][i] for i in range(3)] for j in range(len(pbc))]
        return min([distance(a,b_pbc[i]) for i in range(len(pbc))])

    nelement = len(atomname)
    totatom = sum(natom)
    cutoff = [[1.5 for i in range(nelement)] for j in range(nelement)]

    exatom_position = [[0 for i in range(3)] for j in range(totatom)]
    cellpar = V**(1.0/3.0)
    tot_attempt = 0
    for newatom in range(totatom):
        newatom_type = sum([newatom >= sum(natom[0:j+1]) for j in range(nelement)])
        newatom_position = [cellpar * random.random() for i in range(3)]
        exatom = -1
        while exatom < newatom:
            tot_attempt = tot_attempt + 1
            if tot_attempt > 100000:     # Exit Loop if it takes too long
                sys.exit()
            newatom_position = [cellpar * random.random() for i in range(3)]
            for exatom in range(newatom+1):
               exatom_type = sum([exatom >= sum(natom[0:j+1]) for j in range(nelement)])
               dist = dist_pbc(newatom_position,exatom_position[exatom],cellpar)
               if dist < cutoff[newatom_type][exatom_type]:
                   break
        exatom_position[newatom] = newatom_position

    print(tot_attempt)
    print(cellpar)

    # 6.Writing POSCAR (Converting Fractional to Direct coordinate)
    with open('POSCAR','w') as POSCAR:
        POSCAR.write(''.join(atomname) + '   density: '+str(density)+'   ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        POSCAR.write("\n")
        POSCAR.write("   1.0\n")
        POSCAR.write("{:21.13f}{:19.13f}{:19.13f}\n".format(cellpar,0,0))
        POSCAR.write("{:21.13f}{:19.13f}{:19.13f}\n".format(0,cellpar,0))
        POSCAR.write("{:21.13f}{:19.13f}{:19.13f}\n".format(0,0,cellpar))
        POSCAR.write('   ' + '   '.join(elements)+"\n")
        POSCAR.write('    ' + '    '.join([str(x) for x in natom]) + '  \n')
        POSCAR.write('Selective dynamics\n')
        POSCAR.write('Direct\n')

        exatom_position = [[exatom_position[i][j]/cellpar for j in range(3)] for i in range(totatom)]
        for i in range(totatom-1):
            POSCAR.write("{:19.15f}{:19.15f}{:19.15f}   T   T   T  \n".format(exatom_position[i][0],exatom_position[i][1],exatom_position[i][2]))
        POSCAR.write("{:19.15f}{:19.15f}{:19.15f}   F   F   F  \n".format(exatom_position[totatom-1][0],exatom_position[totatom-1][1],exatom_position[totatom-1][2]))

if __name__ == "__main__":
    #pot_dir =  '/data/vasp4us/pot/PBE54'
    #target_num = 64
    #material = 'Al6Fe12Mo1Ti5'
    #create_rand_structure(material, pot_dir, target_num)
    pass
