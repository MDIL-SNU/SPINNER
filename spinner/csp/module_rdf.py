import numpy as np
import sys
import math

# reading POSCAR and saving all information
# CAUTION: ion numbering starting from zero. (python notation)
class poscar:

    class util_box:   # useful functions
        def dot(self,a,b):
            total = 0
            if len(a) != len(b):
                print("Warning on function 'dot': Vector dimension is not the same.")
                raise NotImplementedError
            for i in range(len(a)):
                total += a[i]*b[i]
            return total
        
        def norm(self,a):
            total = 0
            for i in range(len(a)):
                total += a[i]**2
            return math.sqrt(total)
        
        def cross(self,a,b):
            if len(a) != len(b):
                print("Warning on function 'dot': Vector dimension is not the same.")
                raise NotImplementedError
            elif len(a)+len(b) != 6:
                print("Warning on function 'cross': Only 3-d vector is allowed.")
                raise NotImplementedError
            c = [a[1]*b[2] - a[2]*b[1],
                 a[2]*b[0] - a[0]*b[2],
                 a[0]*b[1] - a[1]*b[0]]
            return c
        
        def triple(self,a):
            return self.dot(a[0],self.cross(a[1],a[2]))
        
        def frac2cart(self,LatticeVectorSet,Fractional):
            Decomposed = []
            for i in [0,1,2]:
                Decomposed.append(np.array(LatticeVectorSet[i])*Fractional[i])
            return list(sum(Decomposed))

        def angle(self,vec1,vec2):
            return math.acos(self.dot(vec1,vec2)/(self.norm(vec1)*self.norm(vec2)))

    def __init__(self,path=False,Lines=False):
        if path:
            self.path = path
            self.RawLines = open(self.path).readlines()
        else:
            self.RawLines = Lines
        self.util = self.util_box()

    def poscar(self):
        self.PoscarLines = self.RawLines
        self.decompose()

   # converting XDATCAR to POSCAR form with configuration selected
   # configNum: N^th configuration which is wanted to be seen
    def xdatcar(self,configNum):
        self.XdatLines = self.RawLines
        self.configNum = configNum
        self.PoscarLines = []
        self.name = self.XdatLines[0]
        self.NumIonList = list(map(lambda x: int(x),self.RawLines[6].split()))
        self.numIon = sum(self.NumIonList)
        self.stopSwitch = False
        for cnt in range(len(self.XdatLines)):
            if self.stopSwitch == True:
                break
            elif self.name == self.XdatLines[cnt]:
                self.PoscarLines = self.XdatLines[cnt:cnt+7]
            elif self.XdatLines[cnt][0].upper() == 'D':
                if int(self.XdatLines[cnt].split('=')[-1]) == self.configNum:
                    self.PoscarLines += self.XdatLines[cnt:cnt+self.numIon+1]
                    self.stopSwitch = True
        self.decompose()

   # extracting all useful information from the object of POSCAR form
   # not necessary to execute this function because it is already included in other functions
    def decompose(self):
        self.name = self.PoscarLines[0].split('\n')[0]
        self.mag = float(self.PoscarLines[1])
        self.LatticeVector = []   # [[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]]
        for i in [2,3,4]:
            self.LatticeVector.append([float(self.PoscarLines[i].split()[x])*self.mag for x in range(3)])
        self.LatticeLength = []   # absolute magnitude of lattice vector
        for i in [0,1,2]:
            self.LatticeLength.append(self.util.norm(self.LatticeVector[i]))
        self.LatticeAngle = []   # alpha,beta,gamma in radian
        for [i,j] in [[1,2],[2,0],[0,1]]:
            self.LatticeAngle.append(self.util.angle(self.LatticeVector[i],self.LatticeVector[j]))
        self.volume = self.util.triple(self.LatticeVector)
        self.LatticePlaneAngle = []   # e.g. angle between a lattice and bc plane
        for [i,j,k] in [[0,1,2],[1,2,0],[2,0,1]]:
            self.LatticePlaneAngle.append(math.pi/2-self.util.angle(self.LatticeVector[i],self.util.cross(self.LatticeVector[j],self.LatticeVector[k])))
        self.PlanePlaneAngle = []   # e.g. angle between ab plane and ac plane
        for [i,j,k] in [[0,1,2],[1,2,0],[2,0,1]]:
            self.PlanePlaneAngle.append(self.util.angle(self.util.cross(self.LatticeVector[i],self.LatticeVector[j]),self.util.cross(self.LatticeVector[i],self.LatticeVector[k])))
        self.ElementList = self.PoscarLines[5].split()
        self.NumIonList = list(map(lambda x: int(x),self.PoscarLines[6].split()))
        self.numIon = sum(self.NumIonList)
        if self.PoscarLines[7][0].upper() == 'S':   # if 'selective dynamics' term is included
            self.posStartNum = 9
            self.selective = True
        else:
            self.posStartNum = 8
            self.selective = False
        self.Relative = []
        self.Cartesian = []
        if self.PoscarLines[self.posStartNum-1][0].upper() == 'D':   # Direct coordinate
            self.dorc = 'D'
            for cnt in range(self.posStartNum,self.posStartNum+self.numIon):
                self.Relative.append([float(self.PoscarLines[cnt].split()[x]) for x in range(3)])
#            for j in range(len(self.Relative)):
#                self.Cartesian.append(self.util.frac2cart(self.LatticeVector,self.Relative[j]))
            self.Cartesian = self.relative2cart(self.Relative)
        elif self.PoscarLines[self.posStartNum-1][0].upper() == 'C':   # Cartesian coordinate
            self.dorc = 'C'
            for cnt in range(self.posStartNum,self.posStartNum+self.numIon):
                self.Cartesian.append([float(self.PoscarLines[cnt].split()[x]) for x in range(3)])
#            Aoriginal = np.array(self.LatticeVector)
#            Atrans = np.transpose(Aoriginal)
#            Ainverse = np.linalg.inv(Atrans)
#            for cnt in range(len(self.Cartesian)):
#                self.Relative.append(list(np.dot(Ainverse,np.transpose(self.Cartesian[cnt]))))
            self.Relative = self.cart2relative(self.Cartesian)
       # labeling atoms with element in sequence
        self.ElementLabel = []
        for cnt in range(len(self.NumIonList)):
            self.ElementLabel[sum(self.NumIonList[0:cnt]):sum(self.NumIonList[0:cnt+1])] = [self.ElementList[cnt]]*self.NumIonList[cnt]
        self.Label = []
        for cnt in range(len(self.NumIonList)):
            self.Label[sum(self.NumIonList[0:cnt]):sum(self.NumIonList[0:cnt+1])] = [cnt]*self.NumIonList[cnt]
       # numbering atoms in sequence, starting from zero (python notation).
        self.Numbering = []
        for cnt in range(len(self.Relative)):
            self.Numbering.append(cnt)
       # calculating interplanar distance
        self.InterSpacing = []
        for [i,j] in [[1,2],[2,0],[0,1]]:
            self.InterSpacing.append(self.volume/np.linalg.norm(np.cross(self.LatticeVector[i],self.LatticeVector[j])))

    def cart2relative(self,CartesianList):
        Aoriginal = np.array(self.LatticeVector)
        Atrans = np.transpose(Aoriginal)
        Ainverse = np.linalg.inv(Atrans)
        Relative = []
        for cnt in range(len(CartesianList)):
            Relative.append(list(np.dot(Ainverse,np.transpose(CartesianList[cnt]))))
        return Relative
        
    def relative2cart(self,RelativeList):
        Cartesian = []
        for j in range(len(RelativeList)):
            Cartesian.append(self.util.frac2cart(self.LatticeVector,RelativeList[j]))
        return Cartesian

   # The copy-and-paste of the result on screen can be used as POSCAR itself.
    def display(self,coord='Direct',showElement=True,addSelective=True):
        printer = ''
        printer += '%s\n'%self.name
        printer += '           1\n'
        for row in range(3):
            printer += ' '
            for col in range(3):
                printer += '%12.6f'%self.LatticeVector[row][col]
            printer += '\n'
        if showElement:
            for elm in self.ElementList:
                printer += '%5s'%elm
            printer += '\n'
        for num in self.NumIonList:
            printer += '%4i'%num
        printer += '\n'
        if addSelective:
            printer += 'Selective Dynamics\n'
        if coord[0].upper() == 'D':
            printer += 'Direct\n'
            for cnt in range(len(self.Relative)):
                printer += ' '
                for subCnt in range(len(self.Relative[cnt])):
                    printer += '%12.8f'%self.Relative[cnt][subCnt]
                if addSelective:
                    printer += '     T T T\n'
                else:
                    printer += '\n'
        if coord[0].upper() == 'C':
            printer += 'Cartesian\n'
            for cnt in range(len(self.Cartesian)):
                for subCnt in range(len(self.Cartesian[cnt])):
                    printer += '%12.8f'%self.Cartesian[cnt][subCnt]
                if addSelective:
                    printer += '     T T T\n'
                else:
                    printer += '\n'
        return printer

   # making a copy of atomic coordinate data, in addition to spatially expand the original ones
   # replicaCutoff: the extent to which you want to expand the original cell from its boundaries
   # outputs:
        # self.ReplicaRelative(Cartesian): coordinate of ions in replicated cell, including original ones
        # self.ReplicaElementLabel: elements listed in the same order to self.ReplicaRelative
        # self.ReplicaNumbering: original ion number in the same order to self.ReplicaRelative
    def replica(self,replicaCutoff):
        self.replicaCutoff = float(replicaCutoff)
        self.ReplicaRelative = []
        self.ReplicaCartesian = []
        self.ReplicaElementLabel = []
        self.ReplicaNumbering = []
        def distance(p1,p2):
            return np.linalg.norm(list(map(lambda x,y:x-y,p1,p2)))
        self.Midpoint = list(map(lambda x,y,z:(x+y+z)/2,self.LatticeVector[0],self.LatticeVector[1],self.LatticeVector[2]))
        CenterVertex = []
        CenterVertex.append(np.linalg.norm(self.Midpoint))
        for i in [0,1,2]:
            #CenterVertex.append(np.linalg.norm(map(lambda x,y:x-y,self.Midpoint,self.LatticeVector[i])))
            CenterVertex.append(distance(self.Midpoint,self.LatticeVector[i]))
        self.circumcircleRadius = max(CenterVertex)
        for j in range(-int(replicaCutoff/self.InterSpacing[0])-1,int(replicaCutoff/self.InterSpacing[0])+2):   # expanding cell
            for k in range(-int(replicaCutoff/self.InterSpacing[1])-1,int(replicaCutoff/self.InterSpacing[1])+2):
                for l in range(-int(replicaCutoff/self.InterSpacing[2])-1,int(replicaCutoff/self.InterSpacing[2])+2):
                    factor = [j,k,l]  # =[0,0,0] if the location is on itself
                    for u in range(len(self.Relative)):  # span all atoms in the unitcell
                        Relative = [self.Relative[u][x]+factor[x] for x in range(3)]
                        temp = []
                        for i in [0,1,2]:
                            temp.append(np.array(self.LatticeVector[i])*Relative[i])
                        Cartesian = list(sum(temp))
                        if factor == [0,0,0]:
                            pass
                        #elif self.circumcircleRadius+self.replicaCutoff < np.linalg.norm(map(lambda x,y:x-y,self.Midpoint,Cartesian)):
                        elif self.circumcircleRadius+self.replicaCutoff < distance(self.Midpoint,Cartesian):
                            continue
                        self.ReplicaRelative.append(Relative)
                        self.ReplicaCartesian.append(Cartesian)
                        self.ReplicaElementLabel.append(self.ElementLabel[u])
                        self.ReplicaNumbering.append(self.Numbering[u])

   # finding neighbors and the distance to them
   # result form: [[1st neighbor, distance, [pos1]],[2nd neighbor, distance, [pos2]],...]
    def environ(self,targetAtom,cutoff):
        cutoff = float(cutoff)
        try:
            if cutoff > self.replicaCutoff:
                self.replica(cutoff)
        except:
            self.replica(cutoff)
        Origin = self.Cartesian[targetAtom]
        output = []
        for cnt in range(len(self.ReplicaCartesian)):
            distance = np.linalg.norm(np.array(self.ReplicaCartesian[cnt])-np.array(Origin))
            if distance < cutoff:
                output.append([self.ReplicaNumbering[cnt],distance,self.ReplicaCartesian[cnt]])
        output.sort(key=lambda x:x[1])
        return output[1:]   # excluding itself

   # listing nearest neighbors for all ions using function self.environ
    def neighbor_analysis(self,cutoff):
        self.replica(cutoff)
        NeighborList = []
        for cnt in range(len(self.Cartesian)):
            NeighborList.append(self.environ(cnt,cutoff))
        return NeighborList

   # finding an atom closest to the point, and returning the distance and its ion number
    def closest(self,Point,cutoff=1.0):
        try:
            if self.replicaCutoff != cutoff:
                self.replica(cutoff)
        except:
            self.replica(cutoff)
        closestNum = -1   # initializing
        closestDist = cutoff+1   # init.
        closestVec = [99,99,99]   # init.
        Point = self.util.frac2cart(self.LatticeVector,Point)
        for i in range(len(self.ReplicaCartesian)):
            tempVec = list(self.ReplicaCartesian[i]-np.array(Point))
            tempDist = self.util.norm(tempVec)
            if closestDist > tempDist:
                closestDist = tempDist
                closestNum = self.ReplicaNumbering[i]
                closestVec = tempVec
            if closestDist < cutoff:
                break
        return [closestNum,closestDist,closestVec]

   # finding a (shortest) distance between the two points in the cell, considering periodic boundary condition
    def distance(self,ionNum1,ionNum2):
        cutoff = max(self.InterSpacing)
        try:
            if cutoff > self.replicaCutoff:
                self.replica(cutoff)
        except:
            self.replica(cutoff)
        Origin = self.Cartesian[ionNum1]
        output = []
        for cnt in range(len(self.ReplicaCartesian)):
            if self.ReplicaNumbering[cnt] == ionNum2:
                distance = np.linalg.norm(np.array(self.ReplicaCartesian[cnt])-np.array(Origin))
                output.append([distance,self.Relative[ionNum1],self.ReplicaRelative[cnt]])
        output.sort(key=lambda x:x[0])
        return output[0]   # returning shortest one

   # finding distance and corresponding coordinates of the two arbitrarily given points, considering PBC.
    def distance_arb(self,X1,X2,coord='D'):
        if coord.upper() == 'C':
            [X1,X2] = self.cart2relative([X1,X2])
        Tmp = self.reduce_pbc([X1,X2])
        X1 = Tmp[0]
        X2 = Tmp[1]
        for i in [0,1,2]:
            if X1[i] - X2[i] > 0.5:
                X2[i] += 1
            if X1[i] - X2[i] < -0.5:
                X2[i] -= 1
        (X1,X2) = self.relative2cart([X1,X2])
        Vec = [X1[i]-X2[i] for i in [0,1,2]]
        return self.util.norm(Vec),X1,X2

   # folding back the vertices to the unitcell
    def reduce_pbc(self,RelativeList):
        ReducedList = []
        for cnt in range(len(RelativeList)):
            ReducedList.append(list(map(lambda x: x-math.floor(x), RelativeList[cnt])))
        return ReducedList


class constants:
    atomic_mass = dict(H=1.01, He=4.00, Li=6.94, Be=9.01, B=10.81, C=12.01,
                   N=14.01, O=16.00, F=19.00, Ne=20.18, Na=22.99, Mg=24.31,
                   Al=26.98, Si=28.09, P=30.97, S=32.07, Cl=35.45, Ar=39.95,
                   K=39.10, Ca=40.08, Sc=44.96, Ti=47.87, V=50.94, Cr=52.00,
                   Mn=54.94, Fe=55.85, Co=58.93, Ni=58.69, Cu=63.55, Zn=65.39,
                   Ga=69.72, Ge=72.61, As=74.92, Se=78.96, Br=79.90, Kr=83.80,
                   Rb=85.47, Sr=87.62, Y=88.91, Zr=91.22, Nb=92.91, Mo=95.94,
                   Tc=98.00, Ru=101.07, Rh=102.91, Pd=106.42, Ag=107.87,
                   Cd=112.41, In=114.82, Sn=118.71, Sb=121.76, Te=127.60,
                   I=126.90, Xe=131.29, Cs=132.91, Ba=137.33, La=138.91,
                   Ce=140.12, Pr=140.91, Nd=144.24, Pm=145.00, Sm=150.36,
                   Eu=151.96, Gd=157.25, Tb=158.93, Dy=162.50, Ho=164.93,
                   Er=167.26, Tm=168.93, Yb=173.04, Lu=174.97, Hf=178.49,
                   Ta=180.95, W=183.84, Re=186.21, Os=190.23, Ir=192.22,
                   Pt=195.08, Au=196.97, Hg=200.59, Tl=204.38, Pb=207.2,
                   Bi=208.98, Po=209.00, At=210.00, Rn=222.00, Fr=223.00,
                   Ra=226.00, Ac=227.00, Th=232.04, Pa=231.04, U=238.03,
                   Np=237.00, Pu=244.00, Am=243.00, Cm=247.00, Bk=247.00,
                   Cf=251.00, Es=252.00, Fm=257.00, Md=258.00, No=259.00,
                   Lr=262.00, Rf=261.00, Db=262.00, Sg=266.00, Bh=264.00,
                   Hs=269.00, Mt=268.00)

# control keys in dicts 'Ge_Sb'
def key_gen(elm1,elm2):
    return '%s_%s'%(elm1,elm2)
def key_rev(string):
    return string.split('_')


# add points
def gaussian_grid(pos,Bins,gaussianSigma):
    '''
    returning gaussian-mapped grid points
    The definition of Sigma follows the one of previous in-house code.
    '''
    f = lambda x: np.exp(-(x-pos)**2/gaussianSigma) / np.sqrt(np.pi*gaussianSigma)
    return f(Bins)


def calculate_rdf(inp_file, atomnamelist, atomnumlist):
    maxDist  = 10.0  # angs
    gridSize      = float(10.0/inp_file['similarity_metric']['rdf_grid'])  # angs
    gaussianSigma = inp_file['similarity_metric']['gaussian_dist']  # angs

    inputDir = 'CONTCAR'

    PC = poscar(inputDir)
    PC.poscar()
    maxDistExtra = maxDist + 0.0

    Bins = np.arange(gridSize,maxDistExtra,gridSize)  # 0.5: extra room


    # initialize grids for RDF
    Accumul = {}
    for elm1 in PC.ElementList:
        for elm2 in PC.ElementList:
            Accumul[key_gen(elm1,elm2)] = np.zeros(np.shape(Bins))
    Accumul['tot'] = np.zeros(np.shape(Bins))


    for atomIdx in range(PC.numIon):
        NN = PC.environ(atomIdx,maxDistExtra)
        for Info in NN:
            Accumul[key_gen(PC.ElementList[PC.Label[atomIdx]],PC.ElementList[PC.Label[Info[0]]])] += gaussian_grid(Info[1],Bins,gaussianSigma)
            Accumul['tot'] += gaussian_grid(Info[1],Bins,gaussianSigma)

    # normalizing by density and 4*pi*r2 / averaging
    Density = {}
    NumIon = {}
    for key in Accumul.keys():
        if key != 'tot':
            Density[key] = PC.NumIonList[PC.ElementList.index(key_rev(key)[1])] / PC.volume
            NumIon[key] = PC.NumIonList[PC.ElementList.index(key_rev(key)[1])]
    Density['tot'] = PC.numIon/PC.volume
    NumIon['tot'] = PC.numIon

    Normaliz = {}
    for key in Accumul.keys():
        Normaliz[key] = Accumul[key] / (Density[key] * 4 * np.pi * Bins**2 * NumIon[key])

    # write output
    atomdic = {}
    for atom in range(len(atomnamelist)):
        atomdic[atomnamelist[atom]] = atomnumlist[atom] 


    weights = []
    total_weight = 0.0
    for i in range(len(atomnamelist)):
        for j in range(len(atomnamelist)):
            total_weight += float(atomnumlist[i]*atomnumlist[j])
            weights.append(float(atomnumlist[i]*atomnumlist[j]))

    for i in range(len(weights)):
        weights[i] = float(weights[i]/total_weight)

    weightnum = -1
    rdfs = []

    for i in range(len(atomnamelist)):
        for j in range(len(atomnamelist)):

            weightnum += 1

            atom1 = atomnamelist[i]
            atom2 = atomnamelist[j]

            key = atom1 + "_" + atom2

            for k in range(len(Normaliz[key])):

                rdf_value = Normaliz[key][k] * np.sqrt(weights[weightnum]) - 1.0

                rdfs.append(rdf_value)

    return rdfs


def write_contcar(coor, latt, atomnumlist, atomnamelist, tot_atom_num):
    contcar = "contcar for rdf\n"
    contcar += "1.000\n"
    contcar += str(latt[0]) + " 0.0 0.0\n"
    contcar += str(latt[3]) + " " + str(latt[1]) + " 0.0\n"
    contcar += str(latt[4]) + " " + str(latt[5]) + " " + str(latt[2]) + "\n"

    for j in range(len(atomnamelist)):
        contcar += atomnamelist[j] + " "
    contcar += "\n"
    for j in range(len(atomnumlist)):
        contcar += str(atomnumlist[j]) + " "
    contcar += "\n"

    contcar += "Cartesian\n"
    for j in range(tot_atom_num):
        contcar += str(coor[j][0])+" "
        contcar += str(coor[j][1])+" "
        contcar += str(coor[j][2])+" "
        contcar += "\n"
   
    with open("CONTCAR","w") as f:
        f.write(contcar) 
