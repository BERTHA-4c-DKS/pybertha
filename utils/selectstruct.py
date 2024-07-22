from rdkit import Chem
from rdkit.Chem.AllChem import GetBestRMS, AlignMol
from rdkit.Chem import Draw
import sys

import matplotlib.pyplot as plt

#########################################################################3

def read_xyzblocks (coord_filepath):

    fp = open(coord_filepath, "r")
    readingmol = False
    readfisrtline = False
    atomread = 0
    natoms = 0
    molnum = 0
    xyzblocks = []
    for line in fp:
        # remove all mutiple spaces in the line
        nsline = ' '.join(line.split())
        snsline = nsline.split()
        if len(snsline) == 1:
            try :
                natoms = int(snsline[0])
                readingmol = True
                readfisrtline = True
                atomread = 0
                xyzblocks.append(line)
            except:
                print("Error: number of atoms is not an integer")
                exit()
        else:
            if readingmol:
                # readfist line of molecule
                if readfisrtline:
                    readfisrtline = False
                    xyzblocks[-1] += line
                else:
        
                    if len(snsline) == 4:
                        atom = snsline[0]
                        x = float(snsline[1])
                        y = float(snsline[2])
                        z = float(snsline[3])
                        #print(atom, x, y, z)
                        xyzblocks[-1] += line
                        atomread += 1
                    
                    if atomread == natoms:
                        readingmol = False
                        molnum += 1
                        #print("End of molecule ", molnum)

    return xyzblocks

#########################################################################

def cmp_rmsd (atoms1, atoms2, checkonly=["C", "H", "O", "Pb"]):

    # get atoms of the same type
    atoms1 = sorted(atoms1, key=lambda x: x[0])
    atoms2 = sorted(atoms2, key=lambda x: x[0])

    totrmsd = 0.0
    for atomtocheck in checkonly:
        as1 = [atom for atom in atoms1 if atom[0] == atomtocheck]
        as2 = [atom for atom in atoms2 if atom[0] == atomtocheck]
        if len(as1) != len(as2):
            print("Error: number of atoms of type ", atomtocheck, " is different")
            exit()
        else:
            for a1 in as1:
                # get nearest atom in as2
                minrmsd = float("inf")
                for a2 in as2:
                    # compute distance between a1 and a2
                    dx = a1[1] - a2[1]
                    dy = a1[2] - a2[2]
                    dz = a1[3] - a2[3]
                    rmsd = (dx**2 + dy**2 + dz**2)**0.5
                    if rmsd < minrmsd:
                        minrmsd = rmsd
                totrmsd += minrmsd

    return totrmsd

#########################################################################

if __name__ == "__main__":
    userdkit = False

    coord1_filepath = "1.xyz"
    coord2_filepath = "1.xyz"

    if len(sys.argv) == 3:
        coord1_filepath = sys.argv[1]
        coord2_filepath = sys.argv[2]

    xyzblocks1 = read_xyzblocks(coord1_filepath)
    xyzblocks2 = read_xyzblocks(coord2_filepath)

    print("Number of molecules in file 1: ", len(xyzblocks1))
    print("Number of molecules in file 2: ", len(xyzblocks2))

    if userdkit:
        mols1 = []
        for xyzblock in xyzblocks1:
            #print(xyzblock)
            mol = Chem.MolFromXYZBlock(xyzblock)
            conn_mol = Chem.Mol(mol)
            mols1.append(conn_mol)
            #img = Draw.MolToImage(conn_mol)
            #fig, ax = plt.subplots()
            #ax.imshow(img)
            #ax.grid(False)
            #ax.axis('off')
            #plt.show()
        
        mols2 = []
        for xyzblock in xyzblocks2:
            #print(xyzblock)
            mol = Chem.MolFromXYZBlock(xyzblock)
            conn_mol = Chem.Mol(mol)
            mols2.append(conn_mol)    
            #img = Draw.MolToImage(conn_mol)
            #fig, ax = plt.subplots()
            #ax.imshow(img)
            #ax.grid(False)
            #ax.axis('off')
            #plt.show()
      
        for i in range(len(mols1)):
            for j in range(len(mols2)):
    
                mol1 = mols1[i]
                mol2 = mols2[j]
                rmsd = GetBestRMS(mol1, mol2)
                print("RMSD: ", rmsd)
                rmsd = AlignMol(mol1, mol2)
                print("RMSD: ", rmsd)

    else:
        similmols1 = []
        for i in range(len(xyzblocks1)):
            xyzblock1 = xyzblocks1[i]
            line1 = xyzblock1.split("\n")[0]
            line1 = ' '.join(line1.split())
            natoms = int(line1.split()[0])
            atoms = []
            for atomidx in range(natoms):
                line = xyzblock1.split("\n")[atomidx+2]
                line = ' '.join(line.split())
                atom = line.split()[0]
                x = float(line.split()[1])
                y = float(line.split()[2])
                z = float(line.split()[3])
                atoms.append([atom, x, y, z])
            similmols1.append(atoms)
        
        similmols2 = []
        for i in range(len(xyzblocks2)):
            xyzblock2 = xyzblocks2[i]
            line2 = xyzblock2.split("\n")[0]
            line2 = ' '.join(line2.split())
            natoms = int(line2.split()[0])
            atoms = []
            for atomidx in range(natoms):
                line = xyzblock2.split("\n")[atomidx+2]
                line = ' '.join(line.split())
                atom = line.split()[0]
                x = float(line.split()[1])
                y = float(line.split()[2])
                z = float(line.split()[3])
                atoms.append([atom, x, y, z])
            similmols2.append(atoms)


        minrmsd = float("inf")
        minidx = [-1, -1]
        for i in range(len(similmols1)):
            for j in range(len(similmols2)):
                atoms1 = similmols1[i]
                atoms2 = similmols2[j]
                rmsd = cmp_rmsd(atoms1, atoms2)
                print("%d %d RMSD: "%(i, j), rmsd)
                if rmsd < minrmsd:
                    minrmsd = rmsd
                    minidx = [i, j]
      
    