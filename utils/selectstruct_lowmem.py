import sys
import time

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
            except:
                print("Error: number of atoms is not an integer")
                exit()
        else:
            if readingmol:
                # readfist line of molecule
                if readfisrtline:
                    readfisrtline = False
                    xyzblocks.append([])
                else:
                    if len(snsline) == 4:
                        atom = snsline[0]
                        x = float(snsline[1])
                        y = float(snsline[2])
                        z = float(snsline[3])
                        #print(atom, x, y, z)
                        xyzblocks[-1].append([atom, x, y, z])
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

    coord1_filepath = "1.xyz"
    coord2_filepath = "1.xyz"
    slicetorun = -1
    totalslices = 0
    if len(sys.argv) == 5:
        coord1_filepath = sys.argv[1]
        coord2_filepath = sys.argv[2]
        slicetorun = int(sys.argv[3])
        totalslices = int(sys.argv[4])
    elif len(sys.argv) == 3:
        coord1_filepath = sys.argv[1]
        coord2_filepath = sys.argv[2]

    similmols1 = read_xyzblocks(coord1_filepath)
    similmols2 = read_xyzblocks(coord2_filepath)

    print("Number of molecules in file 1: ", len(similmols1))
    print("Number of molecules in file 2: ", len(similmols2))
    sys.stdout.flush()

    minrmsd = float("inf")
    minidx = [-1, -1]
    startfrom = 0
    upto = len(similmols1)
    if slicetorun != -1:
        total = totalslices
        step = len(similmols1) // total
        startfrom = (slicetorun - 1)* step
        upto = startfrom + step
        if upto > len(similmols1):
            upto = len(similmols1)
        if slicetorun == total:
            upto = len(similmols1)
        print("Running slice ", slicetorun, " of ",\
               total, " from ", startfrom, " to ", upto)
        sys.stdout.flush()

    for i in range(startfrom, upto):
        localminrmsd = float("inf")
        localminidx = [-1, -1]
        # measure time needed to compare molecule i with all molecules in file 2
        start = time.time()
        for j in range(len(similmols2)):
            atoms1 = similmols1[i]
            atoms2 = similmols2[j]
            rmsd = cmp_rmsd(atoms1, atoms2)
            #print("%d %d RMSD: "%(i, j), rmsd)
            if rmsd < minrmsd:
                minrmsd = rmsd
                minidx = [i, j]

            if rmsd < localminrmsd:
                localminrmsd = rmsd
                localminidx = [i, j]

        end = time.time()
        print("Time taken to compare molecule ", i, " with all molecules in file 2: ", end-start)
        print("Molecule ", i, " in file 1 of ", \
              len(similmols1), " Min RMSD: ", minrmsd, \
                " Molecule ", minidx[1], " in file 2 of ", len(similmols2))
        print("Local Min RMSD: ", localminrmsd, " Molecule ", localminidx[1], \
              " in file 2 of ", len(similmols2))        
        sys.stdout.flush()

    print("Minimum RMSD: ", minrmsd)
    print("Molecules: ", minidx)

    # write the molecules with minimum RMSD to a file
    fp = open("minrmsd1.xyz", "w")
    fp.write(str(len(similmols1[minidx[0]])) + "\n")
    fp.write("\n")
    for atom in similmols1[minidx[0]]:
        fp.write(atom[0] + " " + str(atom[1]) + " " + str(atom[2]) + " " + str(atom[3]) + "\n")
    fp.close()
    fp = open("minrmsd2.xyz", "w")
    fp.write(str(len(similmols2[minidx[1]])) + "\n")
    fp.write("\n")
    for atom in similmols2[minidx[1]]:
        fp.write(atom[0] + " " + str(atom[1]) + " " + str(atom[2]) + " " + str(atom[3]) + "\n")
    fp.close()
