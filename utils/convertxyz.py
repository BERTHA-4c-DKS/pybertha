#################################################################################################

class atom(object): # nuovo stile di classe python subclass di object 
    def __init__(self, at, x, y, z):
        self.__symbol = at
        self.__coordinate = (x, y, z)
        self.__charge = 0
    
    def set_symbol(self, at): 
        self.__symbol = at
    
    def get_symbol(self):
        return self.__symbol

    def set_charge(self, val): 
        self.__charge = val
    
    def get_charge(self):
        return self.__charge
 
    def set_coordinates(self, x, y, z):
        self.__coordinate = (x, y, z)
 
    def get_coordinates(self):
        return self.__coordinate
 
    def get_str(self):
        return '%3s %10.4f %10.4f %10.4f' % (self.__symbol, self.__coordinate[0], \
          self.__coordinate[1], self.__coordinate[2])
 
    def __repr__(self): # overloads printing
        return self.get_str()

#################################################################################################

class molecule(object):
    def __init__ (self, nome = "noname"):
        self.__name = nome
        self.__list_atoms = []
 
    def add_atom (self, atom):
        self.__list_atoms.append(atom)
 
    def get_atoms (self):
        return self.__list_atoms

    def get_num_of_atoms(self):
        return len(self.__list_atoms)
 
    def __repr__ (self):
        str = 'Molecule %s\n' % self.__name
        str = str + 'has %d atoms\n' % len(self.__list_atoms)
      
        for atom in self.__list_atoms:
            str = str + atom.get_str() + '\n'
 
        return str

#################################################################################################

import sys

mol = molecule()

converter = 1.0
xyzfilename = ""

if len(sys.argv) != 3:
    print("usage: " + sys.argv[0] + " xyzfilename coverter")
    exit(1)
else:
    converter = float(sys.argv[2])
    xyzfilename = sys.argv[1]

with open(xyzfilename) as fp:
    dim = int(fp.readline())
    header  = fp.readline()
    for i in range(dim):
        l = fp.readline()
        sl = l.split()

        if len(sl) != 4 and len(sl) != 5:
            print("Error at line "+ l)
            exit(1)

        a = atom (sl[0], float(sl[1]) * converter, \
            float(sl[2]) * converter, float(sl[3]) * converter)

        if len(sl) == 5:
            a.set_charge(int(sl[4])) 

        mol.add_atom(a)
        
print(mol)
