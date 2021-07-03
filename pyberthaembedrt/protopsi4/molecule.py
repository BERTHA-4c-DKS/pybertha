class Molecule():

  def __init__(self,fgeom='geom.xyz'):

      self.geometry = None
      self.gdummy = None
      self.set_geometry(fgeom)
      
  def set_geometry(self,fgeom):

      self.geometry=str()
      with open(fgeom,'r') as data:
         self.natom = int(data.readline()) # natom is the 'proper' number of atoms in the (active) molecule
         next(data)
         for line in data:
            self.geometry += str(line)
      self.internal = self.geometry        # store a 'copy' of the isolated fragment geometry

  def set_ghost(self):
      self.gdummy=str()
      tmp=self.geometry.split('\n')
      tmp.pop()
      for m in tmp:
        self.gdummy +="@"+m.strip()+'\n'

  def append(self,options):                # options is a string like "symmetry c1"+"\n" or a string containing moelcular coordinates
      self.geometry += options
  
  def display_xyz(self):
      print(self.geometry)
