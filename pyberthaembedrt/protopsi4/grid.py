import psi4
import numpy as np



class GridFactoryDensity():
  def __init__(self,mol,points,basis_set,Ca):
      
      self.mol = mol
      self.points = points
      self.basisset = basis_set
      self.Ca = Ca
      self.rho = None
      self.phi = None
      self.phi_builder()  
      self.make_rho()
      
  def phi_builder(self):
      
      xs=psi4.core.Vector.from_array(self.points[:,0])
      ys=psi4.core.Vector.from_array(self.points[:,1])
      zs=psi4.core.Vector.from_array(self.points[:,2])
      ws=psi4.core.Vector.from_array(self.points[:,3])

      delta = 1.0e-2 #private parameter 
 
      basis = psi4.core.BasisSet.build(self.mol, 'ORBITAL',self.basisset)
      basis_extents = psi4.core.BasisExtents(basis,delta)
 
      blockopoints = psi4.core.BlockOPoints(xs, ys, zs, ws,basis_extents)
      npoints = blockopoints.npoints()
      print("n points: %i" % npoints)
 
      nbas = basis.nbf() #number of basis functions
 
      funcs = psi4.core.BasisFunctions(basis,npoints,nbas)
 
      funcs.compute_functions(blockopoints)
 
      phi = np.array(funcs.basis_values()["PHI"])[:npoints, :]
      self.phi = phi

  def make_rho(self):
      MO = np.matmul(self.phi,self.Ca) 
      MO_dens = np.square(MO)
      self.rho = np.einsum('pm->p',MO_dens)
  
  def integrate(self):
      n_el=np.einsum('a,a->',self.rho,self.points[:,3])
      print("n. electrons: %.8f\n" % (2.0*n_el))
