import sys
import os
sys.path.append("/usr/local/xcfun_py3/lib64/python")
#sys.path.append("/usr/local/PyADF-myfork/src/")
import numpy as np
import psi4
from pkg_resources import parse_version
from molecule import Molecule

psi4.set_num_threads(16)
h2o=Molecule('H2O.xyz')
h2o.append('symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
psi4.set_options({'basis': 'aug-cc-pvdz', # can be defined later maybe?
                  'puream': 'True',
                  'DF_SCF_GUESS': 'False',
                  'scf_type': 'direct',
                  'dft_radial_scheme' : 'becke',
                  'dft_radial_points': 89,
                  'dft_spherical_points' : 770,
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8})
h2o_mol=psi4.geometry(h2o.geometry)
ene, h2o_wfn = psi4.energy('blyp', return_wfn=True)
C_h2o=np.array(h2o_wfn.Ca_subset("AO","OCC"))

# 2.1 map A and B density on grid points
from grid import GridFactoryDensity 
#we use the total system grid
Vpot = h2o_wfn.V_potential()
x, y, z, w = Vpot.get_np_xyzw()

points = np.zeros((x.shape[0],4) , dtype=np.float64)
points[: ,0] = x
points[: ,1] = y
points[: ,2] = z
points[: ,3] = w

tempA = GridFactoryDensity(h2o_mol,points,'aug-cc-pvdz',C_h2o) 
rho = 2.0*tempA.rho

#number of electrons from rho
tempA.integrate()

#rho_container = np.zeros((rho.shape[0],10),dtype=np.float_)
#rho_container[:,0] = rho
