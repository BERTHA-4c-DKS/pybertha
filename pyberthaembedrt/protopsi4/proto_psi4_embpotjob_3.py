import sys
import os
#sys.path.append("/usr/local/xcfun_py3/lib64/python")
#sys.path.append("/usr/local/PyADF-myfork/src/")
sys.path.append('/home/mdesantis/local/xcfun_py3/lib64/python')
import numpy as np
import psi4
from pkg_resources import parse_version
from molecule import Molecule
import pyadf 
import pyadf.PyEmbed

from pyadf.Plot.FileWriters import GridWriter, GridFunctionWriter
from pyadf.Plot.FileReaders import GridFunctionReader
from pyadf.Plot.GridFunctions import GridFunctionFactory
from pyadf.Plot.GridFunctions import GridFunctionContainer
import math

#adf_settings = pyadf.adfsettings()
#adf_settings.set_save_tapes([21,10])
#adf_settings.set_functional('BLYP')
#adf_settings.set_convergence(1.0e-8)
#adf_settings.set_integration(accint=4.0)

# other options
#basis_active = "DZP"
#fde_act_opts = {'FULLGRID':'', 'PW91k' : '', 'ENERGY' : '' }

#file_h2o = os.path.join ("./", 'H2O.xyz')
#file_nh3 = os.path.join ("./", 'NH3.xyz')

m_dummy = pyadf.molecule('NH3.xyz')
m_dummy.set_symmetry('NOSYM')
#m_nh3 = pyadf.molecule(file_nh3)
#m_nh3.set_symmetry('NOSYM')

#m_tot = m_h2o + m_nh3
#m_tot.set_symmetry('NOSYM')

h2o=Molecule('H2O.xyz')
nh3=Molecule('NH3.xyz')

tot=Molecule('H2O.xyz')
tot.append(nh3.geometry+'symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
#check
tot.display_xyz()

psi4.set_num_threads(16)
psi4.core.set_output_file('output.dat', False)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Grid generation step
#
# getting a grid from psi4 for the total system

if parse_version(psi4.__version__) >= parse_version('1.3a1'):
    build_superfunctional = psi4.driver.dft.build_superfunctional
else:
    build_superfunctional = psi4.driver.dft_funcs.build_superfunctional


tot_mol=psi4.geometry(tot.geometry)

psi4.set_options({'basis': 'sto-3g', # can be defined later maybe?
                  'puream': 'True',
                  'DF_SCF_GUESS': 'False',
                  'scf_type': 'direct',
                  'dft_radial_scheme' : 'becke',
                  'dft_radial_points': 49,
                  'dft_spherical_points' : 170,
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8})

basis_tot = psi4.core.BasisSet.build(tot_mol, "ORBITAL", "AUG-CC-PVDZ")
sup = build_superfunctional("BLYP", True)[0]
Vpot = psi4.core.VBase.build(basis_tot, sup, "RV")
Vpot.initialize()
x, y, z, w = Vpot.get_np_xyzw()
Vpot.finalize()

points = np.zeros((x.shape[0],4) , dtype=np.float64)
points[: ,0] = x
points[: ,1] = y
points[: ,2] = z
points[: ,3] = w

psi4.core.clean()



### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# FDE calculation
#

# 1. isolated subsystem calculations, get densities and nuclear potentials
h2o.append('symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
h2o_mol=psi4.geometry(h2o.geometry)
ene, h2o_wfn = psi4.energy('blyp', return_wfn=True)
C_h2o=np.array(h2o_wfn.Ca_subset("AO","OCC"))

psi4.core.clean()

nh3.append('symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
nh3_mol=psi4.geometry(nh3.geometry)
ene, nh3_wfn = psi4.energy('blyp', return_wfn=True)

C_nh3=np.array(nh3_wfn.Ca_subset("AO","OCC"))

import psi4.core as p4c
psi4_matrix = p4c.Matrix.from_array(points[:,:3])
nh3_epc = p4c.ESPPropCalc(nh3_wfn)

elpot_nh3 = np.array(nh3_epc.compute_esp_over_grid_in_memory( psi4_matrix ))

#custom grid object (pyadf)
cgrid=pyadf.customgrid(mol=m_dummy,coords=points[:,:3],weights=w)
#DEBUG: dump the psi4 grid
GridWriter.write_xyzw(grid=cgrid,filename=os.path.join("./", 'PSI4GRID'),add_comment=False)


# 2.1 map A and B density on grid points
from grid import GridFactoryDensity 
#we use the total system grid
#the first positional arguents (psi4 molecule object) is used to generate the basis set.
tempA = GridFactoryDensity(h2o_mol,points,'sto-3g',C_h2o) 
rho = 2.0*tempA.rho
rho_container = np.zeros((rho.shape[0],10),dtype=np.float_)
rho_container[:,0] = rho
#DEBUG
print("rhoA length: %i\n" % (tempA.rho).shape[0])
#quick check (n. electrons)
tempA.integrate()

#set the pyadf container
dens_gf_act = GridFunctionFactory.newGridFunction(cgrid,np.ascontiguousarray(rho_container[:,0]), gf_type="density") 
densgrad = GridFunctionFactory.newGridFunction(cgrid, np.ascontiguousarray(rho_container[:, 1:4]))   #for now we have no gradient and hessian of the density
denshess = GridFunctionFactory.newGridFunction(cgrid, np.ascontiguousarray(rho_container[:, 4:10]))  
density_active = GridFunctionContainer([dens_gf_act, densgrad, denshess])
print("check dens_active type %s\n" % type(density_active))
#B

tempB = GridFactoryDensity(nh3_mol,points,'sto-3g',C_nh3) 
rho = 2.0*tempB.rho
#we re-use the previously defined rho_container
rho_container[:,0] = rho
#DEBUG
print("rhoB length: %i\n" % (tempA.rho).shape[0])

dens_gf_enviro = GridFunctionFactory.newGridFunction(cgrid,np.ascontiguousarray(rho_container[:,0]), gf_type="density") 
isolated_dens_enviro =  GridFunctionContainer([dens_gf_enviro, densgrad, denshess])
#check the pyadf object
nel_enviro = isolated_dens_enviro[0].integral()
print("  Integrated number of electrons for subsystem B: %.8f" % nel_enviro)

#the pyadf object for elpot_nh3
#DEBUG

print("elpot_nh3 length: %i\n" % elpot_nh3.shape[0])

elpot_gf_nh3 = GridFunctionFactory.newGridFunction(cgrid,np.ascontiguousarray(elpot_nh3), gf_type="potential") 
isolated_elpot_enviro=GridFunctionContainer([elpot_gf_nh3])

# 2.2 pyembed calculations to get FDE embedding potentials from the isolated subsystems, put these in proper format

embed_settings = pyadf.PyEmbed.EmbedXCFunSettings()
embed_settings.set_fun_nad_xc ({'BeckeX': 1.0, 'LYPC': 1.0})
embed_settings.set_fun_nad_kin({'pw91k' : 1.0})

embed_settings.show_functionals()

embed_eval = pyadf.PyEmbed.EmbedXCFunEvaluator(settings=embed_settings)

nadpot_active = embed_eval.get_nad_pot(density_active,isolated_dens_enviro)
#fde_embpot_h2o = embed_eval.get_emb_pot(density_active, isolated_dens_enviro, isolated_elpot_enviro)
#fde_embpot_nh3 = embed_eval.get_emb_pot(isolated_dens_nh3, isolated_dens_h2o, isolated_elpot_h2o)

#GridFunctionWriter.write_xyzwv(fde_embpot_h2o,filename=os.path.join("./", \
#        'EMBPOT_PYEMBED_ADFGRID_H2O'),add_comment=False)
#GridFunctionWriter.write_xyzwv(fde_embpot_nh3,filename=os.path.join("./", \
#        'EMBPOT_PYEMBED_ADFGRID_NH3'),add_comment=False)
#GridFunctionWriter.write_cube(fde_embpot_h2o,filename=os.path.join("./", \
#        'EMBPOT_PYEMBED_ADFGRID_H2O.cube'))
