import numpy as np
import os
import pyadf 
import pyadf.PyEmbed

from pyadf.Plot.FileWriters import GridWriter, GridFunctionWriter
from pyadf.Plot.FileReaders import GridFunctionReader
from pyadf.Plot.GridFunctions import GridFunctionFactory
from pyadf.Plot.GridFunctions import GridFunctionContainer

from molecule import Molecule
class PyEmbError(Exception):
    """Base class for exceptions in this module."""
    pass


class pyemb:

    # TODO we should use instead a generic molecular data stracuture  
    def __init__ (self, activefname = "", envirofname = "", jobtype= ""):
        self.__activefname = activefname
        self.__envirofname = envirofname
        self.__jobtype = jobtype  #'adf' or 'psi4' as provider of grid and isolated environment density
        self.grid_type = None
        self.enviro_func = None
        self.thresh_conv = None
        self.acc_int = None
        self.basis_frzn = None
        self.agrid = None
        self.isolated_dens_enviro = None
        self.isolated_elpot_enviro = None
        self.__init = False

    def set_active_fname (self, activefname):
        if not isinstance(activefname, str):
            raise TypeError("set_fittcoefffname: input must be a string")

        if (not os.path.isfile(activefname) ):
            raise PyEmbError("File  %s  does not exist" % activefname)

        self.__activefname = activefname
        self.finalize () # we should rerun the init

    def set_enviro_fname (self, envirofname):
        if not isinstance(envirofname, str):
            raise TypeError("set_fittcoefffname: input must be a string")

        if (not os.path.isfile(envirofname) ):
            raise PyEmbError("File  %s  does not exist" % envirofname)

        self.__envirofname = envirofname
        self.finalize () # we should rerun the init

    def get_active_fname (self):
        return self.__activefname

    def get_enviro_fname (self, envirofname):
        return self.__envirofname
    
    def set_options(param,gtype=1, func='BLYP', thresh=1.0e-8):
        self.grid_type = gtype
        self.enviro_func = func
        self.thresh_conv = thresh
        self.acc_int = param      # can be a list of integers or a float
       
        if self.__jobtype == 'adf' :
            adf_settings = pyadf.adfsettings()
            adf_settings.set_save_tapes([21,10])
            adf_settings.set_functional(self.enviro_func)
            adf_settings.set_convergence(self.thresh_conv)
            adf_settings.set_integration(accint=self.acc_int)

        elif self.__jobtype: 'psi4' :
            if len(self.acc_int) != 2  :
               raise PyEmbError ("input must be a list (radial and spherical points) ")

            psi4.set_options({'basis' : self.basis_frzn,
                    'puream' : 'True',
                    'dft_radial_scheme' : 'becke',
                    'dft_radial_points': , self.acc_int[0]
                    'dft_spherical_points' : self.acc_int[1],  #'dft_nuclear_scheme': 'treutler' | default
                    'scf_type' : 'direct',                     #dft_radial_ and spherical_ points determine grid size
                    'DF_SCF_GUESS': 'False',
                    'd_convergence' : self.thresh_conv,
                    'e_convergence' : self.thresh_conv})
        else:
            raise PyEmbError ("incompatible job type")
        
    def print_options(self):

        print("pyad job type : %s, grid type : %s, functional (enviro) : %s e/d thresh : %.2e\n"\
                    % (self.jobtype,self.grid_type,self.enviro_func,self.thresh_conv))
        print("grid specs (accuracy / radial & spherical points)\n")
        print(self.acc_int)
    def initialize (self):
        # Here we can insert some init code: all the basic 
        # needed for the proper get_potential start with checking
        # filename maybe we need to specify here the grid to be used ? 
        
        # TODO 

        self.__init = True

    
        
        if self.__jobype == 'adf'


            m_active = pyadf.molecule(self.__activefname)
            m_active.set_symmetry('NOSYM')
            m_enviro = pyadf.molecule(self.__envirofname)
            m_enviro.set_symmetry('NOSYM')
     
            m_tot = m_active + m_enviro
            m_tot.set_symmetry('NOSYM')
           
            #  i) use adffragmentsjob grid | deafult
            # ii) use the total sys (env+act)  grid
            #iii) use the active sys grid
 
            r_isolated_enviro = pyadf.adfsinglepointjob(m_enviro, basis_frzn, \
               settings=adf_settings, options=['NOSYMFIT']).run()
            if self.grid_type == 1:
                frags = [ pyadf.fragment(None,  [m_active]),
                        pyadf.fragment(r_isolated_enviro, [m_enviro], isfrozen=True) ]
                fde_res = pyadf.adffragmentsjob(frags, self.basis_frzn, settings=adf_settings, options=['NOSYMFIT'])
                fde_res=fde_res.run()
                self.agrid = pyadf.adfgrid(fde_res)
            elif self.grid_type == 2:
                r_tot = pyadf.adfsinglepointjob(m_tot, self.basis_frzn, settings=adf_settings, options=['NOSYMFIT']).run()
                self.agrid = pyadf.adfgrid(r_tot)
            elif self.grid_type == 3:
                r_act = pyadf.adfsinglepointjob(m_active, self.basis_frzn, settings=adf_settings, options=['NOSYMFIT']).run()
                self.agrid = pyadf.adfgrid(r_act)
            
            #elif self.grid_type == 4:
            #    #override 
            #    adf_settings.set_functional("BLYP")
            #    frags = [ pyadf.fragment(None,  [m_active]),
            #            pyadf.fragment(r_isolated_enviro, [m_enviro], isfrozen=True) ]
            #    fde_res = pyadf.adffragmentsjob(frags, basis="AUG/ADZP", settings=adf_settings, fde=fde_act_opts, options=['NOSYMFIT\n EXCITATIONS\n  ONLYSING\n  LOWEST 2\nEND'])
            #    fde_res=fde_res.run()
            #    agrid = pyadf.adfgrid(fde_res) 
            
            self.isolated_dens_enviro  = r_isolated_enviro.get_density(grid=self.agrid, \
                          fit=False, order=2) #grid=grid_active=agrid
            isolated_vnuc_enviro  = r_isolated_enviro.get_potential(grid=self.agrid,\
              pot='nuc')
            isolated_coul_enviro  = r_isolated_enviro.get_potential(grid=self.agrid,\
              pot='coul')
            self.isolated_elpot_enviro = isolated_vnuc_enviro + isolated_coul_enviro
      


        #temporary placeholder
        elif self.__jobtype == 'psi4'
            #dummy pyadf molecule
            m_dummy = pyadf.molecule(self.__envirofname)
            
            #strings for psi4 molecule obj
            m_active=Molecule(self.__activefname)
            m.active.append('symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
            
            m_enviro=Molecule(self.__envirofname)

            tot=Molecule(self.__activefname)
            tot.append(m_enviro.geometry+'symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
            if grid_type == 3:

               raise PyEmbError ("incompatible grid type")
            #psi4 block
             
            if parse_version(psi4.__version__) >= parse_version('1.3a1'):      
                build_superfunctional = psi4.driver.dft.build_superfunctional
            else:
                build_superfunctional = psi4.driver.dft_funcs.build_superfunctional  
            
            tot_mol=psi4.geometry(tot.geometry)
            if self.grid_type == 1
            else :
            
            basis_dummy = psi4.core.BasisSet.build(mol_obj, "ORBITAL", self.basis_frzn)

            sup = build_superfunctional(args.frzn_func, True)[0]

            Vpot = psi4.core.VBase.build(basis_dummy, sup, "RV")
            Vpot.initialize()

            x, y, z, w = Vpot.get_np_xyzw()
            Vpot.finalize()
         
            points = np.zeros((x.shape[0],4)) #dtype?
            points[: ,0] = x
            points[: ,1] = y
            points[: ,2] = z
            points[: ,3] = w
         
            psi4.core.clean()
        
            # prepare a custom pyadf grid object out of x,y,z,w. A mol object is only needed for cube dumping (see documentation)
            self.agrid=pyadf.customgrid(mol=m_dummy,coords=np.ascontiguousarray(points[:,:3],dtype=np.float_),weights=np.ascontiguousarray(w,dtype=np.float_))
            # TODO : code for the calculation of ESP and isolated enviroment density
        else: 
            raise PyEmbError ("incompatible job type")

    def get_grid(self):
        npoints = self.agrid.npoints
        grid = numpy.zeros((npoints, 4))
        grid[:,:3] = self.agrid.get_coordinates(bohr=True)
        grid[:,3] = self.agrid.get_weights()
        grid = numpy.ascontiguousarray(grid, dtype=numpy.double)
       
        return grid

    def get_potential (self, density):
        # TODO given the density on the grid get the potential 
        # remind : density labels the active system

        if not isinstance(density,(np.ndarray)):
            raise TypeError("input must be a numpy.ndarray")

        
        npoints = self.agrid.npoints
        
        if len(density.shape) == 2 and density.shape[1] == 10 :

                if (npoints != density.shape[0]):
                    raise PyEmbError ("incomaptible grid dimension")
        else:
            raise TypeError("input must be a numpy.ndarray npoints X 10")
                
        density_gf = GridFunctionFactory.newGridFunction(self.agrid,np.ascontiguousarray(density[:,0]), gf_type="density")
        densgrad = GridFunctionFactory.newGridFunction(self.agrid, np.ascontiguousarray(density[:, 1:4]))
        denshess = GridFunctionFactory.newGridFunction(self.agrid, np.ascontiguousarray(density[:, 4:10]))  
        
        #wrapper container
        density_act = GridFunctionContainer([density_gf, densgrad, denshess])

        nadpot_active=embed_eval.get_nad_pot(density_act, self.isolated_dens_enviro)
        embpot = self.isolated_elpot_enviro + nadpot_active

        #get the potential as numpy.ndarray
        pot = embpot.get_values()
        pot = np.ascontiguousarray(pot, dtype=np.double)
        
        # TODO 

        return pot

    def finalize (self):
        self.__init = False
        self.__activefname = ""
        self.__envirofname = ""

        # TODO clean all the allocated data 
