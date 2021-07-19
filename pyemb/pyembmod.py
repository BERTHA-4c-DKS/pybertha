import numpy 
import os
import pyadf 
import pyadf.PyEmbed

from pyadf.Plot.FileWriters import GridWriter, GridFunctionWriter
from pyadf.Plot.FileReaders import GridFunctionReader
from pyadf.Plot.GridFunctions import GridFunctionFactory
from pyadf.Plot.GridFunctions import GridFunctionContainer

class Molecule():

  def __init__(self,fgeom='geom.xyz'):

      self.geometry = None
      self.gdummy = None
      self.set_geometry(fgeom)
      
  def set_geometry(self,fgeom):

      self.geometry=str()
      with open(fgeom,'r') as data:
         self.natom = int(data.readline()) # natom is the 'proper' 
                                           # number of atoms in the 
                                           # (active) molecule
         next(data)
         for line in data:
            self.geometry += str(line)
      self.internal = self.geometry        # store a 'copy' of the isolated 
                                           # fragment geometry

  def set_ghost(self):
      self.gdummy=str()
      tmp=self.geometry.split('\n')
      tmp.pop()
      for m in tmp:
        self.gdummy +="@"+m.strip()+'\n'

  def append(self,options):                
      # options is a string like "symmetry c1"+"\n" or a 
      # string containing moelcular coordinates
      self.geometry += options
  
  def display_xyz(self):
      print(self.geometry)


class PyEmbError(Exception):
    """Base class for exceptions in this module."""
    pass


class pyemb:

    # TODO we should use instead a generic molecular data stracuture  
    def __init__ (self, activefname = "", envirofname = "", jobtype= ""):
        self.__activefname = activefname
        self.__envirofname = envirofname
        self.__jobtype = jobtype  #'adf' or 'psi4' as provider of grid and isolated environment density
        self.__agrid = None
        self.__isolated_dens_enviro = None
        self.__isolated_elpot_enviro = None
        self.__init = False
        self.__acc_int = None

        self.__grid_type = 1
        self.__enviro_func = "BLYP"
        self.__basis_frzn = 'AUG/ADZP'
        self.__thresh_conv = 1.0e-8
        # options for xcfun 
        self.__f_nad_xc = 'lda'
        self.__f_nad_kin = 'tkf'

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
    
    def set_grid_type(self, gtype):

        if not isinstance(gtype, int):
            raise TypeError("input must be an integer")

        self.__grid_type = gtype

    def get_grid_type(self):

        return self.__grid_type

    def set_enviro_func(self, func):

        if not isinstance(func, str):
            raise TypeError("input must be a string")

        self.__enviro_func = func

    def get_enviro_func(self):
        
        return self.__enviro_func

    def set_basis_frzn(self, basis):

        if not isinstance(basis, str):
            raise TypeError("input must be a string")

        self.__basis_frzn = basis

    def get_basis_frzn(self):
        
        return self.__basis_frzn

    def set_thresh_conv(self, thresh):

        if not isinstance(thresh, float):
            raise TypeError("input must be a floating-point")

        self.__thresh_conv = thresh

    def get_thresh_conv(self):
        
        return self.__thresh_conv

    def set_f_nad_xc(self, xc):

        if not isinstance(xc, str):
            raise TypeError("input must be a string")

        self.__f_nad_xc = xc

    def get_f_nad_xc(self):
        
        return self.__f_nad_xc

    def set_f_nad_kin(self, kin):

        if not isinstance(kin, str):
            raise TypeError("input must be a string")

        self.__f_nad_kin = kin

    def get_f_nad_kin(self):
        
        return self.__f_nad_kin

    def set_options(self, param, gtype=1, func='BLYP', basis ='AUG/ADZP', \
        f_nad={'xc' : 'lda', 'kin' : 'tfk'}, thresh=1.0e-8):

        self.__acc_int = param      # can be a list of integers or a float

        self.set_grid_type(gtype)
        self.set_enviro_func(func)
        self.set_basis_frzn(basis)
        self.set_thresh_conv (thresh)

        # options for xcfun 
        if not isinstance(f_nad, dict):
            raise TypeError("f_nad (parameter 5) must be a dictionary")

        self.set_f_nad_xc(f_nad['xc'])
        self.set_f_nad_kin(f_nad['kin'])
        
    def get_options(self):

        outstr = "pyad job type : %s, grid type : %s, functional (enviro) : %s, e/d thresh : %.2e\n" \
                    % (self.__jobtype,self.__grid_type,self.__enviro_func,self.__thresh_conv)
        outstr += "grid specs (accuracy / radial & spherical points)\n"
        outstr += str(self.__acc_int)

        return outstr

    def initialize (self):
        # Here we can insert some init code: all the basic 
        # needed for the proper get_potential start with checking
        # filename maybe we need to specify here the grid to be used ? 

        if self.__jobtype == 'adf':

            adf_settings = pyadf.adfsettings()
            adf_settings.set_save_tapes([21,10])
            adf_settings.set_functional(self.__enviro_func)
            adf_settings.set_convergence(self.__thresh_conv)
            adf_settings.set_integration(accint=self.__acc_int)

            m_active = pyadf.molecule(self.__activefname)
            m_active.set_symmetry('NOSYM')
            m_enviro = pyadf.molecule(self.__envirofname)
            m_enviro.set_symmetry('NOSYM')
     
            m_tot = m_active + m_enviro
            m_tot.set_symmetry('NOSYM')
           
            #  i) use adffragmentsjob grid | deafult
            # ii) use the total sys (env+act)  grid
            #iii) use the active sys grid
 
            r_isolated_enviro = pyadf.adfsinglepointjob(m_enviro, self.__basis_frzn, \
               settings=adf_settings, options=['NOSYMFIT']).run()
            if self.__grid_type == 1:
                frags = [ pyadf.fragment(None,  [m_active]),
                        pyadf.fragment(r_isolated_enviro, [m_enviro], isfrozen=True) ]
                fde_res = pyadf.adffragmentsjob(frags, self.__basis_frzn, settings=adf_settings, options=['NOSYMFIT'])
                fde_res=fde_res.run()
                self.__agrid = pyadf.adfgrid(fde_res)
            elif self.__grid_type == 2:
                r_tot = pyadf.adfsinglepointjob(m_tot, self.__basis_frzn, settings=adf_settings, options=['NOSYMFIT']).run()
                self.__agrid = pyadf.adfgrid(r_tot)
            elif self.__grid_type == 3:
                r_act = pyadf.adfsinglepointjob(m_active, self.__basis_frzn, settings=adf_settings, options=['NOSYMFIT']).run()
                self.__agrid = pyadf.adfgrid(r_act)
            
            #elif self.__grid_type == 4:
            #    #override 
            #    adf_settings.set_functional("BLYP")
            #    frags = [ pyadf.fragment(None,  [m_active]),
            #            pyadf.fragment(r_isolated_enviro, [m_enviro], isfrozen=True) ]
            #    fde_res = pyadf.adffragmentsjob(frags, basis="AUG/ADZP", settings=adf_settings, \
            #          fde=fde_act_opts, options=['NOSYMFIT\n EXCITATIONS\n  ONLYSING\n  LOWEST 2\nEND'])
            #    fde_res=fde_res.run()
            #    __agrid = pyadf.adfgrid(fde_res) 
            
            self.__isolated_dens_enviro  = r_isolated_enviro.get_density(grid=self.__agrid, \
                          fit=False, order=2) #grid=grid_active=agrid
            isolated_vnuc_enviro  = r_isolated_enviro.get_potential(grid=self.__agrid,\
              pot='nuc')
            isolated_coul_enviro  = r_isolated_enviro.get_potential(grid=self.__agrid,\
              pot='coul')
            self.__isolated_elpot_enviro = isolated_vnuc_enviro + isolated_coul_enviro
      
            self.__init = True

        elif self.__jobtype == 'psi4':
            #temporary placeholder
            if len(self.__acc_int) != 2  :
               raise PyEmbError ("input must be a list (radial and spherical points) ")

            psi4.set_options({'basis' : self.__basis_frzn,
                    'puream' : 'True',
                    'dft_radial_scheme' : 'becke',
                    'dft_radial_points':  self.__acc_int[0],
                    'dft_spherical_points' : self.__acc_int[1],  #'dft_nuclear_scheme': 'treutler' | default
                    'scf_type' : 'direct',                     #dft_radial_ and spherical_ points determine grid size
                    'DF_SCF_GUESS': 'False',
                    'd_convergence' : self.__thresh_conv,
                    'e_convergence' : self.__thresh_conv})
            #dummy pyadf molecule
            m_dummy = pyadf.molecule(self.__envirofname)
            
            #strings for psi4 molecule obj
            m_active=Molecule(self.__activefname)
            m.active.append('symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
            
            m_enviro=Molecule(self.__envirofname)

            tot=Molecule(self.__activefname)
            tot.append(m_enviro.geometry+'symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
            if __grid_type == 3:

               raise PyEmbError ("incompatible grid type")
            #psi4 block
             
            if parse_version(psi4.__version__) >= parse_version('1.3a1'):      
                build_superfunctional = psi4.driver.dft.build_superfunctional
            else:
                build_superfunctional = psi4.driver.dft_funcs.build_superfunctional  
            
            tot_mol=psi4.geometry(tot.geometry)

            if self.__grid_type == 1:
               #placeholder
               npoints = 999
               points=numpy.zeros((npoints,4))
            else :
               basis_dummy = psi4.core.BasisSet.build(mol_obj, "ORBITAL", self.__basis_frzn)
             
               sup = build_superfunctional(args.frzn_func, True)[0]
             
               Vpot = psi4.core.VBase.build(basis_dummy, sup, "RV")
               Vpot.initialize()
             
               x, y, z, w = Vpot.get_np_xyzw()
               Vpot.finalize()
          
               points = numpy.zeros((x.shape[0],4)) #dtype?
               points[: ,0] = x
               points[: ,1] = y
               points[: ,2] = z
               points[: ,3] = w
          
               psi4.core.clean()
        
            # prepare a custom pyadf grid object out of x,y,z,w. A mol object is only needed for cube dumping (see documentation)
            self.__agrid=pyadf.customgrid(mol=m_dummy,coords=numpy.ascontiguousarray(points[:,:3], \
                 dtype=numpy.float_),weights=numpy.ascontiguousarray(w,dtype=numpy.float_))
            # TODO : code for the calculation of ESP and isolated enviroment density
            self.__init = True
        else: 
            raise PyEmbError ("incompatible job type")

    def get_grid(self):

        if self.__init:
            npoints = self.__agrid.npoints
            grid = numpy.zeros((npoints, 4))
            grid[:,:3] = self.__agrid.get_coordinates(bohr=True)
            grid[:,3] = self.__agrid.get_weights()
            grid = numpy.ascontiguousarray(grid, dtype=numpy.double)

            return grid
        else:
            raise PyEmbError ("Need to be initialized")
       

    def get_potential (self, density):

        if self.__init:
            if not isinstance(density,(numpy.ndarray)):
                raise TypeError("input must be a numpy.ndarray")
            
            npoints = self.__agrid.npoints
            
            if len(density.shape) == 2 and density.shape[1] == 10 :
        
                    if (npoints != density.shape[0]):
                        raise PyEmbError ("incomaptible grid dimension")
            else:
                raise TypeError("input must be a numpy.ndarray npoints X 10")
                    
            density_gf = GridFunctionFactory.newGridFunction(self.__agrid,numpy.ascontiguousarray(density[:,0]), gf_type="density")
            densgrad = GridFunctionFactory.newGridFunction(self.__agrid, numpy.ascontiguousarray(density[:, 1:4]))
            denshess = GridFunctionFactory.newGridFunction(self.__agrid, numpy.ascontiguousarray(density[:, 4:10]))  
            
            #wrapper container
            density_act = GridFunctionContainer([density_gf, densgrad, denshess])
            
            embed_settings = pyadf.EmbedXCFunSettings()
            embed_settings.set_fun_nad_xc ({self.__f_nad_xc  : 1.0})
            embed_settings.set_fun_nad_kin({self.__f_nad_kin : 1.0})
            
            embed_settings.show_functionals()
            
            embed_eval = pyadf.PyEmbed.EmbedXCFunEvaluator(settings=embed_settings)
        
            nadpot_active=embed_eval.get_nad_pot(density_act, self.__isolated_dens_enviro)
            embpot = self.__isolated_elpot_enviro + nadpot_active
        
            #get the potential as numpy.ndarray
            pot = embpot.get_values()
            pot = numpy.ascontiguousarray(pot, dtype=numpy.double)
            
            return pot
        else:
            raise PyEmbError ("Need to be initialized")

    def finalize (self):

        self.__activefname = activefname
        self.__envirofname = envirofname
        self.__jobtype = jobtype  #'adf' or 'psi4' as provider of grid and isolated environment density
        self.__agrid = None
        self.__isolated_dens_enviro = None
        self.__isolated_elpot_enviro = None
        self.__init = False
        self.__acc_int = None

        self.__grid_type = 1
        self.__enviro_func = "BLYP"
        self.__basis_frzn = 'AUG/ADZP'
        self.__thresh_conv = 1.0e-8
        # options for xcfun 
        self.__f_nad_xc = 'lda'
        self.__f_nad_kin = 'tkf'