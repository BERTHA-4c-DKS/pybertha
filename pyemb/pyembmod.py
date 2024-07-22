import numpy
import scipy 
import os
import sys
import pyadf 
import pyadf.PyEmbed
import psi4
import argparse
from pkg_resources import parse_version

from pyadf.Plot.FileWriters import GridWriter, GridFunctionWriter
from pyadf.Plot.FileReaders import GridFunctionReader
from pyadf.Plot.GridFunctions import GridFunctionFactory
from pyadf.Plot.GridFunctions import GridFunctionContainer

import io
import threading

from contextlib import redirect_stdout

########################################################################################################

def drain_pipe():

  global captured_stdout
  while True:
    data = os.read(stdout_pipe[0], 1024)
    if not data:
      break
    captured_stdout.append(data.decode("utf-8"))

########################################################################################################

def init_stdout_redirect ():

    global stdout_fileno 
    global stdout_save 
    global stdout_pipe
    global captured_stdout
    
    captured_stdout = []

    stdout_fileno = sys.stdout.fileno()
    stdout_save = os.dup(stdout_fileno)
    stdout_pipe = os.pipe()
    
    os.dup2(stdout_pipe[1], stdout_fileno)
    os.close(stdout_pipe[1])

########################################################################################################

def finalize_stdout_redirect (fname, writef=1):

    global stdout_fileno 
    global stdout_save 
    global stdout_pipe 
    global captured_stdout
      
    os.close(stdout_pipe[0])
    os.dup2(stdout_save, stdout_fileno)
    os.close(stdout_save)

    if writef != 0:
        fp = None
        if writef == 1:
          fp = open(fname, "w")
        elif writef == 2:
          fp = open(fname, "a")
       
        for line in captured_stdout:
          fp.write(line)
       
        fp.close()

########################################################################################################

class Molecule():

  def __init__(self,fgeom='geom.xyz'):

      self.__geometry = None
      self.__ghosted = None
      self.__natom = None
      self.__geom_init = None
      self.set_geometry(fgeom)
      
  def set_geometry(self,fgeom):

      self.__geometry=str()
      with open(fgeom,'r') as data:
         self.__natom = int(data.readline()) # natom is the 'proper' 
                                           # number of atoms in the 
                                           # (active) molecule
         next(data)
         for line in data:
            self.__geometry += str(line)
      self.geom_init = self.__geometry        # store a 'copy' of the initialized geometry

  def set_ghost(self):
      self.__ghosted=str()
      tmp=self.__geometry.split('\n')
      tmp.pop()
      for m in tmp:
        self.__ghosted +="@"+m.strip()+'\n'

  def append(self,options):                
      # options is a string like "symmetry c1"+"\n" or a 
      # string containing moelcular coordinates
      self.__geometry += options
  
  def display_xyz(self):
      print(self.__geometry)
  
  def geom_str(self):
  
      return self.__geometry
  def geom_sghost(self):
  
      return self.__ghosted

########################################################################################################

class GridDensityFactory():

  def __init__(self,mol,points,basis_set):
      
      self.__mol = mol
      self.__points = points
      self.__basisset = basis_set
      self.__phi = None
      self.__lpos = None # back-compatibility
      self.__nbas = None #
      self.phi_builder()

  def lpos(self):
      return self.__lpos
  
  def nbf(self):
      return self.__nbas
  
  def phi(self):
      return self.__phi
  
  def points(self):
      return self.__points
      
  def phi_builder(self):
      
      xs=psi4.core.Vector.from_array(self.__points[:,0])
      ys=psi4.core.Vector.from_array(self.__points[:,1])
      zs=psi4.core.Vector.from_array(self.__points[:,2])
      ws=psi4.core.Vector.from_array(self.__points[:,3])

      delta = 1.0e-2 #private parameter 
 
      basis = psi4.core.BasisSet.build(self.__mol, 'ORBITAL',self.__basisset,puream=-1)
      basis_extents = psi4.core.BasisExtents(basis,delta)
 
      blockopoints = psi4.core.BlockOPoints(xs, ys, zs, ws,basis_extents)
      npoints = blockopoints.npoints()
      #print("n points: %i" % npoints)
      
      self.__lpos = numpy.array(blockopoints.functions_local_to_global())
      #DEBUG 
      #print("Local basis function mapping")
      #print(lpos) 
 
      self.__nbas = basis.nbf() #number of basis functions
 
      funcs = psi4.core.BasisFunctions(basis,npoints,self.__nbas)
 
      funcs.compute_functions(blockopoints)
 
      self.__phi = numpy.array(funcs.basis_values()["PHI"])[:npoints, :self.__lpos.shape[0]]

  def from_Cocc(self,Cocc):
      if (self.__phi.shape[1] != Cocc.shape[0]):
        raise("Check Cocc and PHI dim")
      MO = numpy.matmul(self.__phi,Cocc)
      MO_dens = numpy.square(MO)
      rho = numpy.einsum('pm->p',MO_dens)
      return rho

  def from_D(self,D,ovap): # to be tested
      temp=numpy.matmul(ovap,numpy.matmul(D.real,ovap))
      try:
        eigvals,eigvecs=scipy.linalg.eigh(temp,ovap,eigvals_only=False)
      except scipy.linalg.LinAlgError:
        print("Error in scipy.linalg.eigh in make_rho from D")

      idx = eigvals.argsort()[::-1]
      eigvals = eigvals[idx]
      eigvecs = eigvecs[:,idx]
      MO = numpy.matmul(self.__phi,eigvecs[:,:self.__ndocc])
      MO_dens = numpy.square(MO)
      rho = numpy.einsum('pm->p',MO_dens)
      return rho

  def print_detail(self):
      print("PHI dim: %i,%i\n" %(self.__phi.shape))
      print("N. basis functions: %i\n" %(self.__nbas))
      print("N. grid points %i\n" %(self.__phi.shape[0]))

########################################################################################################


class PyEmbError(Exception):
    """Base class for exceptions in this module."""
    pass

########################################################################################################

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
        self.__adfoufname = ""
        self.__psioufname = ""

        self.__gridfilename = "grid.dat"

    def set_grid_filename (self, name):
        self.__gridfilename = name

    def set_adf_filenameout (self, name):
        self.__adfoufname = name

    def set_psi4_filenameout (self, name):
        self.__psioufname = name

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

        #if gtype <= 0 or gtype > 3:
        if  not((gtype >= 0 and gtype <= 3) or (gtype == -99)) :
            raise TypeError("input must be an integer i : 0 <= i <= 3 or -99")

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

        if isinstance (param, float) or \
            isinstance (param, list):

            if isinstance (param, list):
                if not all(isinstance(x, int) for x in param):
                    raise TypeError("param (parameter 1) must be a float or a list of integer")

            self.__acc_int = param      # can be a list of integers or a float
        else:
            raise TypeError("param (parameter 1) must be a float or a list of integer")

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

        outstr = "pyemb job type : %s, grid type : %s, functional (enviro) : %s, e/d thresh : %.2e\n" \
                    % (self.__jobtype, self.__grid_type, self.__enviro_func, self.__thresh_conv)
        outstr += "grid specs (accuracy / radial & spherical points)\n"
        outstr += str(self.__acc_int)

        return outstr

    def initialize (self):
        # Here we can insert some init code: all the basic 
        # needed for the proper get_potential start with checking
        # filename maybe we need to specify here the grid to be used ? 

        if self.__jobtype == 'adf':

            if not isinstance (self.__acc_int, float):
              raise TypeError("param (parameter 1) must be a float if jobtype is " + self.__jobtype)
            
            if not ((self.__grid_type > 0 and self.__grid_type <= 3) or \
                     (self.__grid_type == -99)) :
              raise TypeError("gridtype must be 1,2 or 3 or -99 if jobtype is " + self.__jobtype)

            f = io.StringIO()
            with redirect_stdout(f):
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
                #PLACEHOLDER
                elif self.__grid_type == -99:
                    idx = 0

                    try:
                       with open(self.__gridfilename, "r") as fgrid:
                          nlines = int(next(fgrid))
                          points = numpy.zeros((nlines,4),dtype=numpy.float_)
                          for line in fgrid:
                             raw = line.split()
                             points[idx,:]=raw
                             idx += 1
                    except IOError:
                       raise IOError("File grid.dat does not exist")
                
                    self.__agrid=pyadf.customgrid(mol=m_tot,coords=numpy.ascontiguousarray(points[:,:3], \
                         dtype=numpy.float_),weights=numpy.ascontiguousarray(points[:,3],dtype=numpy.float_))
                #GridWriter.write_xyzw(grid=self.__agrid,filename='grid.check',add_comment=False)
                #elif self.__grid_type == 4:
                #    #override 
                #    adf_settings.set_functional("BLYP")
                #    frags = [ pyadf.fragment(None,  [m_active]),
                #            pyadf.fragment(r_isolated_enviro, [m_enviro], isfrozen=True) ]
                #    fde_res = pyadf.adffragmentsjob(frags, basis="AUG/ADZP", settings=adf_settings, \
                #          fde=fde_act_opts, options=['NOSYMFIT\n EXCITATIONS\n  ONLYSING\n  LOWEST 2\nEND'])
                #    fde_res=fde_res.run()
                #    __agrid = pyadf.adfgrid(fde_res) 
                
                self.__isolated_dens_enviro  = r_isolated_enviro.get_density(grid=self.__agrid,\
                   fit=False, order=2) #grid=grid_active=agrid
                isolated_vnuc_enviro  = r_isolated_enviro.get_potential(grid=self.__agrid,\
                  pot='nuc')
                isolated_coul_enviro  = r_isolated_enviro.get_potential(grid=self.__agrid,\
                  pot='coul')
                self.__isolated_elpot_enviro = isolated_vnuc_enviro + isolated_coul_enviro
            
            if self.__adfoufname != "":
                fp = open(self.__adfoufname, "w")
                fp.write(f.getvalue())
                fp.close()
    
            self.__init = True

        elif self.__jobtype == 'psi4':
            
            if len(self.__acc_int) != 2  :
               raise PyEmbError ("input must be a list (radial and spherical points) ")

            if self.__grid_type == 1 or self.__grid_type > 3:
              raise TypeError("gridtype must be 2 or 3 if jobtype is " + self.__jobtype)

            if not isinstance (self.__acc_int, list):
              raise TypeError("param (parameter 1) must be a list of integer if jobtype " +
               self.__jobtype )

            if not all(isinstance(x, int) for x in self.__acc_int):
              raise TypeError("param (parameter 1) must be a list of integer when jobtype " + 
                       self.__jobtype )

            init_stdout_redirect ()
            t = threading.Thread(target=drain_pipe)
            t.start()
            
            #psi4.core.set_output_file('debug_psi4.dat', False)
            psi4.set_options({'basis' : self.__basis_frzn,
                    'puream' : 'True',
                    'dft_radial_scheme' : 'becke',
                    'dft_radial_points':  self.__acc_int[0],
                    'dft_spherical_points' : self.__acc_int[1]})  #'dft_nuclear_scheme': 'treutler' | default
            #dummy pyadf molecule
            m_dummy = pyadf.molecule(self.__envirofname)
            
            #strings for psi4 molecule obj
            m_active=Molecule(self.__activefname)
            m_active.append('symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
            
            m_enviro=Molecule(self.__envirofname)

            tot=Molecule(self.__activefname)
            tot.append(m_enviro.geom_str()+'symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')

            #psi4 block
             
            if parse_version(psi4.__version__) >= parse_version('1.3a1'):      
                build_superfunctional = psi4.driver.dft.build_superfunctional
            else:
                build_superfunctional = psi4.driver.dft_funcs.build_superfunctional  
            
            if self.__grid_type == 2:
               mol_obj=psi4.geometry(tot.geom_str())
            elif self.__grid_type == 3:
               mol_obj=psi4.geometry(m_active.geom_str())
            elif self.__grid_type == -99:
               print("Grid from extern file\n")
               

            else:
               raise PyEmbError ("incompatible grid type")

            if self.__grid_type == -99:
               idx = 0
               with open("grid.dat","r") as fgrid:
                nlines = int(next(fgrid))
                gridbin = numpy.zeros((nlines,4),dtype=numpy.float_)
                for line in fgrid:
                 raw = line.split()
                 gridbin[idx,:]=raw
                 idx += 1

               fgrid.close() 
               x = gridbin[:,0]
               y = gridbin[:,1]
               z = gridbin[:,2]
               w = gridbin[:,3]

            else:
               basis_dummy = psi4.core.BasisSet.build(mol_obj, "ORBITAL", self.__basis_frzn)
               
               sup = build_superfunctional(self.__enviro_func, True)[0]
               
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
            GridWriter.write_xyzw(grid=self.__agrid,filename='grid.txt',add_comment=False)
            # TODO : code for the calculation of ESP and isolated enviroment density
            #the quality of the grid for the environment gs calculation is set to 'good' (see Psi4 man)
            psi4.set_options({'dft_radial_scheme' : 'becke',    #dft_radial_ and spherical_ points determine grid size
                              'dft_radial_points':  75,
                              'dft_spherical_points' : 434, #'dft_nuclear_scheme': 'treutler' | default
                              'puream' : 'True',            
                              'scf_type' : 'direct',                     
                              'DF_SCF_GUESS': 'False',
                              'd_convergence' : self.__thresh_conv,
                              'e_convergence' : self.__thresh_conv})
            m_enviro.append('symmetry c1' +'\n' + 'no_com' + '\n' + 'no_reorient' + '\n')
            enviro_obj=psi4.geometry(m_enviro.geom_str())
            ene, enviro_wfn = psi4.energy(self.__enviro_func,return_wfn=True)
            
            psi4_matrix = psi4.core.Matrix.from_array(points[:,:3])
            enviro_epc = psi4.core.ESPPropCalc(enviro_wfn)

            elpot_enviro = numpy.array(enviro_epc.compute_esp_over_grid_in_memory( psi4_matrix ))
            psi4.core.clean()

           
            #represent the environ density on the grid
            environment=GridDensityFactory(enviro_obj,points,self.__basis_frzn)
            environment.print_detail()
            enviro_dens = environment.from_Cocc(numpy.asarray(enviro_wfn.Ca_subset("AO","OCC") ))
            density=numpy.zeros((x.shape[0],10))
            density[:,0] = 2.0*enviro_dens

            #cast into pyadf containers
            dens_gf_enviro = GridFunctionFactory.newGridFunction(self.__agrid,numpy.ascontiguousarray(density[:,0],dtype=numpy.float_),gf_type="density")
            densgrad = GridFunctionFactory.newGridFunction(self.__agrid, numpy.ascontiguousarray(density[:, 1:4],dtype=numpy.float_))
            denshess = GridFunctionFactory.newGridFunction(self.__agrid, numpy.ascontiguousarray(density[:, 4:10],dtype=numpy.float_))  
            self.__isolated_dens_enviro =  GridFunctionContainer([dens_gf_enviro, densgrad, denshess])

            self.__isolated_elpot_enviro = GridFunctionFactory.newGridFunction(self.__agrid,numpy.ascontiguousarray(elpot_enviro,dtype=numpy.float_), gf_type="potential") 

            os.close(stdout_fileno)
            t.join()
            
            towrite = 1
            if self.__psioufname == "":
                towrite = 0

            finalize_stdout_redirect(self.__psioufname, towrite)
 
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

            pot = None

            f = io.StringIO()
            with redirect_stdout(f):
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

            if self.__adfoufname != "":
                fp = open(self.__adfoufname, "w")
                fp.write(f.getvalue())
                fp.close()
 
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
        self.__adfoufname = ""
        self.__psioufname = ""

        self.__gridfilename = "grid.dat"
