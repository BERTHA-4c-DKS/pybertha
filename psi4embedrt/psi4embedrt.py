import os
import sys
import time
import shutil
import os.path

#sys.path.append("/home/matteod/build/xcfun/build/lib/python")
#sys.path.append("/home/matteod/psi4conda/lib/python3.7")
#sys.path.append("/home/matteod/pybertha/psi4rt")
#sys.path.append("/home/matteod/pybertha/src")
#sys.path.append("/home/matteod/build/pyadf/src")

#os.environ['PSIPATH']="/home/redo/BERTHAEmb/psi4conda/share/psi4/basis"
#os.environ['PYBERTHAROOT'] = "/home/matteod/pybertha/"
#os.environ['RTHOME'] = "/home/matteod/pybertha/psi4rt"
#sys.path.append(os.environ['PYBERTHAROOT']+"/src")
#sys.path.append(os.environ['RTHOME'])
#sys.path.append(os.environ['PSIPATH'])

import psi4
import numpy as np
import argparse

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


def finalize_stdout_redirect (fname, writef=False):


    global stdout_fileno 
    global stdout_save 
    global stdout_pipe 
    global captured_stdout
      
    os.close(stdout_pipe[0])
    os.dup2(stdout_save, stdout_fileno)
    os.close(stdout_save)

    fp = None
    if writef:
      fp = open(psioufname, "w")
    else:
      fp = open(psioufname, "a")

    for line in captured_stdout:
      fp.write(line)
    fp.close()

########################################################################################################

def scfiterations (args, maxiter, jk, H, Cocc, func, wfn, D, vemb, E, Eold, \
  Fock_list, DIIS_error, phi, agrid, densgrad, denshess, \
    isolated_elpot_enviro, E_conv, D_conv, restricted = True):

    #outer loop
    for OUT_ITER in range(0,maxiter):
        
        #inner loop
        for SCF_ITER in range(1, maxiter + 1):
    
            # Compute JK
            jk.C_left_add(Cocc)
            jk.compute()
            jk.C_clear()
       
            # Build Fock matrix
            F = H.clone()
            F.axpy(2.0, jk.J()[0])
            #XC potential
            #sup = psi4.driver.dft_funcs.build_superfunctional(func, restricted)[0]
            sup = psi4.driver.dft.build_superfunctional(func, restricted)[0]
            sup.set_deriv(2)
            sup.allocate()
            vname = "RV"
            if not restricted:
                vname = "UV"
            potential=psi4.core.VBase.build(wfn.basisset(),sup,vname)
            potential.initialize()
            potential.set_D([D])
            V=psi4.core.Matrix(nbf,nbf)
            potential.compute_V([V])
            potential.finalize()
            Exc= potential.quadrature_values()["FUNCTIONAL"]
            # hybrid update
            if sup.is_x_hybrid():
              alpha = sup.x_alpha()
              K = jk.K()[0]
              F.axpy(-alpha,K)
              Exc += -alpha*np.trace(np.matmul(D,K))
            ####### 
            F.axpy(1.0, V)
            F.axpy(1.0, vemb)
            if args.static_field:
                fmax = args.fmax 
                fdir = args.fdir
                F.axpy(fmax,dipole[fdir])
                print("fdir = %4.16f "%(fdir))
                print("fmax = %4.16f "%(fmax))

            twoel = 2.00*np.trace(np.matmul(np.array(D),np.array(jk.J()[0])))
            # DIIS error build and update
            diis_e = psi4.core.Matrix.triplet(F, D, S, False, False, False)
            diis_e.subtract(psi4.core.Matrix.triplet(S, D, F, False, False, False))
            diis_e = psi4.core.Matrix.triplet(A, diis_e, A, False, False, False)
       
            diis.add(F, diis_e)
       
            # SCF energy and update
            SCF_E = 2.0*H.vector_dot(D) + Enuc + Exc + twoel
            if args.static_field:
                     SCF_E = SCF_E + 2.0*np.trace(np.matmul(np.array(dipole[fdir]),np.array(D)))*fmax 
       
            dRMS = diis_e.rms()
       
            ztmp= np.matmul(np.array(D),np.array(dipole[2]))
            diptmp = np.trace(ztmp)
            print("SCF Iteration %3d: Energy = %4.16f "%(SCF_ITER, SCF_E))
            print("   dE = %1.5E"%(SCF_E - Eold))
            print(" dRMS = %1.5E"%(dRMS))
            print("   Pz = %1.5E"%(diptmp))
            
            if (abs(SCF_E - Eold) < E_conv) and (dRMS < D_conv):
                diffD= D.np-Dold
                norm_D=np.linalg.norm(diffD,'fro')
                print("norm_D at OUT_ITER(%i): %.12f\n" % (OUT_ITER,norm_D))
                break
       
            Eold = SCF_E
            Dold = np.copy(D)
       
            F = psi4.core.Matrix.from_array(diis.extrapolate())
       
            # Diagonalize Fock matrix
            orbene, C, Cocc, D = build_orbitals(F)
       
            if SCF_ITER == maxiter:
                psi4.core.clean()
                raise Exception("Maximum number of SCF cycles exceeded.\n")
            #end inner loop
        
        if ( ((not args.nosscf) and (OUT_ITER > 1)) and (not args.nofde) ) :
             
             diffv= vemb_in - vemb.np
             diffD= D.np-D_in
             norm_D=np.linalg.norm(diffD,'fro')
             norm_v=np.linalg.norm(diffv,'fro')
             if (norm_D<(1.0e-6) and norm_v<(1.0e-8)):
                 break
             else:
                print("norm_D : %.12f\n" % norm_D)
                print("norm_v : %.12f\n" % norm_v)

        if ( (not args.nosscf) and (not args.nofde) ):
            # calc new emb potential
     
            # Here the calculation of vemb(D_inner) starts
            # if needed elpot of the active system can be built on the fly
      
            #phi matrix already built
      
            #export D to grid repres.
      
            #fde_util.dens2grid(phi,D,xs,ys,zs,ws,0)
            #temp = 2.0 * np.einsum('pm,mn,pn->p', phi, np.copy(D), phi)
            temp = 2.0 * fde_util.denstogrid( phi, np.copy(D), S,ndocc)
            rho[:,0] = temp
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
            # call PyADF
            #
      
            #print("Subsystem A:")
      
            dens_gf = GridFunctionFactory.newGridFunction(agrid,np.ascontiguousarray(rho[:,0]), gf_type="density") 
            density_active = GridFunctionContainer([dens_gf, densgrad, denshess]) #refresh container
            #nel_active = density_active[0].integral()
            #print("Integrated number of electrons for subsystem A: %.8f" % nel_active)
            print() 
      
            nadpot_active = embed_eval.get_nad_pot(density_active, isolated_dens_enviro)
      
            print("Adding the electrostatic potential")
      
            embpot_A_2 = isolated_elpot_enviro + nadpot_active
      
            if args.fdecorr : 
              embpot_A_2 = fde_util.fcorr(embpot_A_2,density_active[0],isolated_dens_enviro[0])
            #print("   c. Exporting the potential to file")
      
            pot = embpot_A_2.get_values()
            #GridFunctionWriter.write_xyzwv(embpot_A_2,filename=os.path.join("./", 'EMBPOT_PYEMBED_ADFGRID_H2O'),add_comment=False)
      
            #copy vemb and D of the internal loop
            vemb_in = np.copy(vemb)
            D_in = np.copy(D)
      
            #transform EMBPOT_PYEMB_ADFGRID_H2O in basis set representation
            vemb_new = fde_util.embpot2mat(phi,nbas,pot,ws,lpos)
            vemb = psi4.core.Matrix.from_array(vemb_new) 
            if OUT_ITER == maxiter:
                raise Exception("Maximum number of SCF cycles exceeded.\n")
            print("Outer iteration %i : DONE\n" % OUT_ITER)
        else:
             break # for args.fde = False, quit the outer loop of splitSCF scheme
    #end outer loop

    return D, C, Cocc, F, SCF_E, twoel, Exc, orbene

    # Diagonalize routine

########################################################################################################

def build_orbitals(diag):

    Fp = psi4.core.Matrix.triplet(A, diag, A, True, False, True)
    
    Cp = psi4.core.Matrix(nbf, nbf)
    eigvals = psi4.core.Vector(nbf)
    Fp.diagonalize(Cp, eigvals, psi4.core.DiagonalizeOrder.Ascending)
    
    C = psi4.core.Matrix.doublet(A, Cp, False, False)
    
    Cocc = psi4.core.Matrix(nbf, ndocc)
    Cocc.np[:] = C.np[:, :ndocc]
    
    D = psi4.core.Matrix.doublet(Cocc, Cocc, False, True)
      
    return eigvals, C, Cocc, D

########################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-gA","--geom_act", help="Specify geometry file for active subsystem", required=True, 
            type=str, default="geomA.xyz")
    parser.add_argument("-gB","--geom_env", help="Specify geometry file for environment", required=True, 
            type=str, default="geomB.xyz")
    parser.add_argument("--ghostB", help="Add ghosted B atoms to geomA", 
            required=False, action="store_true", default=False)
    parser.add_argument("-d", "--debug", help="Debug on, prints debug info to err.txt", required=False,
            default=False, action="store_true")
    parser.add_argument("--modpaths", help="set berthamod and all other modules path [\"path1;path2;...\"] (default = ../src)", 
            required=False, type=str, default="../src")
    parser.add_argument("--jumprt", help="jump the RT ropagation", 
            required=False, action="store_true", default=False)

    parser.add_argument("--env_obs", help="Specify the orbital basis set for the enviroment (default: AUG/ADZP)", required=False, 
            type=str, default="AUG/ADZP")

    parser.add_argument("--static_field", help="Add a static field to the SCF (default : False)", required=False,
            default=False, action="store_true")
    parser.add_argument("--fmax", help="Static field amplitude (default : 1.0e-5)", required=False,
            type=np.float64, default=1.0e-5)
    parser.add_argument("--fdir", help="External field direction (cartesian)  (default: 2)",
            required=False, type=int, default=2)

    parser.add_argument("--env_func", help="Specify the function for the environment density (default: BLYP)", required=False, 
            type=str, default="BLYP")
    parser.add_argument("--act_obs", help="Specify the orbital basis set for the active system (deafult: aug-cc-pvdz)", required=False, 
            type=str, default="aug-cc-pvdz")
    parser.add_argument("--act_func", help="Specify the function for the active density (default: blyp)", required=False, 
            type=str, default="blyp")
    parser.add_argument("--grid_opts", help="Set the type of integration grid (1,2,3) 1) adffragmentsjob grid, 2) supramolecular grid (default), 3) active system grid",
            default=2, type = int)

    parser.add_argument("-f", "--nofde", help="FDE off", required=False,
            default=False, action="store_true")
    parser.add_argument("--nosscf", help="No SplitSCF for ground state FDE on", required=False,
            default=False, action="store_true")

    parser.add_argument("-fcorr", "--fdecorr", help="FDE long range correction on", required=False,
            default=False, action="store_true")
    parser.add_argument("-i", "--iterative", help="Vemb is updated during propagation", required=False,
            default=False, action="store_true")
 
    parser.add_argument("-a", "--axis", help="The axis of  electric field direction (x = 0, y = 1, z = 2, default 2)",
            default=2, type = int)
    parser.add_argument("-p", "--period", help="Time period of Vemb update during propagation",
            default=0.1, type = float)
    parser.add_argument("-w","--wkd", help="Specify the working dir", required=False, 
            type=str, default="./")
    parser.add_argument("--select", help="Specify the occ-virt MO weighted dipole moment. (-2; 0 & 0 to activate) (default: 0; 0 & 0)",
            default="0; 0 & 0", type=str)
    parser.add_argument("--acc_int", help="Set integration accuracy (default 4.0)",
            default=4.0, type = float)
    parser.add_argument("--inputfile", help="Set input filename, [default = input.inp]",
            default="input.inp")
    #more option to be added
    
    args = parser.parse_args() #temporary
    rt_nthreads = int(os.getenv('OMP_NUM_THREADS', 1))

    for path in args.modpaths.split(";"):
        sys.path.append(path)

    modpaths = os.environ.get('PYBERTHA_MOD_PATH')

    for path in modpaths.split(";"):
        sys.path.append(path)

    import util
    import rtutil
    import fde_util

    Dir = args.axis
    basis_set = args.act_obs
    geomA = args.geom_act
    geomB = args.geom_env
    WD = args.wkd
    debug = args.debug 
    func = args.act_func
    
    #basis_set : defined from input

    print("Read input ... ")
    
    if not args.jumprt:
      imp_opts, calc_params = util.set_params(args.inputfile)
      func = calc_params['func_type'] # from input.inp. default : blyp
    
    geom, mol = fde_util.set_input(geomA, basis_set,geomB,args.ghostB)
    psi4.set_num_threads(rt_nthreads)

    ene = None 
    wfn_scf = None
    adfoufname = "./adf.out"
    psioufname = "./psi4.out"

    if not args.nofde:
    
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      #ADF part

      import pyadf 
      import pyadf.PyEmbed
    
      from pyadf.Plot.FileWriters import GridWriter, GridFunctionWriter
      from pyadf.Plot.FileReaders import GridFunctionReader
      from pyadf.Plot.GridFunctions import GridFunctionFactory
      from pyadf.Plot.GridFunctions import GridFunctionContainer
    
      import math
    
      adf_settings = pyadf.adfsettings()
      adf_settings.set_save_tapes([21,10])
      adf_settings.set_functional(args.env_func)
      adf_settings.set_convergence(1.0e-8)
      adf_settings.set_integration(accint=args.acc_int)
    
      # other options
      basis_env = args.env_obs
      fde_act_opts = {'TNAD' : 'THOMASFERMI' } #
    
      file_active = os.path.join (WD, geomA)
      file_enviro = os.path.join (WD, geomB)
    
      m_active = pyadf.molecule(file_active)
      m_active.set_symmetry('NOSYM')
      m_enviro = pyadf.molecule(file_enviro)
      m_enviro.set_symmetry('NOSYM')
    
      m_tot = m_active + m_enviro
      m_tot.set_symmetry('NOSYM')
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      # Grid generation step
      #
      # getting a grid from adf. Alternative options:  i) use adffragmentsjob grid | deafult
      #                                               ii) use the total sys (env+act)  grid
      #                                              iii) use the active sys grid
      print("Running ADF single-point...")
      if os.path.isfile(adfoufname):
        print("  Removing "+ adfoufname )
        os.remove(adfoufname)
      for resdir in ["./resultfiles", "./jobtempdir"]:
        if os.path.isdir(resdir):
          print ("  Removing "+ resdir )
          shutil.rmtree(resdir)

      f = io.StringIO()
      with redirect_stdout(f):
        r_isolated_enviro = pyadf.adfsinglepointjob(m_enviro, basis_env, \
          settings=adf_settings, options=['NOSYMFIT']).run()

        if args.grid_opts == 1:
            frags = [ pyadf.fragment(None,  [m_active]),
                    pyadf.fragment(r_isolated_enviro, [m_enviro], isfrozen=True) ]
            fde_res = pyadf.adffragmentsjob(frags, basis_env, settings=adf_settings, options=['NOSYMFIT'])
            fde_res=fde_res.run()
            agrid = pyadf.adfgrid(fde_res)
        elif args.grid_opts == 2:
            r_tot = pyadf.adfsinglepointjob(m_tot, basis_env, settings=adf_settings, options=['NOSYMFIT']).run()
            agrid = pyadf.adfgrid(r_tot)
        elif args.grid_opts == 3:
            r_act = pyadf.adfsinglepointjob(m_active, basis_env, settings=adf_settings, options=['NOSYMFIT']).run()
            agrid = pyadf.adfgrid(r_act)
        elif args.grid_opts == 4:
            adf_settings.set_functional("BLYP")
            frags = [ pyadf.fragment(None,  [m_active]),
                    pyadf.fragment(r_isolated_enviro, [m_enviro], isfrozen=True) ]
            fde_res = pyadf.adffragmentsjob(frags, basis="AUG/ADZP", settings=adf_settings, \
		fde=fde_act_opts, options=['NOSYMFIT\n EXCITATIONS\n  ONLYSING\n  LOWEST 2\nEND'])
            fde_res=fde_res.run()
            agrid = pyadf.adfgrid(fde_res)
        
        GridWriter.write_xyzw(grid=agrid,filename=os.path.join("./", 'ADFGRID'),add_comment=False)

      fp = open(adfoufname, "w")
      fp.write(f.getvalue())
      fp.close()
      print("ADF out  " + adfoufname)
         
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      # Preliminary step for psi4
      #
    
      xs,ys,zs,ws,prop = fde_util.read_gridfunction("ADFGRID")
    
      print("Getting gs density")
      #Vvals, wfn_scf = fde_util.get_elpot(xs,ys,zs)
      if os.path.isfile(psioufname):
        print("  Removing "+ psioufname )
        os.remove(psioufname)
      
      init_stdout_redirect ()
      t = threading.Thread(target=drain_pipe)
      t.start()

      ene, wfn_scf = psi4.energy(func,return_wfn=True)
      psi4.cubeprop(wfn_scf)

      os.close(stdout_fileno)
      t.join()
      finalize_stdout_redirect(psioufname, True)
      print("PSI4 out  " + psioufname)

      os.rename("Da.cube","Da0.cube")
      os.rename("Db.cube","Db0.cube")
      os.rename("Dt.cube","Dt0.cube")
      os.rename("Ds.cube","Ds0.cube")
    else :
      init_stdout_redirect ()
      t = threading.Thread(target=drain_pipe)
      t.start()

      ene, wfn_scf = psi4.energy(func,return_wfn=True)

      os.close(stdout_fileno)
      t.join()
      finalize_stdout_redirect(psioufname, True)
      print("PSI4 out  " + psioufname)

    init_stdout_redirect ()
    t = threading.Thread(target=drain_pipe)
    t.start()

    D = np.array(wfn_scf.Da())
    C0 = np.array(wfn_scf.Ca()) #unpolarized MOs - for later use
    Dref = np.copy(D)
    Eref = ene
    mints = psi4.core.MintsHelper(wfn_scf.basisset())
    S = mints.ao_overlap()
    ndocc = wfn_scf.nalpha() 

    os.close(stdout_fileno)
    t.join()
    finalize_stdout_redirect(psioufname)

    if not args.nofde:
      #build phi, the matrix containing the values of basis set on grid
      print("Building phi matrix")

      init_stdout_redirect ()
      t = threading.Thread(target=drain_pipe)
      t.start()

      phi, lpos, nbas = fde_util.phi_builder(mol,xs,ys,zs,ws,basis_set)

      os.close(stdout_fileno)
      t.join()
      finalize_stdout_redirect(psioufname)
      print("PSI4 out  " + psioufname)

      #export D to grid repres.
      #fde_util.dens2grid(phi,D,xs,ys,zs,ws,0)
      # density on grid
      #temp = 2.0 * np.einsum('pm,mn,pn->p', phi, D, phi)
      temp = 2.0 * fde_util.denstogrid( phi, D, S,ndocc)
      rho = np.zeros((temp.shape[0],10),dtype=np.float_)
      rho[:,0] = temp
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      # PyADF part
      #
      f = io.StringIO()
      with redirect_stdout(f):

        embed_settings = pyadf.EmbedXCFunSettings()
        embed_settings.set_fun_nad_xc ({'lda' : 1.0})
        embed_settings.set_fun_nad_kin({'tfk' : 1.0})
    
        embed_settings.show_functionals()
    
        embed_eval = pyadf.PyEmbed.EmbedXCFunEvaluator(settings=embed_settings)

      fp = open(adfoufname, "a")
      fp.write(f.getvalue())
      fp.close()
      print("  ADF out  " + adfoufname)

      print()
      print("Subsystem A:")
      start_map=time.time()
      dens_gf = GridFunctionFactory.newGridFunction(agrid,np.ascontiguousarray(rho[:,0]), gf_type="density") 
      densgrad = GridFunctionFactory.newGridFunction(agrid, np.ascontiguousarray(rho[:, 1:4]))
      denshess = GridFunctionFactory.newGridFunction(agrid, np.ascontiguousarray(rho[:, 4:10]))  
      density_active = GridFunctionContainer([dens_gf, densgrad, denshess])
      end_map=time.time()
      print("  time to map active density on grid : %.3f" % (end_map-start_map))
      nel_active = density_active[0].integral()
      print("  Integrated number of electrons for subsystem A: %.8f" % nel_active)

      #grid_active, isolated_elpot_active, density_active = \
      #    GridFunctionReader.read_density_elpot_xyzwv(os.path.join("./", 'FRZDNS.a.txt'))
      #nel_active = density_active[0].integral()
      #check that grid_active and agrid should be the same
      #GridWriter.write_xyzw(grid=grid_active,filename=os.path.join("./", 'ADFGRID_H20'),add_comment=False)
      #isolated ammonia get densities and nuclear potentials
  
      f = io.StringIO()
      with redirect_stdout(f):
        isolated_dens_enviro  = r_isolated_enviro.get_density(grid=agrid, \
          fit=False, order=2) #grid=grid_active=agrid
        isolated_vnuc_enviro  = r_isolated_enviro.get_potential(grid=agrid, \
          pot='nuc')
        isolated_coul_enviro  = r_isolated_enviro.get_potential(grid=agrid, \
          pot='coul')
        isolated_elpot_enviro = isolated_vnuc_enviro + isolated_coul_enviro
      
      fp = open(adfoufname, "a")
      fp.write(f.getvalue())
      fp.close()
      print("  ADF out  " + adfoufname)

      print("Evaluation of the embedding potential in two steps:")
      print ("  a. getting the non-additive potential")

      f = io.StringIO()
      with redirect_stdout(f):
        #nadkin_active_dd = embed_eval.get_nad_pot_kin(density_active,isolated_dens_enviro)
        #nadxc_active_dd  = embed_eval.get_nad_pot_xc(density_active, isolated_dens_enviro)
        #embpot_active_dd = nadkin_active_dd+nadxc_active_dd+isolated_elpot_enviro
        nadpot_active = embed_eval.get_nad_pot(density_active, isolated_dens_enviro)
        #nadpot_enviro = embed_eval.get_nad_pot(isolated_dens_enviro, density_active) #not used if the density of environment is not propagated
      
      fp = open(adfoufname, "a")
      fp.write(f.getvalue())
      fp.close()
      print("    ADF out  " + adfoufname)

      print("  b. adding the electrostatic potential")
    
      embpot_A_2 = isolated_elpot_enviro + nadpot_active
      if args.fdecorr : 
        embpot_A_2 = fde_util.fcorr(embpot_A_2, \
          density_active[0],isolated_dens_enviro[0])
      #embpot_B_2 = isolated_elpot_active + nadpot_enviro #not used if the density of environment is not propagated
      #print("empot_A_2 is %s\n" % type(embpot_A_2))
      pot = embpot_A_2.get_values()
      
      #print("   c. Exporting the potential to file")
      #GridFunctionWriter.write_xyzwv(embpot_A_2,filename=os.path.join("./", 'EMBPOT_PYEMBED_ADFGRID_H2O'),add_comment=False)
      #GridFunctionWriter.write_xyzwv(embpot_B_2,filename=os.path.join("./", 'EMBPOT_PYEMBED_ADFGRID_NH3'),add_comment=False)
      print()
      #print "2. Evaluation of the embedding potential in one step, exporting to file"
      #embpot_A_1 = embed_eval.get_emb_pot(density_A, density_B, ve_B)
      #embpot_B_1 = embed_eval.get_emb_pot(density_B, density_A, ve_A)
      #transform EMBPOT_PYEMB_ADFGRID_H2O in basis set representation
      res = fde_util.embpot2mat(phi,nbas,pot,ws,lpos)
    else:
      res = np.zeros_like(D)
      phi = None
      agrid = None
      densgrad = None 
      denshess = None
      isolated_elpot_enviro = None


    #########################################################
    """
    Restricted Kohn-Sham code using the Psi4 JK class for the
    4-index electron repulsive integrals and Vbase class for
    exchange-correlation potential. Based on:
    """
    
    """
    A restricted Hartree-Fock code using the Psi4 JK class for the 
    4-index electron repulsion integrals.
    
    References:
    - Algorithms from [Szabo:1996], [Sherrill:1998], and [Pulay:1980:393]
    """
    
    __authors__ = "Daniel G. A. Smith"
    __credits__ = ["Daniel G. A. Smith"]
    
    __copyright__ = "(c) 2014-2018, The Psi4NumPy Developers"
    __license__ = "BSD-3-Clause"
    __date__ = "2017-9-30"
    
    np.set_printoptions(precision=5, linewidth=200, suppress=True)
    import helper_HF
    
    ftime = open("time_stats.txt", "w")
    # Set tolerances
    maxiter = 200
    E_conv = 1.0E-8
    if not args.nosscf:
      D_conv = 1.0E-8
    else :
      D_conv = 1.0E-12

    init_stdout_redirect ()
    t = threading.Thread(target=drain_pipe)
    t.start()

    # Integral generation from Psi4's MintsHelper
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('BASIS')) #use basis_set instead
    #mints = psi4.core.MintsHelper(wfn.basisset()) # moved upward
    #S = mints.ao_overlap()                        #
    
    # Get nbf and ndocc for closed shell molecules
    nbf = wfn.nso()
    #ndocc = wfn_scf.nalpha() # moved upward
    if wfn.nalpha() != wfn.nbeta():
        raise PsiException("Only valid for RHF wavefunctions!")
    
    os.close(stdout_fileno)
    t.join()
    finalize_stdout_redirect(psioufname)
    print("PSI4 out  " + psioufname)
    
    print('Number of occupied orbitals: %d' % ndocc)
    print('Number of basis functions:   %d' % nbf)
    
    # read from vemb.txt
    #temp = np.zeros((nbf,nbf),dtype=np.float_)
    #j = 0
    #with open("vemb.txt","r") as f:
    # for line in f:
    #  raw = line.split()
    #  temp[j,:]=raw
    #  j += 1
    
    # check vemb matrix
    vemb = psi4.core.Matrix.from_array(res)
    
    # Build H_core
    V = mints.ao_potential()
    T = mints.ao_kinetic()
    H = T.clone()
    H.add(V)
    dipole=mints.ao_dipole()
    #dipole is  a list of psi4.core.Matrix objects
    # dipole[0] -> x, dipole[1] -> y, dipole[2] -> z
    Ndip= mol.nuclear_dipole()
    
    # Orthogonalizer A = S^(-1/2)
    A = mints.ao_overlap()
    A.power(-0.5, 1.e-16)

    # Build diis
    diis = helper_HF.DIIS_helper(max_vec=6)
    
    # Build core orbitals
    #C, Cocc, D = build_orbitals(H)
    #as initial guess we may use the ground state isolated density
    Cocc = wfn_scf.Ca_subset("AO","OCC")
    D = wfn_scf.Da()
    
    # Setup data for DIIS
    t = time.time()
    E = 0.0
    Enuc = mol.nuclear_repulsion_energy()
    Eold = 0.0
    #Dold = psi4.core.Matrix(nbf, nbf)
    Fock_list = []
    DIIS_error = []
    
    # Initialize the JK object
    jk = psi4.core.JK.build(wfn.basisset())
    jk.set_memory(int(1.25e8))  # 1GB
    jk.initialize()
    #jk.print_header()
    
    ftime.write('\nTotal time taken for setup: %.3f seconds\n' % (time.time() - t))
    
    print('\nStart SCF iterations:\n\n')
    start = time.time()
    cstart = time.process_time()

    D, C, Cocc, F, SCF_E, twoel, Exc, orb_eigv = scfiterations (args, maxiter, jk, H, Cocc, func, \
      wfn, D, vemb, E, Eold, Fock_list, DIIS_error,  phi, agrid, densgrad, denshess, \
        isolated_elpot_enviro, E_conv, D_conv)
    
    end = time.time()
    cend = time.process_time()
    ftime.write('Total time[time()] for SCF iterations: %.3f seconds \n\n' % (-start + end))
    ftime.write('Total time[clock()] for SCF iterations: %.3f seconds \n\n' % (-cstart + cend))
    
    print('Final DFT energy: %.16f hartree' % SCF_E)
    print('    twoel energy: %.16f hartree' % twoel)
    print('       xc energy: %.16f hartree' % Exc)
    
    print('Orbital Energies [Eh]\n')
    print('Doubly Occupied:\n')
    for k in range(ndocc):
         print('%iA : %.6f' %(k+1,orb_eigv.np[k]))
    print('Virtual:\n')

    for k in range(ndocc,nbf):
         print('%iA : %.6f'% (k+1,orb_eigv.np[k]))

    dipz = np.matmul(np.array(D),np.array(dipole[2]))
    dipvalz = np.trace(2.0*dipz) #+ Ndip[2]
    dipy = np.matmul(np.array(D),np.array(dipole[1]))
    dipvaly = np.trace(2.0*dipy) #+ Ndip[1]
    dipx = np.matmul(np.array(D),np.array(dipole[0]))
    dipvalx = np.trace(2.0*dipx) #+ Ndip[0]
    fout = open("emb_res.txt", "w")
    conv = 2.541765 # 1 unit of electric dipole moment (au) = 2.541765 Debye
    conv = 1.0 # for atomic units
    fout.write('polarized (embedding) dipole:  \n')
    if args.static_field:
     fout.write('External Field direction: %.16f \n' % args.fdir)
     fout.write('External Filed strength (a.u.): %.7f \n' %args.fmax)
    fout.write('polarized electron dipole: %.16f, %.16f, %.16f a.u\n' %  (dipvalx*conv,dipvaly*conv,dipvalz*conv))
    fout.write('Nuclear dipole:            %.16f, %.16f, %.16f a.u\n' %  (Ndip[0]*conv,Ndip[1]*conv,Ndip[2]*conv))
    fout.write('total dipole:              %.16f, %.16f, %.16f a.u\n' %  ((dipvalx+Ndip[0])*conv,(dipvaly+Ndip[1])*conv,(dipvalz+Ndip[2])*conv))
    fout.write('Final DFT energy: %.16f hartree\n' % (SCF_E)) 

    dipz = np.matmul(np.array(Dref),np.array(dipole[2]))
    dipvalz = np.trace(2.0*dipz) #+ Ndip[2]
    dipy = np.matmul(np.array(Dref),np.array(dipole[1]))
    dipvaly = np.trace(2.0*dipy) #+ Ndip[1]
    dipx = np.matmul(np.array(Dref),np.array(dipole[0]))
    dipvalx = np.trace(2.0*dipx) #+ Ndip[0]

    fout.write('unpolarized (no embedding and no external field) dipole: \n')
    fout.write('unpolarized electron dipole: %.16f, %.16f, %.16f a.u\n' %  (dipvalx*conv,dipvaly*conv,dipvalz*conv))
    fout.write('Nuclear dipole:              %.16f, %.16f, %.16f a.u\n' %  (Ndip[0]*conv,Ndip[1]*conv,Ndip[2]*conv))
    fout.write('total unpolarized  dipole:   %.16f, %.16f, %.16f a.u\n' %  ((dipvalx+Ndip[0])*conv,(dipvaly+Ndip[1])*conv,(dipvalz+Ndip[2])*conv))
    print("Ndip : %.8f, %.8f, %.8f\n" % (Ndip[0],Ndip[1],Ndip[2]))
    fout.write('Final DFT energy: %.16f hartree' % (Eref)) 
    fout.close()

    init_stdout_redirect ()
    t = threading.Thread(target=drain_pipe)
    t.start()

    wfn_scf.Da().copy(D)
    wfn_scf.Db().copy(D)
    psi4.cubeprop(wfn_scf)

    os.close(stdout_fileno)
    t.join()
    finalize_stdout_redirect(psioufname)
    print("PSI4 out  " + psioufname)

    if (args.jumprt):
      exit()

    # check of D and Cocc
    ovap = np.array(S)
    import scipy.linalg
    R = np.matmul(ovap,np.matmul(D,ovap))
    eigvals, eigvecs = scipy.linalg.eigh(R,ovap, eigvals_only=False)
    np.savetxt("eigvals.txt", eigvals)
    Cocc_test = eigvecs[:,-ndocc:]
    Dtest = np.matmul(Cocc_test,np.conjugate(Cocc_test.T))
    print("D and Dtest : %s" % (np.allclose(D,Dtest,atol=1.0e-12)))
    #for free-field propagation Fmax = 0.0
    
    molist = args.select.split("&")
    occlist = molist[0].split(";")
    occlist = [int(m) for m in occlist]
    do_weighted = occlist.pop(0)
    virtlist = molist[1].split(";")
    virtlist = [int(m) for m in virtlist]
    
    if (do_weighted == -2):
      if debug:
        print("Selected transitions from %s to %s MOs"% (str(occlist), str(virtlist)))
    
    numpy_memory = 8
    #initialize mints object
    mints = psi4.core.MintsHelper(wfn.basisset())
    S=np.array(mints.ao_overlap())
    #get T,V,dipole_z
    T=np.array(mints.ao_kinetic())
    V=np.array(mints.ao_potential())
    dipole=mints.ao_dipole()
    #dipole is  a list of psi4.core.Matrix objects
    # dipole[0] -> x, dipole[1] -> y, dipole[2] -> z
    direction = args.axis
    dipz_mat=np.copy(dipole[direction])
    H=T+V
    #internal defined functional 
    #svwn5_func = {
    #    "name": "SVWN5",
    #    "x_functionals": {
    #        "LDA_X": {}
    #    },
    #    "c_functionals": {
    #        "LDA_C_VWN": {}
    #    }
    #}
    #nalpha, ndocc, nbf have already been defined
    #from psi4-rt.py
    print('Number of occupied orbitals: %d' % ndocc)
    print('Number of basis functions: %d' % nbf)
    
    # Run a quick check to make sure everything will fit into memory
    I_Size = (nbf**4) * 8.e-9
    print("Size of the ERI tensor will be %4.2f GB." % I_Size)
    
    # Estimate memory usage
    memory_footprint = I_Size * 1.5
    if I_Size > numpy_memory:
        psi4.core.clean()
        raise Exception("Estimated memory utilization (%4.2f GB) exceeds numpy_memory limit of %4.2f GB." % (memory_footprint, numpy_memory))
    #Get Eri (2-electron repulsion integrals)
    I=np.array(mints.ao_eri())
    
    #D_0 is the density in the reference molecular basis
    D_0=np.zeros((nbf,nbf))
    for num in range(int(ndocc)):
        D_0[num,num]=1.0
    Dp_0 = D_0
    
    #containers
    ene_list = []
    dip_list = []
    imp_list=[]
    weighted_dip = []
    
    #unpolarized MOs
    #C0 and C0_inv 
    C0_inv=np.linalg.inv(C0)
    dipz_mo0 = np.matmul(np.conjugate(C0.T),np.matmul(dipz_mat,C0))
    
    #C_inv used to backtransform D(AO)
    C = np.array(C)
    C_inv=np.linalg.inv(C)
    #check inversion
    test = np.matmul(C,C_inv)
    Id = np.eye(nbf)
    print("CC^-1 = 1 : %s\n" % np.allclose(Id,test,atol=1.0e-12))
    #check C orthornmality
    test = np.matmul(np.conjugate(C.T),np.matmul(S,C))
    print("C(.H)SC = 1 : %s\n" % np.allclose(Id,test,atol=1.0e-12))
    #set propagation params
    
    
    if imp_opts['imp_type'] == 'analytic' :
        analytic = True
    else:
        analytic = False
    
    dt =  calc_params['delta_t']
    #time_int in atomic unit
    time_int=calc_params['time_int']
    niter=int(time_int/dt)
    
    basisset=wfn.basisset()
    
    #rt summary 
    fo  = open("err.txt", "w")
    fo.write("total time :  %5.2f \ndt : %.3f \nFmax : %.8f\nnoFDE : %s\nno Split_scf : %s\nVemb update : %s\n\
      period : %5.2f\n" % (time_int,dt,imp_opts['Fmax'],\
      args.nofde, args.nosscf, args.iterative,args.period))
    #
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    # analytic pertubation goes here
    #
    
    if (analytic):
       print('Perturb density with analytic delta')
       # set the perturbed density -> exp(-ikP)D_0exp(+ikP)
       k = imp_opts['Fmax']
       
       #dipz_mat is transformed to the reference MO basis
       dip_mo=np.matmul(np.conjugate(C.T),np.matmul(dipz_mat,C))
       u0=util.exp_opmat(dip_mo,np.float_(-k))
       Dp_init= np.matmul(u0,np.matmul(Dp_0,np.conjugate(u0.T)))
       func_t0=k
       #backtrasform Dp_init
       D_init=np.matmul(C,np.matmul(Dp_init,np.conjugate(C.T)))
       D = D_init
       Dp_0 = Dp_init
     
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    # get a new vemb matrix
    #
    #get new density and elpot
    if not args.nofde:
      print("Start calculation of vemb(D_pol) starts")
      start=time.time()
      cstart =time.process_time()
      #if needed elpot of the active system can be built on the fly
      #end_v=time.time()
      #cend_v=time.process_time()
      #cdiff_v = end_v -start
      #diff_v = cend_v -cstart
      #phi matrix already built
    
      #export D to grid repres.
      start_d2g=time.time()
      cstart_d2g=time.process_time()
      #temp = 2.0 * np.einsum('pm,mn,pn->p', phi, D, phi)
      #rho[:,0] = temp
      end_d2g=time.time()
      cend_d2g=time.process_time()
      diff_d2g = end_d2g- start_d2g 
      cdiff_d2g = cend_d2g- cstart_d2g 
      #
      #def denstogrid(phi,D,S):
      start_dg=time.time()
      cstart_dg=time.process_time()
      temp = 2.0 * fde_util.denstogrid( phi, np.array(D), S,ndocc)
      rho[:,0] = temp
      end_dg=time.time()
      cend_dg=time.process_time()
      diff_dg = end_dg- start_dg 
      cdiff_dg = cend_dg- cstart_dg 
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
      # call PyADF
      #
    
      #print("Subsystem A:")
      start_pyin=time.time()
      cstart_pyin=time.process_time()
      dens_gf = GridFunctionFactory.newGridFunction(agrid,np.ascontiguousarray(rho[:,0]), gf_type="density") 
      density_active = GridFunctionContainer([dens_gf, densgrad, denshess]) #refresh container
      #nel_active = density_active[0].integral()
      #print("Integrated number of electrons for subsystem A: %.8f" % nel_active)
      #print() 
      end_pyin =time.time()
      cend_pyin =time.process_time()
      diff_pyin = end_pyin -start_pyin
      cdiff_pyin = cend_pyin -cstart_pyin
    
      start_xfun=time.time()
      cstart_xfun=time.process_time()
      nadpot_active = embed_eval.get_nad_pot(density_active, isolated_dens_enviro)
    
      print("   b. adding the electrostatic potential")
    
      embpot_A_2 = isolated_elpot_enviro + nadpot_active
      end_xfun=time.time()
      cend_xfun=time.process_time()
      diff_xfun = end_xfun-start_xfun
      cdiff_xfun = cend_xfun-cstart_xfun
    
      if args.fdecorr : 
       embpot_A_2 = fde_util.fcorr(embpot_A_2,density_active[0],isolated_dens_enviro[0])
      #print("   c. Exporting the potential to file")
      pot = embpot_A_2.get_values()
      #start_v2f=time.time()
      #cstart_v2f=time.process_time()
      #GridFunctionWriter.write_cube(embpot_A_2,filename=os.path.join("./", 'EMBPOT_PYEMBED_ADFGRID_H2O.cube'),add_comment=False)
      #end_v2f=time.time()
      #cend_v2f=time.process_time()
      #diff_v2f = end_v2f-start_v2f
      #cdiff_v2f = cend_v2f-cstart_v2f
      #transform EMBPOT_PYEMB_ADFGRID_H2O in basis set representation
      start_2mat=time.time()
      cstart_2mat=time.process_time()
      vemb = fde_util.embpot2mat(phi,nbas,pot,ws,lpos)
      
      #vemb = np.copy(res) #in case we want to use the previous vemb
      end = time.time()
      cend = time.process_time()
      ftime.write('Time[time()] for pre-propagation vemb calculation : %.3f seconds \n\n' % (end-start))
      ftime.write('Time[clock()] for pre-propagation vemb calculation : %.3f seconds \n\n' % (cend-cstart))
      #ftime.write('time and ctime for calculation of electrostatic potential (grid_esp) : %.3f , %.3f seconds \n\n' % (diff_v,cdiff_v))
      ftime.write('time and ctime to map D on grid: %.3f , %.3f seconds \n\n' % (diff_d2g,cdiff_d2g))
      ftime.write('time and ctime to map D on grid through MOs: %.3f , %.3f seconds \n\n' % (diff_dg,cdiff_dg))
      ftime.write('time and ctime to import D in PyADF : %.3f , %.3f seconds \n\n' % (diff_pyin,cdiff_pyin))
      ftime.write('time and ctime for Xfun : %.3f , %.3f seconds \n\n' % (diff_xfun,cdiff_xfun))
      #ftime.write('time and ctime for Vemb export to grid (file) : %.3f , %.3f seconds \n\n' % (diff_v2f,cdiff_v2f))
      ftime.write('time and ctime for Vemb matricial repr. : %.3f , %.3f seconds \n\n' % (end-start_2mat,cend-cstart_2mat))
      
    else:
      vemb = np.zeros_like(D)
    
    print('Entering in the first step of propagation')
    
    J0,Exc0,func_t0,F_t0,fock_mid_init=util.mo_fock_mid_forwd_eval(np.array(D),F,\
                         0,np.float_(dt),H,I,dipz_mat,C,C_inv,S,nbf,imp_opts,func,fo,basisset,vemb)
    
    #check hermicity of fock_mid_init
    
    Ah=np.conjugate(fock_mid_init.T)
    fo.write('Fock_mid hermitian: %s\n' % np.allclose(fock_mid_init,Ah))
    
    #propagate D_t0 -->D(t0+dt)
    #
    #fock_mid_init is transformed in the MO ref basis
    fockp_mid_init=np.matmul(np.conjugate(C.T),np.matmul(fock_mid_init,C))
    
    #u=scipy.linalg.expm(-1.j*fockp_mid_init*dt)
    u=util.exp_opmat(fockp_mid_init,np.float_(dt))
    
    temp=np.matmul(Dp_0,np.conjugate(u.T))
    
    Dp_t1=np.matmul(u,temp)
    
    #check u if unitary
    test_u=np.matmul(u,np.conjugate(u.T))
    fo.write('U is unitary :%s\n' % np.allclose(test_u,np.eye(u.shape[0])))
    
    fock_mid_backwd=np.copy(fock_mid_init)
    
    #backtrasform Dp_t1
    
    D_t1=np.matmul(C,np.matmul(Dp_t1,np.conjugate(C.T)))
    
    if (func == 'hf'):
        ene_list.append(np.trace(np.matmul(np.array(D),(H+F_t0))))
    else:    
        ene_list.append(2.00*np.trace(np.matmul(np.array(D),H))+J0+Exc0)
    dip_list.append(np.trace(np.matmul(np.array(D),dipz_mat)))
    
    #weighted dipole
    if (do_weighted == -2):
      D_mo0 = np.matmul(C0_inv,np.matmul(D,np.conjugate(C0_inv.T)))
      res = util.dipoleanalysis(dipz_mo0,D_mo0,ndocc,occlist,virtlist,debug,False)
      weighted_dip.append(res)
    
    fock_mid_backwd=np.copy(fock_mid_init) #prepare the fock at the previous midpint
    D_ti=D_t1
    Dp_ti=Dp_t1
    #Enuc_list.append(-func_t0*Ndip_dir+Nuc_rep) #just in case of non-zero nuclear dipole
    #
    imp_list.append(func_t0)
    if debug :  
      #trace of D_t1
      fo.write('%.8f\n' % np.trace(Dp_ti).real)
      fo.write('Trace of DS %.8f\n' % np.trace(np.matmul(S,D_ti)).real)
      fo.write('Trace of SD.real %.14f\n' % np.trace(np.matmul(S,D_ti.real)))
      fo.write('Trace of SD.imag %.14f\n' % np.trace(np.matmul(S,D_ti.imag)))
      fo.write('Dipole %.8f %.15f\n' % (0.000, 2.00*dip_list[0].real))
    ct_acc= 0.0
    t_acc= 0.0
    start_tot = time.time()
    cstart_tot =time.process_time()
    vcount = 0.0

    print("Startinr propagation...")

    for j in range(1,niter+1):
        # now we have to update vemb :(
        if ((not args.nofde) and args.iterative):
          if ( ( j % int(args.period/dt) ) == 0.0 ):
            start_fde = time.time()
            cstart_fde = time.process_time()
            #get elpot
            #if needed elpot of active system can be built on the fly
            #phi matrix already built
            #export D to grid repres.
            #temp = 2.0 * np.einsum('pm,mn,pn->p', phi, D_ti, phi)
            temp = 2.0 * fde_util.denstogrid( phi, D_ti, S,ndocc)
            rho[:,0] = temp
     
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
            # call PyADF
            #print("Subsystem A:")
     
            dens_gf = GridFunctionFactory.newGridFunction(agrid,np.ascontiguousarray(rho[:,0]), gf_type="density") 
            density_active = GridFunctionContainer([dens_gf, densgrad, denshess]) #refresh container
            #nel_active = density_active[0].integral()
            #print("Integrated number of electrons for subsystem A: %.8f" % nel_active)
            #print() 
     
            nadpot_active = embed_eval.get_nad_pot(density_active, isolated_dens_enviro)
     
            #print("   b. adding the electrostatic potential")
     
            embpot_A_2 = isolated_elpot_enviro + nadpot_active
            if args.fdecorr : 
              embpot_A_2 = fde_util.fcorr(embpot_A_2,density_active[0],isolated_dens_enviro[0])
            pot = embpot_A_2.get_values()
            vcount+=1
            #print("   c. Exporting the potential to file")
            #GridFunctionWriter.write_xyzwv(embpot_A_2,filename=os.path.join("./", 'EMBPOT_PYEMBED_ADFGRID_H2O'),add_comment=False)
            #transform EMBPOT_PYEMB_ADFGRID_H2O in basis set representation
            vemb = fde_util.embpot2mat(phi,nbas,pot,ws,lpos)
            end_fde = time.time()
            cend_fde = time.process_time()
            diff = -start_fde  + end_fde
            cdiff = -cstart_fde + cend_fde
            t_acc += diff
            ct_acc += cdiff
        #print("here, type(vemb) : %s\n" % type(vemb))
        #
        J_i,Exc_i,func_ti,F_ti,fock_mid_tmp=util.mo_fock_mid_forwd_eval(np.copy(D_ti),fock_mid_backwd,\
                               j,np.float_(dt),H,I,dipz_mat,C,C_inv,S,nbf,imp_opts,func,fo,basisset,vemb)
        fo.write('%.8f\n' % np.trace(np.matmul(S,D_ti)).real)
        Ah=np.conjugate(fock_mid_tmp.T)
        fo.write('Fock_mid hermitian: %s\n' % np.allclose(fock_mid_tmp,Ah))
        #transform fock_mid_init in MO basis
        fockp_mid_tmp=np.matmul(np.conjugate(C.T),np.matmul(fock_mid_tmp,C))
        u=util.exp_opmat(np.copy(fockp_mid_tmp),np.float_(dt))
        #u=scipy.linalg.expm(-1.0j*fockp_mid_tmp*dt)
        #check u is unitary
        test_u=np.matmul(u,np.conjugate(u.T))
        if (not np.allclose(np.eye(u.shape[0]),test_u)):
            print('U is not unitary\n')
        
        #check the trace of density to evolve
        fo.write('tr of density to evolve: %.8f\n' % np.trace(Dp_ti).real)
        
        #evolve the density in orthonormal basis
        temp=np.matmul(Dp_ti,np.conjugate(u.T))
        Dp_ti_dt=np.matmul(u,temp)
    
        #backtransform Dp_ti_dt
        D_ti_dt=np.matmul(C,np.matmul(Dp_ti_dt,np.conjugate(C.T)))
        fo.write('%.8f\n' % np.trace(Dp_ti_dt).real)
        #dipole expectation for D_ti
        dip_list.append(np.trace(np.matmul(dipz_mat,D_ti)))
        if (do_weighted == -2):
          #weighted dipole 
          D_ti_mo0 = np.matmul(C0_inv,np.matmul(D_ti,np.conjugate(C0_inv.T)))
          res = util.dipoleanalysis(dipz_mo0,D_ti_mo0,ndocc,occlist,virtlist,debug,False)
          weighted_dip.append(res)
        #Energy expectation value at t = t_i 
        #Enuc_list.append(-func_ti*Ndip_dir+Nuc_rep) #just in case of non-zero nuclear dipole
        if (func=='hf'):
            ene_list.append(np.trace(np.matmul(D_ti,(H+F_ti))))
        else:
            ene_list.append(2.00*np.trace(np.matmul(D_ti,H))+J_i+Exc_i)
        imp_list.append(func_ti)
       
        #update D_ti and Dp_ti for the next step
        
        if debug :
          fo.write('here I update the matrices Dp_ti and D_ti\n')
        D_ti=np.copy(D_ti_dt)
        Dp_ti=np.copy(Dp_ti_dt)
        #for debug
        if debug:
          fo.write('%.8f\n' % np.trace(Dp_ti).real)
          fo.write('%.8f\n' % np.trace(np.matmul(S,D_ti)).real)
          fo.write('Dipole  %.8f %.15f\n' % (j*dt, 2.00*dip_list[j].real))
          fo.write('Vemb up : %3.3f %s\n' % (vcount/(time_int/args.period)*100.,'%'))
    
        fo.flush()
    
        #update fock_mid_backwd for the next step
        fock_mid_backwd=np.copy(fock_mid_tmp)
    
    end_tot = time.time()
    cend_tot =time.process_time()
    ftime.write('Cumulative time[time()] due to vemb calculation : %.3f seconds \n\n' % (t_acc))
    ftime.write('Cumulative time[clock()] due to vemb calculation : %.3f seconds \n\n' % (ct_acc))
    ftime.write('Cumulative time[time()] due to %i iter : %.8e seconds \n\n' % ((niter),(end_tot-start_tot)))
    ftime.write('Cumulative time[clock()] due to %i iter : %.8e seconds \n\n' % ((niter),(cend_tot-cstart_tot)))
    fo.close()
    t_point=np.linspace(0.0,niter*dt,niter+1)
    
    if (do_weighted == -2):
      wd_dip=2.00*np.array(weighted_dip).real
      np.savetxt('weighteddip.txt', np.c_[t_point,wd_dip], fmt='%.12e')
    print('Total number of Vemb eval %i\n' % vcount)
    np.savetxt('dipole.txt', np.c_[t_point,np.array(dip_list).real*2.0], fmt='%.14e')
    np.savetxt('ene.txt', np.c_[t_point,np.array(ene_list).real], fmt='%.14e')
