import argparse
import shutil
import numpy
import sys
import os

import os.path
import time

from dataclasses import dataclass
from pathlib import Path

MAXIT = 100

@dataclass
class pyberthaembedoption:
    activefile: str
    envirofile :str 
    fitcoefffile: str
    vctfile: str
    ovapfile: str
    dumpfiles: bool
    debug: bool
    linemb: bool
    verbosity: int
    thresh: numpy.float64
    wrapperso: str
    berthamodpath: str
    eda_nocv_info: bool
    eda_nocv_frag_file: str
    thresh_conv: numpy.float64 = 1.0e-8
    inputfile: str = "input.inp"
    fittfile: str = "fitt2.inp"
    gtype: int = 2
    param: float = 4.0
    basis: str = 'AUG/ADZP'
    excfuncenv : str = "BLYP"
    static_field : bool = False
    fmax : numpy.float64 = 1.0e-5
    fdir: int = 2
    densityzero: str = ""
    density : str = ""
    drx: float = 0.1
    dry: float = 0.1
    drz: float = 0.1
    margin: float = 5.5
    restartfname : str = "restart.npy.npz" 
    gridfname : str = "grid.dat"

##########################################################################################

def runspberthaembed (pberthaopt, restart = False, stdoutprint = True):

    ovapm = None 
    eigem = None 
    fockm = None  
    eigen = None

    verbosity = -1
    verbosity = pberthaopt.verbosity

    if not stdoutprint:
        verbosity = 0

    import berthamod
    import pyembmod
    
    if stdoutprint:
    	print("Options: ")
    	for att in [a for a in dir(pberthaopt) if not a.startswith('__')]:
        	print(att, " = ", getattr(pberthaopt, att)) 
    	print("")
    	print("")

    	print("")
    
    if not os.path.isfile(pberthaopt.wrapperso):
        raise Exception("SO File ", pberthaopt.wrapperso, " does not exist")

    fittcoefffname = pberthaopt.fitcoefffile
    vctfilename = pberthaopt.vctfile
    ovapfilename = pberthaopt.ovapfile
    
    fnameinput = pberthaopt.inputfile
    if not os.path.isfile(fnameinput):
        raise Exception("File ", fnameinput, " does not exist")
    
    fittfname = pberthaopt.fittfile
    if not os.path.isfile(fittfname):
        raise Exception("File ", fittfname , " does not exist")
    
    dumpfiles = int(pberthaopt.dumpfiles)

    static_field = pberthaopt.static_field
    fmax = pberthaopt.fmax 
    fdir = pberthaopt.fdir
    nofde= pberthaopt.nofde
    eigem = None
    rho = None
    Da0 = None
    etotal = 0.0
    dipx_ref = 0.0
    dipy_ref = 0.0
    dipz_ref = 0.0

    bertha = berthamod.pybertha(pberthaopt.wrapperso)

    if not restart:
        try:
            os.remove(pberthaopt.restartfname)
        except OSError:
            pass

        bertha.set_fittcoefffname(fittcoefffname)
        bertha.set_ovapfilename(ovapfilename)
        bertha.set_vctfilename(vctfilename)
        bertha.set_fnameinput(fnameinput)
        bertha.set_fittfname(fittfname)
        bertha.set_thresh(pberthaopt.thresh)
    
        bertha.set_verbosity(verbosity)
        bertha.set_dumpfiles(dumpfiles)
    
        bertha.set_densitydiff(1)
    
        bertha.init()
    
        ndim = bertha.get_ndim()
        nshift = bertha.get_nshift()
        nocc = bertha.get_nocc()
        sfact = bertha.get_sfact()
        nopen = bertha.get_nopen()

        if stdoutprint:
            print("Verbosity       : ", verbosity)
            print("Dumpfiles       : ", dumpfiles)
            print("")
            print("Matrix dimension: ", ndim)
            print("            nocc: ", nocc)
            print("          nshift: ", nshift)
            print("           nopen: ", nopen)
            print("     level shift: ", sfact)
            print("")
    
        ovapm, eigem, fockm, eigen = bertha.run()

        if pberthaopt.dumpfiles :
            os.rename(vctfilename,'unpert_vct.txt') 

        if (fockm is None) or (eigen is None) or (fockm is None) \
                or (eigen is None):
            raise Exception("Error in bertha run")
    
        if stdoutprint:
            sys.stdout.flush()
            print("")
            print("Final results ")
            for i in range(nocc+nopen):
                print("eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact))
            print("      lumo       %20.8f"%(eigen[i+nshift+1]))
    
        erep = bertha.get_erep()
        etotal = bertha.get_etotal()
        ecoul = bertha.get_eecoul()
        exc = bertha.get_eexc()
    
        if stdoutprint:
            print("")
            print("total electronic energy  = %30.15f"%(etotal-(sfact*nocc)))
            print("nuclear repulsion energy = %30.15f"%(erep))
            print("total energy             = %30.15f"%(etotal+erep-(sfact*nocc)))
            print("coulomb energy           = %30.15f"%(ecoul))
            print("Exc     energy           = %30.15f"%(exc))
    else:
        if stdoutprint:
            print("Restart: ")

        my_file = Path(pberthaopt.restartfname)
        if not my_file.is_file():
            raise Exception(pberthaopt.restartfname + " does not exist restart=on")

        alldata = numpy.load(pberthaopt.restartfname) 

        eigem = alldata["eigem"]
        Da0 = alldata["Da0"]
        rho = alldata["rho"]
        etotalnpa = alldata["etotalnpa"]

        etotal = etotalnpa[0]
        dipx_ref = etotalnpa[1]
        dipy_ref = etotalnpa[2] 
        dipz_ref = etotalnpa[3] 

    activefname = pberthaopt.activefile
    if not os.path.isfile(activefname):
        raise Exception("File ", activefname , " does not exist")
    
    envirofname = pberthaopt.envirofile
    if not os.path.isfile(envirofname):
        raise Exception("File ", envirofname , " does not exist")

    embfactory = pyembmod.pyemb(activefname,envirofname,'adf') #jobtype='adf' is default de facto
    #grid_param =[50,110] # psi4 grid parameters (see Psi4 grid table)
    #embfactory.set_options(param=grid_param, \
    embfactory.set_options(param=pberthaopt.param, \
       gtype=pberthaopt.gtype, basis=pberthaopt.basis) 
    embfactory.set_enviro_func(pberthaopt.excfuncenv) 
    embfactory.set_thresh_conv(pberthaopt.thresh_conv)
    embfactory.set_grid_filename(pberthaopt.gridfname)
    # several paramenters to be specified in input- e.g AUG/ADZP for ADF, aug-cc-pvdz for psi4

    if (pberthaopt.debug):
       import uuid
       filename = "adfout_" + str(uuid.uuid4()) + ".txt"
       embfactory.set_adf_filenameout(filename)
       filename = "psi4out_" + str(uuid.uuid4()) + ".txt"
       embfactory.set_psi4_filenameout(filename)
   
    if stdoutprint:
        print("embfactory Options:")
        print(embfactory.get_options())

    embfactory.initialize()
    grid = embfactory.get_grid() 

    if (pberthaopt.debug):
       numpy.savetxt ("grid.txt", grid)

    if not restart:
        rho = bertha.get_density_on_grid(grid)

    if stdoutprint:
        print("Grid dimensions : " , grid.shape)
        print(" Rho dimensions : " , rho.shape)

    density=numpy.zeros((rho.shape[0],10))
    density[:,0] = rho

    fpot = None
    pot = embfactory.get_potential(density) 
    if (pberthaopt.debug):
        numpy.savetxt ("initialpot.txt", pot)
    #DEBUG
    #pot=numpy.zeros(rho.shape[0])
    if static_field:
      fpot = grid[:,fdir]*fmax 
      fpot = numpy.ascontiguousarray(fpot, dtype=numpy.double)
      if (pberthaopt.debug):
          numpy.savetxt ("initialfpot.txt", fpot)

    if not restart:

        if stdoutprint:
            #print("TEST density on grid")
            #print("Type density", type(density), density.shape)
            print("Scalar product" , "density.weigt", numpy.dot(density[:,0],grid[:,3]))
            print("Dip x" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,0]))
            print("Dip y" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,1]))
            print("Dip z" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,2]))

        bertha.realtime_init()

        normalise = 1

        dipx_mat, dipy_mat, dipz_mat = \
            bertha.get_realtime_dipolematrix (0, normalise)

        occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)
        iocc = 0

        for i in range(ndim):
            if i >= nshift and iocc < nocc:
                for j in range(ndim):
                    occeigv[j, iocc] = eigem[j, i]
                iocc = iocc + 1

        Da0 = numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))
  
        if stdoutprint:
            print("Dump ground state unperturbed density " + pberthaopt.densityzero)

        if pberthaopt.densityzero != "": 
            bertha.density_to_cube(Da0.T, pberthaopt.densityzero, \
                drx=pberthaopt.drx, dry=pberthaopt.dry, drz=pberthaopt.drz, \
                margin=pberthaopt.margin)

        dipx_ref = numpy.trace(numpy.matmul(Da0,dipx_mat)).real
        dipy_ref = numpy.trace(numpy.matmul(Da0,dipy_mat)).real
        dipz_ref = numpy.trace(numpy.matmul(Da0,dipz_mat)).real

        bertha.finalize()

        etotalnpa = numpy.zeros((4))
        etotalnpa[0] = etotal 
        etotalnpa[1] = dipx_ref
        etotalnpa[2] = dipy_ref
        etotalnpa[3] = dipz_ref

        numpy.savez(pberthaopt.restartfname, eigem=eigem, \
            Da0=Da0, rho=rho, etotalnpa=etotalnpa) 

    if stdoutprint:
        print("unperturbed Dip x    ",dipx_ref)
        print("unperturbed Dip y    ",dipy_ref)
        print("unperturbed Dip z    ",dipz_ref)

    #if lin_emb=True, a single scf is performed at constant Vemb
    maxiter = 20
    Dold = Da0 
    Da = Da0
    Eold = etotal
    lin_emb = pberthaopt.linemb

    dipx_val = 0.0
    dipy_val = 0.0
    dipz_val = 0.0

    dipx_mat = None
    dipy_mat = None
    dipz_mat = None

    for out_iter in range (maxiter):
        time_iter_start = time.time()
        process_time_iter_start = time.process_time()


        bertha.set_fnameinput(fnameinput)
        bertha.set_fittfname(fittfname)

        bertha.set_thresh(pberthaopt.thresh)
        
        bertha.set_verbosity(verbosity)
        bertha.set_dumpfiles(dumpfiles)
        
        bertha.set_densitydiff(1)
  
        bertha.init()

        ndim = bertha.get_ndim()
        nshift = bertha.get_nshift()
        nocc = bertha.get_nocc()
        sfact = bertha.get_sfact()
        nopen = bertha.get_nopen()
        
        # run with Vemb included
        if static_field:
           totpot = pot+fpot
           if nofde:
              totpot = fpot
           totpot = numpy.ascontiguousarray(totpot, dtype=numpy.double)

           if (pberthaopt.debug):
           	    numpy.savetxt ("fullpot%d.txt"%(out_iter), totpot)

           bertha.set_embpot_on_grid(grid, totpot)
        else:
           if (pberthaopt.debug):
                numpy.savetxt ("fullpot%d.txt"%(out_iter), pot)
           
           if nofde:
              pot=numpy.zeros_like(pot)
 
           bertha.set_embpot_on_grid(grid, pot)
        
        ovapm, eigem, fockm, eigen = bertha.run(eigem)
        etotal2 = bertha.get_etotal()

        # save the extern pot used in the scf run
        # pot_in = pot
        # the avg associated to pot
        emb_avg_in = numpy.dot(density[:,0]*grid[:,3],pot)
        #print

        if stdoutprint:
            print("")
            print("total electronic energy  = %30.15f ... outer iteration :%i"%((etotal2-(sfact*nocc)), (out_iter +1) ))
            print("mean value of embed pot  = %30.15f ... outer iteration :%i"%((emb_avg_in), (out_iter +1) ))

        erep = bertha.get_erep()
        etotal = bertha.get_etotal()
        ecoul = bertha.get_eecoul()
        exc = bertha.get_eexc()
        
        if stdoutprint:
            print("outer iteration : ", out_iter +1)
            print("total electronic energy  = %30.15f"%(etotal-(sfact*nocc)))
            print("nuclear repulsion energy = %30.15f"%(erep))
            print("total energy             = %30.15f"%(etotal+erep-(sfact*nocc)))
            print("coulomb energy           = %30.15f"%(ecoul))
            print("Exc     energy           = %30.15f"%(exc))

        if lin_emb :
            #rho = bertha.get_density_on_grid(grid)
            #density=numpy.zeros((rho.shape[0],10))
            #density[:,0] = rho
            occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)

            iocc = 0
            for i in range(ndim):
                if i >= nshift and iocc < nocc:
                    for j in range(ndim):
                        occeigv[j, iocc] = eigem[j, i]
                    iocc = iocc + 1
    
            Da = numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))
      
            if stdoutprint:
                print("Dump ground state perturbed density " + pberthaopt.density)

            if pberthaopt.density != "": 
                bertha.density_to_cube(Da.T, pberthaopt.density, drx=pberthaopt.drx, \
                    dry=pberthaopt.dry, drz=pberthaopt.drz, margin=pberthaopt.margin)

            bertha.realtime_init()

            normalise = 1

            dipx_mat, dipy_mat, dipz_mat = \
              bertha.get_realtime_dipolematrix (0, normalise)

            numofatom = bertha.get_natoms()
            tmpcoords = []
            tmpz      = []
    
            for i in range(numofatom):
    
                 tmpcoords.append(list(bertha.get_coords(i))[0:3])
                 tmpz.append(list(bertha.get_coords(i))[3])
    
            matcoords  = numpy.array(tmpcoords)
            vecZ       = numpy.array(tmpz)
            vec_nuc_dip= numpy.matmul(vecZ.T,matcoords)
    
            nucdipx = vec_nuc_dip[0]
            nucdipy = vec_nuc_dip[1]
            nucdipz = vec_nuc_dip[2]


            bertha.finalize()
            break

        #if out_iter == (maxiter - 1):
        #    bertha.finalize()
        #    raise Exception("Maximum number of SCF cycles reached.\n")

        # calculate the embedding potential corresponding to the new density

        time_get_density_on_grid_start = time.time()
        process_time_get_density_on_grid_start = time.process_time()

        rho = bertha.get_density_on_grid(grid)

        time_get_density_on_grid_end= time.time()
        process_time_get_density_on_grid_end= time.process_time()


        print("Fitted density on grid time: %15.5f"%(time_get_density_on_grid_end - time_get_density_on_grid_start), \
                    " (CPU time: %15.5f"%(process_time_get_density_on_grid_end - process_time_get_density_on_grid_start), ") s ")

        density=numpy.zeros((rho.shape[0],10))
        density[:,0] = rho
        pot_old=pot 

        time_start = time.time()
        process_time_start = time.process_time()

        pot = embfactory.get_potential(density)

        time_end= time.time()
        process_time_end= time.process_time()

        print("Vemb pot on grid: embfactory.get_potential: %15.5f"%(time_end - time_start), \
                    " (CPU time: %15.5f"%(process_time_end - process_time_start), ") s ")



        #DEBUG
        #pot=numpy.zeros_like(rho)
   
        # the avg associated to new  embed pot
        emb_avg = numpy.dot(density[:,0]*grid[:,3],pot)

        # form the density matrix
        occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)
        iocc = 0
        #debug
        if stdoutprint:
            print("ndim : %i\n" % ndim) 
            print("nshift : %i\n" % nshift) 
        
        for i in range(ndim):
            if i >= nshift and iocc < nocc:
                for j in range(ndim):
                    occeigv[j, iocc] = eigem[j, i]
                iocc = iocc + 1
 
        Da = numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))
        diffD = Da - Dold
        diffE = etotal2 -Eold
        norm_pot = numpy.sqrt(numpy.sum((pot-pot_old)**2))
        norm_D = numpy.linalg.norm(diffD,'fro') 
        
        if stdoutprint:
            print("2-norm of diffD  = ", norm_D, " ... outer iteration :%i"%(out_iter +1))
            print("E(actual)-E(prev)= ", diffE, " ... outer iteration :%i"%(out_iter +1))
            print("DE_emb(actual-prev)  = %30.15f ... outer iteration :%i"%((emb_avg_in-emb_avg), (out_iter +1) ))
            print("norm_pot(actual-prev)  = %30.15f ... outer iteration :%i"%((norm_pot), (out_iter +1) ))

        if ( norm_D<(1.0e-6) and norm_pot <(1.0e-4)):
            #not needed?
            #iocc = 0
            #for i in range(ndim):
            #    if i >= nshift and iocc < nocc:
            #        for j in range(ndim):
            #            occeigv[j, iocc] = eigem[j, i]
            #        iocc = iocc + 1
    
            #Da = numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))
      
            if stdoutprint:
                print("Dump ground state perturbed density density.cube")

            bertha.realtime_init()

            dipx_mat, dipy_mat, dipz_mat = bertha.get_realtime_dipolematrix (0, normalise)

            if pberthaopt.density != "":
               bertha.density_to_cube(Da.T, pberthaopt.density,  drx=pberthaopt.drx, \
                    dry=pberthaopt.dry, drz=pberthaopt.drz, margin=pberthaopt.margin)

            bertha.finalize()
            break

        time_iter_end= time.time()
        process_time_iter_end= time.process_time()

        print(" Time for each extrenal iteration: %15.5f"%(time_iter_end - time_iter_start), \
                    " (CPU time: %15.5f"%(process_time_iter_end - process_time_iter_start), ") s ")


 


        if out_iter == maxiter-1:
            if pberthaopt.density != "":
               bertha.density_to_cube(Da.T, pberthaopt.density,  drx=pberthaopt.drx, \
                    dry=pberthaopt.dry, drz=pberthaopt.drz, margin=pberthaopt.margin)

        Dold = Da
        Eold = etotal2

        bertha.realtime_init()

        normalise = 1

        dipx_mat, dipy_mat, dipz_mat = \
            bertha.get_realtime_dipolematrix (0, normalise)
        
        """
        Calculation of the Nuclear dipole (a.u)
        """
#####    print("Compute nuclear dipole (a.u.)")
    
        numofatom = bertha.get_natoms()
        tmpcoords = []
        tmpz      = []
    
        for i in range(numofatom):
    
             tmpcoords.append(list(bertha.get_coords(i))[0:3])
             tmpz.append(list(bertha.get_coords(i))[3])
    
        matcoords  = numpy.array(tmpcoords)
        vecZ       = numpy.array(tmpz)
        vec_nuc_dip= numpy.matmul(vecZ.T,matcoords)
    
        nucdipx = vec_nuc_dip[0]
        nucdipy = vec_nuc_dip[1]
        nucdipz = vec_nuc_dip[2]

        bertha.finalize()


    dipx_val = numpy.trace(numpy.matmul(Da,dipx_mat)).real
    dipy_val = numpy.trace(numpy.matmul(Da,dipy_mat)).real
    dipz_val = numpy.trace(numpy.matmul(Da,dipz_mat)).real

    if stdoutprint:
        print("")
        print("Dip x    ",dipx_val)
        print("Dip y    ",dipy_val)
        print("Dip z    ",dipz_val)

        print("MOLECULAR PROPERTIES of ACTIVE SYSTEM IN EMBEDDING:...")
        print("  ")
        print("Dipole Moment: ")
        print("Axis   electronic dipole (a.u.)    nuclear dipole (a.u.)      Total dipole (a.u.)      Total dipole (Debye) ")
        print(" x     %20.10f    %20.10f     %20.10f      %20.10f "%(dipx_val, nucdipx, dipx_val+nucdipx, (dipx_val+nucdipx)*2.541580))
        print(" y     %20.10f    %20.10f     %20.10f      %20.10f "%(dipy_val, nucdipy, dipy_val+nucdipy, (dipy_val+nucdipy)*2.541580))
        print(" z     %20.10f    %20.10f     %20.10f      %20.10f "%(dipz_val, nucdipz, dipz_val+nucdipz, (dipz_val+nucdipz)*2.541580))
        print("  ")
        tot_dip_vec = numpy.array([dipx_val+nucdipx,dipy_val+nucdipy,dipz_val+nucdipz])
        print("Molecular dipole module:  %20.10f (a.u.)  %20.10f (Debye) "%(numpy.linalg.norm(tot_dip_vec), numpy.linalg.norm(tot_dip_vec)*2.541580))

        print("MOLECULAR PROPERTIES of ACTIVE SYSTEM WITHOUT EMBEDDING (free active system):...")
        print("  ")
        print("Dipole Moment: ")
        print("Axis   electronic dipole (a.u.)    nuclear dipole (a.u.)      Total dipole (a.u.)      Total dipole (Debye) ")
        print(" x     %20.10f    %20.10f     %20.10f      %20.10f "%(dipx_ref, nucdipx, dipx_ref+nucdipx, (dipx_ref+nucdipx)*2.541580))
        print(" y     %20.10f    %20.10f     %20.10f      %20.10f "%(dipy_ref, nucdipy, dipy_ref+nucdipy, (dipy_ref+nucdipy)*2.541580))
        print(" z     %20.10f    %20.10f     %20.10f      %20.10f "%(dipz_ref, nucdipz, dipz_ref+nucdipz, (dipz_ref+nucdipz)*2.541580))
        print("  ")
        tot_dip_vec = numpy.array([dipx_ref+nucdipx,dipy_ref+nucdipy,dipz_ref+nucdipz])
        print("Molecular dipole module:  %20.10f (a.u.)  %20.10f (Debye) "%(numpy.linalg.norm(tot_dip_vec), numpy.linalg.norm(tot_dip_vec)*2.541580))

    if stdoutprint:
        print("Dipole moment analitical: Tr(D dip_mat)")
        print("dipX_mat dim: %i,%i\n" % (dipx_mat.shape))
        print("dipY_mat dim: %i,%i\n" % (dipy_mat.shape))
        print("dipZ_mat dim: %i,%i\n" % (dipz_mat.shape))
        print("Da dim: %i,%i\n" % (Da.shape))

        print("TEST dipole moment from density on grid numerical integration")
        print("  ")
        #print("Type density", type(density), density.shape)
        print("Scalar product" , "density.weigt", numpy.dot(density[:,0],grid[:,3]))
        print("Dip x" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,0]))
        print("Dip y" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,1]))
        print("Dip z" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,2]))
        
        print("  ")
        sys.stdout.flush()
        print("Final results of SplitSCF")
        for i in range(nocc+nopen):
            print("eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact))
        print("      lumo       %20.8f"%(eigen[i+nshift+1]))

    return ovapm, eigem, fockm, eigen, pot 

##########################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-gA","--geom_act", help="Specify active system (Angstrom) geometry (default: geomA.xyz)", required=False, 
            type=str, default="geomA.xyz")
    parser.add_argument("-gB","--geom_env", help="Specify frozen system (Angstrom) geometry (default: geomB.xyz)", required=False, 
            type=str, default="geomB.xyz")
    parser.add_argument("--embthresh", help="set EMB threshold (default = 1.0e-8)", required=False, 
            type=numpy.float64, default=1.0e-8)
    parser.add_argument("-d", "--debug", help="Debug on, prints debug info to debug_info.txt", required=False, 
            default=False, action="store_true")
    parser.add_argument("--modpaths", help="set berthamod and all other modules path [\"path1;path2;...\"] (default = ../src)", 
        required=False, type=str, default="../src")

    parser.add_argument("--act_obs", \
        help="Specify BERTHA (Active system) basisset \"atomname1:basisset1,atomname2:basisset2,...\"", \
        required=True, type=str, default="")
    parser.add_argument("--act_fittset", \
        help="Specify BERTHA (Active system) fitting set \"atomname1:fittset1,atomname2:fittset2,...\"", \
        required=True, type=str, default="")
    parser.add_argument("--act_func", 
	help="Specify exchange–correlation energy functional for active system available: LDA,B88P86,HCTH93,BLYP (default=BLYP)", \
        type=str, default="BLYP")
    parser.add_argument("--env_obs", help="Specify the orbital basis set for the enviroment (default: AUG/ADZP)", required=False,
            type=str, default="AUG/ADZP")
    parser.add_argument("--env_func", help="Specify the function for the environment density (default: BLYP)", required=False,
            type=str, default="BLYP")
    parser.add_argument("--grid_opts", help="set gridtype (default: 2)",
        required=False, type=int, default=2)
    parser.add_argument("--gridfname", help="set grid filename (default = grid.dat)",
        required=False, type=str, default="grid.dat")
    parser.add_argument("-l", "--linemb", help="Linearized embedding on: the outer loop is skipped", required=False,
            default=False, action="store_true")
    parser.add_argument("--nofde", help="embedding off: just for test", required=False, 
            default=False, action="store_true")
    parser.add_argument("--static_field", help="Add a static field to the SCF (default : False)", required=False, 
            default=False, action="store_true")
    parser.add_argument("--fmax", help="Static field amplitude (default : 1.0e-5)", required=False, 
            type=numpy.float64, default=1.0e-5)
    parser.add_argument("--fdir", help="External field direction (cartesian)  (default: 2)",
            required=False, type=int, default=2)
    parser.add_argument("--restart", help="Force restart from a previous initial single-point",
        required=False, action="store_true", default=False)
    parser.add_argument("-c","--fitcoefffile", help="Specify BERTHA fitcoeff output file (default: fitcoeff.txt)",
        required=False, type=str, default="fitcoeff.txt")
    parser.add_argument("-e","--vctfile", help="Specify BERTHA vct output file (default: vct.txt)", required=False, 
        type=str, default="vct.txt")
    parser.add_argument("-p","--ovapfile", help="Specify BERTHA ovap output file (default: ovap.txt)", required=False, 
        type=str, default="ovap.txt")
    parser.add_argument("-s", "--dumpfiles", help="Dumpfile on, default is off", required=False,
        default=False, action="store_true")
    parser.add_argument("-v", "--verbosity", help="Verbosity level 0 = minim, -1 = print iteration info, " + 
        "1 = maximum (defaul -1)", required=False, default=-1, type=int)
    parser.add_argument("--thresh", help="set bertha threshold (default = 1.0e-11)", required=False, 
        type=numpy.float64, default=1.0e-11)
    parser.add_argument("--wrapperso", help="set wrapper SO (default = ../lib/bertha_wrapper.so)", 
        required=False, type=str, default="../lib/bertha_wrapper.so")
    parser.add_argument("--eda_nocv_info", help="set to dump info useful for py_eda_nocv",action='store_true',default=False)
    parser.add_argument("--eda_nocv_frag_file", help="set a file (default: info_eda_nocv_fragX.json)",
        required=False, type=str, default="info_eda_nocv_fragX.json")
    parser.add_argument("-j","--jsonbasisfile", \
        help="Specify BERTHA JSON file for fitting and basis (default: fullsets.json)", \
       required=False, type=str, default="fullsets.json")
    parser.add_argument("--convertlengthunit", help="Specify a length converter [default=1.0]", \
        type=float, default=1.0)
    parser.add_argument("--berthamaxit", help="set bertha maxiterations (default = %d)"%(MAXIT), 
        required=False, type=int, default=MAXIT)

    parser.add_argument("--cube", help="Specify if density need to be saved in cube file format (default: 0)", required=False,
               default=False, action="store_true")
    parser.add_argument("--deltax", help="cube dx step (default = 0.2)", required=False,
            type=numpy.float64, default=0.2)
    parser.add_argument("--deltay", help="cube dy step (default = 0.2)", required=False,
            type=numpy.float64, default=0.2)
    parser.add_argument("--deltaz", help="cube dz step (default = 0.2)", required=False,
            type=numpy.float64, default=0.2)
    parser.add_argument("--lmargin", help="cube margin parameter (default = 10.0)", required=False,
            type=numpy.float64, default=10.0)

   
    args = parser.parse_args()
  
    for path in args.modpaths.split(";"):
        sys.path.append(path)

    modpaths = os.environ.get('PYBERTHA_MOD_PATH')

    if modpaths is not None :
        for path in modpaths.split(";"):
            sys.path.append(path)

    for resdir in ["./resultfiles", "./jobtempdir"]:
      if os.path.isdir(resdir):
         print ("  Removing "+ resdir )
         shutil.rmtree(resdir)

    import pybgen

    pygenoption_fraga = pybgen.berthainputoption
    pygenoption_fraga.inputfile = args.geom_act
    pygenoption_fraga.jsonbasisfile = args.jsonbasisfile
    pygenoption_fraga.fittset = args.act_fittset
    pygenoption_fraga.basisset = args.act_obs
    pygenoption_fraga.functxc = args.act_func
    pygenoption_fraga.convertlengthunit = args.convertlengthunit
    pygenoption_fraga.maxit = args.berthamaxit

    for filename in ["input.inp", "fitt2.inp"]:
      try:
        os.remove(filename)
      except OSError:
        pass

    pybgen.generateinputfiles (pygenoption_fraga)

    pberthaopt = pyberthaembedoption

    pberthaopt.gridfname = args.gridfname

    pberthaopt.fitcoefffile = args.fitcoefffile
    pberthaopt.vctfile = args.vctfile
    pberthaopt.ovapfile = args.ovapfile
    pberthaopt.dumpfiles = args.dumpfiles
    pberthaopt.debug = args.debug
    pberthaopt.linemb = args.linemb
    pberthaopt.nofde = args.nofde
    pberthaopt.verbosity = args.verbosity
    pberthaopt.thresh = args.thresh
    pberthaopt.static_field = args.static_field
    pberthaopt.fmax = args.fmax
    pberthaopt.fdir = args.fdir
    pberthaopt.wrapperso = args.wrapperso
    pberthaopt.eda_nocv_info = args.eda_nocv_info
    pberthaopt.eda_nocv_frag_file = args.eda_nocv_frag_file
    pberthaopt.activefile = args.geom_act
    pberthaopt.envirofile = args.geom_env
    pberthaopt.gtype = args.grid_opts
    pberthaopt.thresh_conv = args.embthresh
    pberthaopt.basis = args.env_obs
    pberthaopt.excfuncenv = args.env_func

    if args.cube:
        pberthaopt.densityzero = "density0.cube"
        pberthaopt.density = "density.cube"
        pberthaopt.drx = args.deltax 
        pberthaopt.dry = args.deltay
        pberthaopt.drz = args.deltaz
        pberthaopt.margin = args.lmargin

    runspberthaembed (pberthaopt, args.restart)
