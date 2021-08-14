import argparse
import ctypes
from os import X_OK
import numpy
import sys
import re

import os.path
sys.path.append('/home/redo/BERTHA/pybertha/pyemb')
sys.path.append("/home/redo/BERTHA/xcfun/build/lib/python")
sys.path.append("/home/redo/BERTHA/pybertha/src")
sys.path.append("/home/redo/BERTHA/pyadf/src")
os.environ['PYBERTHAROOT'] = "/home/redo/BERTHA/pybertha/"
os.environ['RTHOME'] = "/home/redo/BERTHA/pybertha/psi4rt"
sys.path.append(os.environ['PYBERTHAROOT']+"/src")
sys.path.append(os.environ['RTHOME'])

import pyembmod

from numpy.linalg import eigvalsh
from scipy.linalg import eigh

import json
from json import encoder

import time

from dataclasses import dataclass

@dataclass
class pyberthaembedoption:
    inputfile: str
    fittfile: str
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
    gtype: int = 2
    param: float = 4.0
    basis: str = 'AUG/ADZP'
    denistyzero: str = "density0.cube"
    density : str = "density.cube"
    drx: float = 0.1
    dry: float = 0.1
    drz: float = 0.1
    margin: float = 5.5

##########################################################################################

def runspberthaembed (pberthaopt, stdoutprint = True):

    ovapm = None 
    eigem = None 
    fockm = None  
    eigen = None

    verbosity = -1
    verbosity = pberthaopt.verbosity

    if not stdoutprint:
        verbosity = 0

    sys.path.insert(0, pberthaopt.berthamodpath)
    import berthamod
    
    if stdoutprint:
    	print("Options: ")
    	for att in [a for a in dir(pberthaopt) if not a.startswith('__')]:
        	print(att, " = ", getattr(pberthaopt, att)) 
    	print("")
    	print("")

    	print("")
    
    if not os.path.isfile(pberthaopt.wrapperso):
        raise Exception("SO File ", pberthaopt.wrapperso, " does not exist")

    bertha = berthamod.pybertha(pberthaopt.wrapperso)
    
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
    

    #start = time.time()
    #cstart = time.process_time() 

    ovapm, eigem, fockm, eigen = bertha.run()

    if pberthaopt.dumpfiles :
       os.rename(vctfilename,'unpert_vct.txt') 

    #end = time.time()
    #cend = time.process_time()

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
    ecoul  = bertha.get_eecoul()
    exc    = bertha.get_eexc()
    
    if stdoutprint:
        print("")
        print("total electronic energy  = %30.15f"%(etotal-(sfact*nocc)))
        print("nuclear repulsion energy = %30.15f"%(erep))
        print("total energy             = %30.15f"%(etotal+erep-(sfact*nocc)))
        print("coulomb energy           = %30.15f"%(ecoul))
        print("Exc     energy           = %30.15f"%(exc))
    
    #initialize pyembed instance

    activefname = pberthaopt.activefile
    if not os.path.isfile(activefname):
        raise Exception("File ", activefname , " does not exist")
    
    envirofname = pberthaopt.envirofile
    if not os.path.isfile(envirofname):
        raise Exception("File ", envirofname , " does not exist")

    embfactory = pyembmod.pyemb(activefname,envirofname,'adf') #jobtype='adf' is default de facto
    #grid_param =[50,110] # psi4 grid parameters (see Psi4 grid table)
    embfactory.set_options(param=pberthaopt.param, \
       gtype=pberthaopt.gtype, basis=pberthaopt.basis) 
    # several paramenters to be specified in input- e.g AUG/ADZP for ADF, aug-cc-pvdz for psi4
   
    if stdoutprint:
        print("embfactory Options:")
        print(embfactory.get_options())

    embfactory.initialize()
    grid = embfactory.get_grid() 
    
    #DEBUG : quick check of grid
    #if stdoutprint:
    #    print("Type grid", type(grid), grid.shape)
    
    rho = bertha.get_density_on_grid(grid)
    density=numpy.zeros((rho.shape[0],10))
    density[:,0] = rho
   
    pot = embfactory.get_potential(density)    

    #TEST density on grid

    if stdoutprint:
        #print("TEST density on grid")
        #print("Type density", type(density), density.shape)
        print("Scalar product" , "density.weigt", numpy.dot(density[:,0],grid[:,3]))
        print("Dip x" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,0]))
        print("Dip y" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,1]))
        print("Dip z" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,2]))

    bertha.realtime_init()

    normalise = 1
    dip_mat = None

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
        print("Dump ground state unperturbed density " + pberthaopt.denistyzero)
    
    bertha.density_to_cube(Da0.T, pberthaopt.denistyzero, \
        drx=pberthaopt.drx, dry=pberthaopt.dry, drz=pberthaopt.drz, \
        margin=pberthaopt.margin)

    dipx_ref = numpy.trace(numpy.matmul(Da0,dipx_mat)).real
    dipy_ref = numpy.trace(numpy.matmul(Da0,dipy_mat)).real
    dipz_ref = numpy.trace(numpy.matmul(Da0,dipz_mat)).real

    if stdoutprint:
        print("unperturbed Dip x    ",dipx_ref)
        print("unperturbed Dip y    ",dipy_ref)
        print("unperturbed Dip z    ",dipz_ref)

    bertha.finalize()

    #bertha.set_fittcoefffname(fittcoefffname)
    #bertha.set_ovapfilename(ovapfilename)
    #bertha.set_vctfilename(vctfilename)
    
    #if lin_emb=True, a single scf is performed at constant Vemb
    maxiter = 10 
    Dold = Da0 
    Eold = etotal
    lin_emb = pberthaopt.linemb

    for out_iter in range (maxiter):

        bertha.set_fnameinput(fnameinput)
        bertha.set_fittfname(fittfname)

        bertha.set_thresh(pberthaopt.thresh)
        
        bertha.set_verbosity(verbosity)
        bertha.set_dumpfiles(dumpfiles)
        
        bertha.set_densitydiff(1)
  
        bertha.init()
        
        # run with Vemb included
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
        ecoul  = bertha.get_eecoul()
        exc    = bertha.get_eexc()
        
        if stdoutprint:
            print("outer iteration : ", out_iter +1)
            print("total electronic energy  = %30.15f"%(etotal-(sfact*nocc)))
            print("nuclear repulsion energy = %30.15f"%(erep))
            print("total energy             = %30.15f"%(etotal+erep-(sfact*nocc)))
            print("coulomb energy           = %30.15f"%(ecoul))
            print("Exc     energy           = %30.15f"%(exc))

        if lin_emb :

            rho = bertha.get_density_on_grid(grid)
            density=numpy.zeros((rho.shape[0],10))
            density[:,0] = rho
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
            
            bertha.density_to_cube(Da.T, pberthaopt.density, drx=pberthaopt.drx, \
              dry=pberthaopt.dry, drz=pberthaopt.drz, margin=pberthaopt.margin)
            bertha.finalize()
            break

        if out_iter == maxiter:
            bertha.finalize()
            raise Exception("Maximum number of SCF cycles exceeded.\n")

        # calculate the embedding potential corresponding to the new density

        rho = bertha.get_density_on_grid(grid)
        density=numpy.zeros((rho.shape[0],10))
        density[:,0] = rho
       
        pot = embfactory.get_potential(density) 
   
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
        norm_D=numpy.linalg.norm(diffD,'fro') 
        
        if stdoutprint:
            print("2-norm of diffD  = ", norm_D, " ... outer iteration :%i"%(out_iter +1))
            print("E(actual)-E(prev)= ", diffE, " ... outer iteration :%i"%(out_iter +1))
            print("DE_emb(actual-new)  = %30.15f ... outer iteration :%i"%((emb_avg_in-emb_avg), (out_iter +1) ))

        if ( norm_D<(1.0e-3) and diffE <(1.0e-6)):
            iocc = 0
            for i in range(ndim):
                if i >= nshift and iocc < nocc:
                    for j in range(ndim):
                        occeigv[j, iocc] = eigem[j, i]
                    iocc = iocc + 1
    
            Da = numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))
      
            if stdoutprint:
                print("Dump ground state perturbed density density.cube")

            bertha.density_to_cube(Da.T, "density.cube",drx=0.1,dry=0.1,drz=0.1,margin=5.5)
            bertha.finalize()
            break

        Dold = Da
        Eold = etotal2

        bertha.finalize()

    if stdoutprint:
        print("Dipole moment analitical: Tr(D dip_mat)")

    dipx_val = numpy.trace(numpy.matmul(Da,dipx_mat)).real
    dipy_val = numpy.trace(numpy.matmul(Da,dipy_mat)).real
    dipz_val = numpy.trace(numpy.matmul(Da,dipz_mat)).real

    if stdoutprint:
        print("Dip x    ",dipx_val)
        print("Dip y    ",dipy_val)
        print("Dip z    ",dipz_val)
        
        print("TEST dipole moment from density on grid numerical integration")
        print("  ")
        #print("Type density", type(density), density.shape)
        print("Scalar product" , "density.weigt", numpy.dot(density[:,0],grid[:,3]))
        print("Dip x" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,0]))
        print("Dip y" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,1]))
        print("Dip z" , "density.weigt", -1.*numpy.dot(density[:,0]*grid[:,3],grid[:,2]))
        
        print("  ")
        sys.stdout.flush()

    return ovapm, eigem, fockm, eigen 

##########################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--inputfile", help="Specify BERTHA input file (default: input.inp)", required=False, 
            type=str, default="input.inp")
    parser.add_argument("-gA","--geomA", help="Specify active system (Angstrom) geometry (default: geomA.xyz)", required=False, 
            type=str, default="geomA.xyz")
    parser.add_argument("-gB","--geomB", help="Specify frozen system (Angstrom) geometry (default: geomB.xyz)", required=False, 
            type=str, default="geomB.xyz")
    parser.add_argument("-t","--fittfile", help="Specify BERTHA fitting input file (default: fitt2.inp)", required=False, 
            type=str, default="fitt2.inp")
    parser.add_argument("-c","--fitcoefffile", help="Specify BERTHA fitcoeff output file (default: fitcoeff.txt)",
            required=False, type=str, default="fitcoeff.txt")
    parser.add_argument("-e","--vctfile", help="Specify BERTHA vct output file (default: vct.txt)", required=False, 
            type=str, default="vct.txt")
    parser.add_argument("-p","--ovapfile", help="Specify BERTHA ovap output file (default: ovap.txt)", required=False, 
            type=str, default="ovap.txt")
    parser.add_argument("-s", "--dumpfiles", help="Dumpfile on, default is off", required=False,
            default=False, action="store_true")
    parser.add_argument("-d", "--debug", help="Debug on, prints debug info to debug_info.txt", required=False, 
            default=False, action="store_true")
    parser.add_argument("-l", "--linemb", help="Linearized embedding on: the outer loop is skipped", required=False, 
            default=False, action="store_true")
    parser.add_argument("-v", "--verbosity", help="Verbosity level 0 = minim, -1 = print iteration info, " + 
            "1 = maximum (defaul -1)", required=False, default=-1, type=int)
    parser.add_argument("--thresh", help="det threshold (default = 1.0e-11)", required=False, 
            type=numpy.float64, default=1.0e-11)
    parser.add_argument("--wrapperso", help="set wrapper SO (default = ../../lib/bertha_wrapper.so)", 
            required=False, type=str, default="../lib/bertha_wrapper.so")
    parser.add_argument("--berthamodpath", help="set berthamod path (default = ../src)", 
            required=False, type=str, default="../src")
    parser.add_argument("--eda_nocv_info", help="set to dump info useful for py_eda_nocv",action='store_true',default=False)
    parser.add_argument("--eda_nocv_frag_file", help="set a file (default: info_eda_nocv_fragX.json)",
            required=False, type=str, default="info_eda_nocv_fragX.json")
    parser.add_argument("--gridtype", help="set gridtype (default: 2)",
            required=False, type=int, default=2)
    
    args = parser.parse_args()

    pberthaopt = pyberthaembedoption

    pberthaopt.inputfile = args.inputfile
    pberthaopt.fittfile = args.fittfile
    pberthaopt.fitcoefffile = args.fitcoefffile
    pberthaopt.vctfile = args.vctfile
    pberthaopt.ovapfile = args.ovapfile
    pberthaopt.dumpfiles = args.dumpfiles
    pberthaopt.debug = args.debug
    pberthaopt.linemb = args.linemb
    pberthaopt.verbosity = args.verbosity
    pberthaopt.thresh = args.thresh
    pberthaopt.wrapperso = args.wrapperso
    pberthaopt.berthamodpath = args.berthamodpath
    pberthaopt.eda_nocv_info = args.eda_nocv_info
    pberthaopt.eda_nocv_frag_file = args.eda_nocv_frag_file
    pberthaopt.activefile = args.geomA
    pberthaopt.envirofile = args.geomB
    pberthaopt.gtype = args.gridtype

    runspberthaembed (pberthaopt)
