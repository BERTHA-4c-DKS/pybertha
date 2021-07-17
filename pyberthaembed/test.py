import argparse
import ctypes
from os import X_OK
import numpy
import sys
import re

import os.path
sys.path.append('/home/matteod/pybertha/pyemb')
sys.path.append("/home/matteod/build/xcfun/build/lib/python")
sys.path.append("/home/matteod/pybertha/psi4rt")
sys.path.append("/home/matteod/pybertha/src")
sys.path.append("/home/matteod/build/pyadf/src")

os.environ['PYBERTHAROOT'] = "/home/matteod/pybertha/"
os.environ['RTHOME'] = "/home/matteod/pybertha/psi4rt"
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
    verbosity: int
    thresh: numpy.float64
    wrapperso: str
    berthamodpath: str
    eda_nocv_info: bool
    eda_nocv_frag_file: str

##########################################################################################

def runspberthaembed (pberthaopt):

    sys.path.insert(0, pberthaopt.berthamodpath)
    import berthamod
    
    print("Options: ")
    for att in [a for a in dir(pberthaopt) if not a.startswith('__')]:
        print(att, " = ", getattr(pberthaopt, att)) 
    print("")
    print("")

    print("")
    
    if not os.path.isfile(pberthaopt.wrapperso):
        print("SO File ", pberthaopt.wrapperso, " does not exist")
        exit(1)
    
    bertha = berthamod.pybertha(pberthaopt.wrapperso)
    
    fittcoefffname = pberthaopt.fitcoefffile
    vctfilename = pberthaopt.vctfile
    ovapfilename = pberthaopt.ovapfile
    
    fnameinput = pberthaopt.inputfile
    if not os.path.isfile(fnameinput):
        print("File ", fnameinput, " does not exist")
        exit(1)
    
    fittfname = pberthaopt.fittfile
    if not os.path.isfile(fittfname):
        print("File ", fittfname , " does not exist")
        exit(1)
    
    verbosity = pberthaopt.verbosity
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
    
    print("Verbosity       : ", verbosity)
    print("Dumpfiles       : ", dumpfiles)
    print("")
    print("Matrix dimension: ", ndim)
    print("            nocc: ", nocc)
    print("          nshift: ", nshift)
    print("           nopen: ", nopen)
    print("     level shift: ", sfact)
    print("")
    

    #generate the unperturbed ground state density
    """
    npoints = 10
    grid = numpy.zeros((npoints, 4))
    grid = numpy.ascontiguousarray(grid, dtype=numpy.double)
    pot = numpy.zeros(npoints)
    pot = numpy.ascontiguousarray(pot, dtype=numpy.double)

    x = -1.0
    y = -100.0
    z = -1000.0
    w = 1.0
    for i in range(npoints):
        grid[i,0] = x
        grid[i,1] = y
        grid[i,2] = z
        grid[i,3] = w

        pot[i] = x*y*z*w

        x += 1.0
        y += 1.0
        z += 1.0
        w += 0.1
    """

    start = time.time()
    cstart = time.process_time() 

    
#main run here

    ovapm, eigem, fockm, eigen = bertha.run()
    
    end = time.time()
    cend = time.process_time()

#   initialize pyembed instance

    activefname = pberthaopt.activefile
    if not os.path.isfile(activefname):
        print("File ", activefname , " does not exist")
        exit(1)
    
    envirofname = pberthaopt.envirofile
    if not os.path.isfile(envirofname):
        print("File ", envirofname , " does not exist")
        exit(1)

    embfactory = pyembed.pyemb(activefname,envirofname,'adf') #jobtype='adf' is default de facto
    embfactory.set_options(4.0)  # only the integration accuracy parameter is specified 
   
    embfactory.print_options()

    embfactory.initialize()
    grid = embfactory.get_grid() 
    
    #DEBUG : quick check of grid
    print("Type grid", type(grid), grid.shape)
    
    density = bertha.get_density_on_grid(grid)
    
    pot = embfactory.get_potential(density)    

#  TEST density on grid

#    print("TEST density on grid")
#    print("Type density", type(density), density.shape)
#    print("Scalar product" , "density.weigt", numpy.dot(density,grid[:,3]))
#    print("Dip x" , "density.weigt", numpy.dot(density*grid[:,3],grid[:,0]))
#    print("Dip y" , "density.weigt", numpy.dot(density*grid[:,3],grid[:,1]))
#    print("Dip z" , "density.weigt", numpy.dot(density*grid[:,3],grid[:,2]))

    """
    for i in range(npoints):
        val = grid[i,0] * grid[i,1]  * grid[i,2] * grid[i,3]
        print("Python L: %15.5f vs %15.5f"%(density[i], val))
    """


    print("Totaltime:    ", end - start, " (CPU time: " , cend - cstart, ") s ")
    print("MainRun Time: ", bertha.get_mainruntime() , \
            " (CPU time: ", bertha.get_mainrunctime(), ") s ")

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

    Da = numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))

    dipx_ref = numpy.trace(numpy.matmul(Da,dipx_mat)).real
    dipy_ref = numpy.trace(numpy.matmul(Da,dipy_mat)).real
    dipz_ref = numpy.trace(numpy.matmul(Da,dipz_mat)).real

    print("Dipx    ",dipx_ref)
    print("Dipy    ",dipy_ref)
    print("Dipz    ",dipz_ref)




    bertha.finalize()

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

#   run with Vemb included
    bertha.set_embpot_on_grid(grid, pot)
    ovapm, eigem, fockm, eigen = bertha.run()
    etotal2 = bertha.get_etotal()

    dipx = (etotal2 - etotal1)/(2*field)


    print("Dipole moment analitical: Tr(D dip_mat)")

    print("Dip x    ",dipx_ref)
    print("Dip y    ",dipy_ref)
    print("Dip z    ",dipz_ref)

    print("TEST dipole moment from density on grid numerical integration")
    print("  ")
    print("Type density", type(density), density.shape)
    print("Scalar product" , "density.weigt", numpy.dot(density,grid[:,3]))
    print("Dip x" , "density.weigt", numpy.dot(density*grid[:,3],grid[:,0]))
    print("Dip y" , "density.weigt", numpy.dot(density*grid[:,3],grid[:,1]))
    print("Dip z" , "density.weigt", numpy.dot(density*grid[:,3],grid[:,2]))

    print("  ")
    print("TEST one electron potential (embedding) in g-spinor")
    print("Potential used -x*E and x*E evaluating dipole via finite field")
    print("Dip x" , "from finite field", dipx)
    
    sys.stdout.flush()

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
    
    args = parser.parse_args()

    pberthaopt = pyberthaembedoption

    pberthaopt.inputfile = args.inputfile
    pberthaopt.fittfile = args.fittfile
    pberthaopt.fitcoefffile = args.fitcoefffile
    pberthaopt.vctfile = args.vctfile
    pberthaopt.ovapfile = args.ovapfile
    pberthaopt.dumpfiles = args.dumpfiles
    pberthaopt.debug = args.debug
    pberthaopt.verbosity = args.verbosity
    pberthaopt.thresh = args.thresh
    pberthaopt.wrapperso = args.wrapperso
    pberthaopt.berthamodpath = args.berthamodpath
    pberthaopt.eda_nocv_info = args.eda_nocv_info
    pberthaopt.eda_nocv_frag_file = args.eda_nocv_frag_file
    pberthaopt.activefile = args.geomA
    pberthaopt.envirofile = args.geomB

    runspberthaembed (pberthaopt)
