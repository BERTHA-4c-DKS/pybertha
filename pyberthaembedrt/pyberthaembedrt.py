import argparse
import ctypes
from os import X_OK
import numpy
import sys
import re

import os.path

from numpy.linalg import eigvalsh
from scipy.linalg import eigh

import json
from json import encoder

import time

from dataclasses import dataclass

@dataclass
class pyberthaembedrtoption:
    inputfile: str
    fittfile: str
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

def runspberthaembedrt (pberthaopt):

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
    

    #generate a sample grid
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

        x += 1.0
        y += 1.0
        z += 1.0
        w += 0.1

    start = time.time()
    cstart = time.process_time() 

    bertha.set_embpot_on_grid(grid, pot)
    
    #main run here
    ovapm, eigem, fockm, eigen = bertha.run()

    density = bertha.get_density_on_grid(grid)

    for i in range(npoints):
        val = grid[i,0] * grid[i,1]  * grid[i,2] * grid[i,3]
        print("Python L: %15.5f vs %15.5f"%(density[i], val))

    end = time.time()
    cend = time.process_time()

    print("Totaltime:    ", end - start, " (CPU time: " , cend - cstart, ") s ")
    print("MainRun Time: ", bertha.get_mainruntime() , \
            " (CPU time: " , bertha.get_mainrunctime(), ") s ")
    
    sys.stdout.flush()

##########################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--inputfile", help="Specify BERTHA input file (default: input.inp)", required=False, 
            type=str, default="input.inp")
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

    pberthaopt = pyberthaembedrtoption

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

    runspberthaembedrt (pberthaopt)
