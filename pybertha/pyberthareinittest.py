import argparse
import ctypes
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
class pyberthaoption:
    inputfile: str = "input.inp"
    fittfile: str = "fitt2.inp"
    fitcoefffile: str = "fitcoeff.txt"
    vctfile: str ="vct.txt"
    ovapfile: str ="ovap.txt"
    dumpfiles: bool =False
    debug: bool =False
    verbosity: int =-1
    thresh: numpy.float64 = 1.0e-11
    wrapperso: str = "../lib/bertha_wrapper.so"
    berthamodpath: str ="../src"
    eda_nocv_info: bool =False
    eda_nocv_frag_file: str = "info_eda_nocv_fragX.json"


##########################################################################################

def get_json_data(pberthaopt, etotal, erep, ecoul, exc, ndim, nocc, occeigv):

    json_data = {}
    for arg in vars(pberthaopt):
        #print(type(arg), arg, getattr(pberthaopt, arg))
        if arg.find("__") < 0:
            json_data[arg] = getattr(pberthaopt, arg)

    othervals = {
            'etotal': etotal,
            'erep'  : erep, 
            'ecoul' : ecoul, 
            'exc'  : exc, 
            'ndim' : ndim,
            'nocc' : nocc, 
            'occeigv_REAL' : numpy.real(occeigv).tolist(),  
            'occeigv_IMAG' : numpy.imag(occeigv).tolist() 
            }

    json_data.update(othervals)

    return json_data

##########################################################################################

def runspbertha (pberthaopt):

    sys.path.insert(0, pberthaopt.berthamodpath)
    import berthamod
    
    print("Options: ")
    for att in [a for a in dir(pberthaopt) if not a.startswith('__')]:
        print(att, " = ", getattr(pberthaopt, att)) 
    print("")
    print("")
    
    if not os.path.isfile(pberthaopt.wrapperso):
        print("SO File ", pberthaopt.wrapperso, " does not exist")
        exit(1)

    verbosity = pberthaopt.verbosity
    dumpfiles = int(pberthaopt.dumpfiles)
 
    
    bertha = berthamod.pybertha(pberthaopt.wrapperso)

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
    
    start = time.time()
    cstart = time.process_time() 
    
    ovapm, eigem, fockm, eigen = bertha.run()
    
    end = time.time()
    cend = time.process_time()
    
    print("Totaltime:    ", end - start, " (CPU time: " , cend - cstart, ") s ")
    print("MainRun Time: ", bertha.get_mainruntime() , \
            " (CPU time: " , bertha.get_mainrunctime(), ") s ")
    
    sys.stdout.flush()
    
    
    if (fockm is None) or (eigen is None) or (fockm is None) \
            or (eigen is None):
        print("Error in bertha run")
        exit(-1)
    
    
    print("")
    print("Final results ")
    for i in range(nocc+nopen):
        print("eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact))
        
    print("      lumo       %20.8f"%(eigen[i+nshift+1]))
    
    erep = bertha.get_erep()
    etotal = bertha.get_etotal()
    ecoul  = bertha.get_eecoul()
    exc    = bertha.get_eexc()
    
    
    print("")
    print("total electronic energy  = %30.15f"%(etotal-(sfact*nocc)))
    print("nuclear repulsion energy = %30.15f"%(erep))
    print("total energy             = %30.15f"%(etotal+erep-(sfact*nocc)))
    print("coulomb energy           = %30.15f"%(ecoul))
    print("Exc     energy           = %30.15f"%(exc))
    
    
    bertha.finalize()

    print("\n\nSTART second RUN \n\n")

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
    
    start = time.time()
    cstart = time.process_time() 
    
    ovapm, eigem, fockm, eigen = bertha.run()
    
    end = time.time()
    cend = time.process_time()
    
    print("Totaltime:    ", end - start, " (CPU time: " , cend - cstart, ") s ")
    print("MainRun Time: ", bertha.get_mainruntime() , \
            " (CPU time: " , bertha.get_mainrunctime(), ") s ")
    
    sys.stdout.flush()
    
    
    if (fockm is None) or (eigen is None) or (fockm is None) \
            or (eigen is None):
        print("Error in bertha run")
        exit(-1)
    
    
    print("")
    print("Final results ")
    for i in range(nocc+nopen):
        print("eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact))
        
    print("      lumo       %20.8f"%(eigen[i+nshift+1]))
    
    erep = bertha.get_erep()
    etotal = bertha.get_etotal()
    ecoul  = bertha.get_eecoul()
    exc    = bertha.get_eexc()
    
    
    print("")
    print("total electronic energy  = %30.15f"%(etotal-(sfact*nocc)))
    print("nuclear repulsion energy = %30.15f"%(erep))
    print("total energy             = %30.15f"%(etotal+erep-(sfact*nocc)))
    print("coulomb energy           = %30.15f"%(ecoul))
    print("Exc     energy           = %30.15f"%(exc))
    
    
    bertha.finalize()
 
    

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

    pberthaopt = pyberthaoption

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

    runspbertha (pberthaopt)
