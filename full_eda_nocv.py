import subprocess
import argparse
import numpy
import shlex
import sys
import os


def execute(cmd):

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    while True:
        output = process.stdout.readline()
        if process.poll() is not None:
            break
        if output:
            print(output.strip().decode("utf-8"), flush=True)

    rc = process.poll()
            
    return rc

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--externalprocess", help="Specify to use an external process for pybertha run",
            required=False, default=False, action="store_true")
    

    parser.add_argument("--fragA", help="Specify the fragA XYZ file", required=True, 
            type=str, default="")
    parser.add_argument("--fragB", help="Specify the fragB XYZ file", required=True, 
            type=str, default="")
    parser.add_argument("--molecule", help="Specify the molecule XYZ file", required=True, 
            type=str, default="")

    parser.add_argument("--energyconverter", help="Specify energy converter (default: 1.0)", required=False, 
            type=float, default=1.0)
    parser.add_argument("--cube", help="Specify if nocv orbital and def. density need to be saved in cube file format (default: False)", required=False,
               default=False, action="store_true")
    parser.add_argument("--deltax", help="cube dx step (default = 0.2)", required=False, 
            type=numpy.float64, default=0.2)
    parser.add_argument("--deltay", help="cube dy step (default = 0.2)", required=False, 
            type=numpy.float64, default=0.2)
    parser.add_argument("--deltaz", help="cube dz step (default = 0.2)", required=False, 
            type=numpy.float64, default=0.2)
    parser.add_argument("--lmargin", help="cube margin parameter (default = 10.0)", required=False, 
            type=numpy.float64, default=10.0)
    parser.add_argument("-np", "--npairs", help="Specify the numerber of nocv-pair density (default: 0)", required=False,
               default=0, type=int)
 
    parser.add_argument("-d", "--debug", help="Debug on, prints debug info to debug_info.txt", required=False, 
            default=False, action="store_true")
    parser.add_argument("-v", "--verbosity", help="Verbosity level 0 = minim, -1 = print iteration info, " + 
            "1 = maximum (defaul -1)", required=False, default=-1, type=int)
    parser.add_argument("--wrapperso", help="set wrapper SO (default = ../../lib/bertha_wrapper.so)", 
            required=False, type=str, default="../lib/bertha_wrapper.so")
    parser.add_argument("--berthamodpaths", help="set berthamod and all other modules path [\"path1;path2;...\"] (default = ../src)", 
            required=False, type=str, default="../src")
          
    parser.add_argument("-j","--jsonbasisfile", \
        help="Specify BERTHA JSON file for fitting and basis (default: fullsets.json)", \
        required=False, type=str, default="fullsets.json")
    parser.add_argument("-b","--basisset", \
        help="Specify BERTHA basisset \"atomname1:basisset1,atomname2:basisset2,...\"", \
        required=True, type=str, default="")
    parser.add_argument("-t","--fittset", \
        help="Specify BERTHA fitting set \"atomname1:fittset1,atomname2:fittset2,...\"", \
        required=True, type=str, default="")
    parser.add_argument("--convertlengthunit", help="Specify a length converter [default=1.0]", \
        type=float, default=1.0)
    parser.add_argument("--functxc", help="EX-POTENTIAL available: LDA,B88P86,HCTH93,BLYP (default=BLYP)", \
        type=str, default="BLYP")

    args = parser.parse_args()

    MAXIT = 100

    for path in args.berthamodpaths.split(";"):
        sys.path.append(path)

    import py_eda_nocv
    import pybertha
    import pybgen

    pyberthaoption_fraga = pybertha.pyberthaoption
    pygenoption_fraga = pybgen.berthainputoption

    pygenoption_fraga.inputfile = args.fragA
    pygenoption_fraga.jsonbasisfile = args.jsonbasisfile
    pygenoption_fraga.fittset = args.fittset
    pygenoption_fraga.basisset = args.basisset
    pygenoption_fraga.functxc = args.functxc
    pygenoption_fraga.convertlengthunit = args.convertlengthunit

    pygenoption_fraga.maxit = MAXIT

    pybgen.generateinputfiles (pygenoption_fraga)

    pyberthaoption_fraga.eda_nocv_info = True
    pyberthaoption_fraga.eda_nocv_frag_file = "info_eda_nocv_fragA.json"

    if args.externalprocess:
       toexe = "python3 ./pybertha.py --eda_nocv_info --eda_nocv_frag_file " \
               + pyberthaoption_fraga.eda_nocv_frag_file
       print("Start pybertha fragA")
       execute(shlex.split(toexe))
       #results  = subprocess.run(toexe, shell=True, check=True, \
       #           stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
       #           universal_newlines=True)
    else:
        pybertha.runspbertha (pyberthaoption_fraga)


    if (os.path.isfile("input.inp"):
        os.remove("input.inp") 
    if (os.path.isfile("fitt2.inp"):
        os.remove("fitt2.inp") 
    if (os.path.isfile("eps.txt"):
        os.remove("eps.txt")
    if (os.path.isfile("fockmtx.txt"):
        os.remove("fockmtx.txt")

    pyberthaoption_fragb = pybertha.pyberthaoption
    pygenoption_fragb = pybgen.berthainputoption

    pygenoption_fragb.inputfile = args.fragB
    pygenoption_fragb.jsonbasisfile = args.jsonbasisfile
    pygenoption_fragb.fittset = args.fittset
    pygenoption_fragb.basisset = args.basisset
    pygenoption_fraga.functxc = args.functxc
    pygenoption_fragb.convertlengthunit = args.convertlengthunit

    pygenoption_fragb.maxit = MAXIT

    pybgen.generateinputfiles (pygenoption_fragb)

    pyberthaoption_fragb.eda_nocv_info = True
    pyberthaoption_fragb.eda_nocv_frag_file = "info_eda_nocv_fragB.json"

    if args.externalprocess:
       toexe = "python3 ./pybertha.py --eda_nocv_info --eda_nocv_frag_file " + \
               pyberthaoption_fragb.eda_nocv_frag_file
       print("Start pybertha fragB")
       execute(shlex.split(toexe))
       #results  = subprocess.run(toexe, shell=True, check=True, \
       #           stdout=subprocess.PIPE, stderr=subprocess.PIPE, \
       #           universal_newlines=True)

    else:
        pybertha.runspbertha (pyberthaoption_fragb)

    if (os.path.isfile("input.inp"):
        os.remove("input.inp") 
    if (os.path.isfile("fitt2.inp"):
        os.remove("fitt2.inp") 
    if (os.path.isfile("eps.txt"):
        os.remove("eps.txt")
    if (os.path.isfile("fockmtx.txt"):
        os.remove("fockmtx.txt")

    pygenoption_mol = pybgen.berthainputoption

    pygenoption_mol.inputfile = args.molecule
    pygenoption_mol.jsonbasisfile = args.jsonbasisfile
    pygenoption_mol.fittset = args.fittset
    pygenoption_mol.basisset = args.basisset
    pygenoption_mol.convertlengthunit = args.convertlengthunit

    pygenoption_mol.maxit = MAXIT

    pybgen.generateinputfiles (pygenoption_mol)
 
    py_eda_nocvoption = py_eda_nocv.ppynocvedaoption

    py_eda_nocvoption.npairs = args.npairs
    py_eda_nocvoption.energyconverter = args.energyconverter

    py_eda_nocvoption.cube = args.cube
    py_eda_nocvoption.deltax = args.deltax
    py_eda_nocvoption.deltay = args.deltay
    py_eda_nocvoption.deltaz = args.deltaz
    py_eda_nocvoption.lmargin = args.lmargin

    py_eda_nocv.runnocveda (py_eda_nocvoption)

    if (os.path.isfile("input.inp"):
        os.remove("input.inp") 
    if (os.path.isfile("fitt2.inp"):
        os.remove("fitt2.inp") 
    if (os.path.isfile("eps.txt"):
        os.remove("eps.txt")
    if (os.path.isfile("fockmtx.txt"):
        os.remove("fockmtx.txt")

    if (os.path.isfile("info_eda_nocv_fragA.json"):
        os.remove("info_eda_nocv_fragA.json")
    if (os.path.isfile("info_eda_nocv_fragB.json"):
        os.remove("info_eda_nocv_fragB.json")

