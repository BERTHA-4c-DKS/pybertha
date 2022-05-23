import argparse
import numpy
import uuid
import sys

import os.path

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
    cubebox: bool = False
    cube: bool = False
    xmin: numpy.float64 = -10.0
    ymin: numpy.float64 = -10.0
    zmin: numpy.float64 = -10.0
    xmax: numpy.float64 = 10.0
    ymax: numpy.float64 = 10.0
    zmax: numpy.float64 = 10.0
    deltax: numpy.float64 = 0.2
    deltay: numpy.float64 = 0.2
    deltaz: numpy.float64 = 0.2
    lmargin: numpy.float64 = 10.0
    gridfilename: str = ""
    potfilename: str =  ""

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
    
    #import os, psutil

    #process = psutil.Process(os.getpid())
    #print(process.memory_info())  # in bytes 

    
    print("Options: ")
    for att in [a for a in dir(pberthaopt) if not a.startswith('__')]:
        print(att, " = ", getattr(pberthaopt, att)) 
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
    
    bertha.set_densitydiff(0)
    
    bertha.init()

    #print(process.memory_info())  # in bytes 

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

    #print(process.memory_info())  # in bytes

    if pberthaopt.gridfilename != "" and \
        pberthaopt.potfilename != "":

        if os.path.isfile(pberthaopt.gridfilename) and \
            os.path.isfile(pberthaopt.potfilename):
            grid = berthamod.read_sgrid_file (pberthaopt.gridfilename)
            pot = berthamod.read_pot_file (pberthaopt.potfilename)

            if (pot is not None) and (grid is not None):
                if grid.shape[0] == pot.shape[0]:
                    bertha.set_embpot_on_grid(grid, pot)
                else:
                    print("ERROR: grid and pot files are not compatible")
                    exit(1)
        else:
            print("ERROR: "+ pberthaopt.gridfilename + " and/or " + 
               pberthaopt.potfilename + " do not exist")
            exit(1)


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
    
    
    """
    counter = 0
    for i in range(ndim):
          print "i ==> ", i+1, eigen[i]
          for j in range(ndim):
              sys.stdout.write("(%20.10f %20.10fi) \n"%(
                    eigenvctbu[counter], eigenvctbu[counter+1]))
              counter = counter + 2
          print ""
    """
    
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
    
    occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)
    
    etotal = etotal+erep-(sfact*nocc)
    
    
    iocc = 0
    for i in range(ndim):
        if i >= nshift and iocc < nocc:
            for j in range(ndim):
                occeigv[j, iocc] = eigem[j, i]
            iocc = iocc + 1
    
    if pberthaopt.eda_nocv_info :
       encoder.FLOAT_REPR = lambda o: format(o, '.25E')
       json_data = get_json_data(pberthaopt, etotal, erep, ecoul, exc, ndim, nocc, occeigv)
       with open(pberthaopt.eda_nocv_frag_file, 'w') as fp:
           json.dump(json_data, fp, sort_keys=True, indent=4)
    
    """
    for i in range(nocc):
          print "i ==> ", i+1, eigen[i+nshift]
          for j in range(ndim):
              sys.stdout.write("(%20.10f %20.10fi)\n"%(
                  occeigv[j, i].real, occeigv[j, i].imag))
          print ""
    """
    
    print("")
    print("Compute density matrix ")
    density = numpy.matmul(occeigv, numpy.conjugate(occeigv.transpose()), out=None)
    density = numpy.matmul(density, ovapm)
    print("Trace  ")
    trace = density.trace()
    print("(%20.10f, %20.10fi)"%(trace.real, trace.imag))

    drx = pberthaopt.deltax
    dry = pberthaopt.deltay
    drz = pberthaopt.deltaz
    margin = pberthaopt.lmargin


    if (pberthaopt.cube == True):
       bertha.density_to_cube((density).T, "density.cube", margin, drx, dry, drz )

    if (pberthaopt.cubebox == True):
       bertha.density_to_cube_limit((density).T, "densitybox.cube", \
               (pberthaopt.xmin, pberthaopt.ymin, \
               pberthaopt.zmin), (pberthaopt.xmax, \
               pberthaopt.ymax, pberthaopt.zmax), \
               drx, dry, drz )
    
    bertha.realtime_init()
    normalise = 1
    
   
    numofatom = bertha.get_natoms()
    print ("Num of atoms: ",numofatom)
    for i in range(numofatom):
        print(bertha.get_coords(i))
    

    """
    outep = open("eps.txt", "w")
    
    griddim = 50
    
    xmin = 9.0
    xmax = 12.0
    
    ymin = 6.0
    ymax = 10.0
    
    zmax =  7.0
    zmin = -3.0
 
    dx = (xmax - xmin)/float(griddim)
    dy = (ymax - ymin)/float(griddim)
    x = xmin
    for i in range(griddim):
        y = ymin
        for j in range(griddim):
            z = 0.0
    
            eps = bertha.get_eps (x, y, z)
    
            outep.write (" %10.5e %10.5e %10.5e \n"%(x, y, eps))
    
            y = y + dy
        x = x + dx
    
    dz = (zmax - zmin)/float(griddim)
    z = zmin
    for i in range(griddim):
    
        eps = bertha.get_eps (0.0, 0.0, z)
    
        outep.write (" %10.5e %10.5e \n"%(z, eps))
    
        z = z + dz

    """
    
    bertha.finalize()
    
    """
    # to check if needed
    eigvals, eigvecs = eigh(fockm, ovapm, eigvals_only=False)
    
    iocc = 0 
    for i in range(ndim): 
        if i >= nshift and iocc < nocc:
            print eigvals[i] - sfact
            iocc = iocc + 1
    
    
    ovapcmp = berthamod.read_ovapfile ("ovap.txt")
    
    for i in range(ndim):
        for j in range(ndim):
            if (ovapcmp[i, j] != ovapm[i, j]):
              sys.stdout.write("(%20.10f %20.10fi) -> (%20.10f %20.10fi) \n"%(ovapm[i, j].real, ovapm[i, j].imag,
                  ovapcmp[i, j].real, ovapcmp[i, j].imag))
        #sys.stdout.write("\n")
    
    eigecmp = berthamod.read_vctfile ("vct.txt")
    
    for i in range(ndim):
          print "i ==> ", i+1, eigen[i]
          for j in range(ndim):
            if (eigecmp[i, j] != eigem[i, j]):
              sys.stdout.write("(%20.10f %20.10fi) -> (%20.10f %20.10fi) \n"%(
                  eigem[i, j].real, eigem[i, j].imag,
                  eigecmp[i, j].real, eigecmp[i, j].imag))
          print ""
    """


##########################################################################################

if __name__ == "__main__":

    MAXIT = 100

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--inputfile", help="Specify BERTHA input file (default: input.inp)", required=False, 
            type=str, default="input.inp")
    parser.add_argument("-t","--fittfile", help="Specify BERTHA fitting input file (default: fitt2.inp)", required=False, 
            type=str, default="fitt2.inp")

    parser.add_argument("-g","--geometry", help="Specify system (Angstrom) geometry", required=False, 
        type=str, default="")
    parser.add_argument("--obs", \
        help="Specify BERTHA (Active system) basisset \"atomname1:basisset1,atomname2:basisset2,...\"", \
        required=False, type=str, default="")
    parser.add_argument("--fittset", \
        help="Specify BERTHA (Active system) fitting set \"atomname1:fittset1,atomname2:fittset2,...\"", \
        required=False, type=str, default="")
    parser.add_argument("--func", 
	help="Specify exchange–correlation energy functional for active system available: LDA,B88P86,HCTH93,BLYP (default=BLYP)", \
        type=str, default="BLYP")
    parser.add_argument("-j","--jsonbasisfile", \
        help="Specify BERTHA JSON file for fitting and basis (default: fullsets.json)", \
       required=False, type=str, default="fullsets.json")
    parser.add_argument("--convertlengthunit", help="Specify a length converter [default=1.0]", \
        type=float, default=1.0)
    parser.add_argument("--berthamaxit", help="set bertha maxiterations (default = %d)"%(MAXIT), 
        required=False, type=int, default=MAXIT)

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
    parser.add_argument("--modpaths", help="set berthamod path (default = ../src)", 
            required=False, type=str, default="../src")
    parser.add_argument("--eda_nocv_info", help="set to dump info useful for py_eda_nocv",action='store_true',default=False)
    parser.add_argument("--eda_nocv_frag_file", help="set a file (default: info_eda_nocv_fragX.json)",
            required=False, type=str, default="info_eda_nocv_fragX.json")
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

    parser.add_argument("--cubebox", help="Specify if need to dump density on a cubeBOX(default: 0)", required=False,
               default=False, action="store_true")
    parser.add_argument("--xmin", help="cube xmin (default = -10.0)", required=False,
            type=numpy.float64, default=-10.0)
    parser.add_argument("--ymin", help="cube ymin (default = -10.0)", required=False,
            type=numpy.float64, default=-10.0)
    parser.add_argument("--zmin", help="cube zmin (default = -10.0)", required=False,
            type=numpy.float64, default=-10.0)
    parser.add_argument("--xmax", help="cube xmax (default =  10.0)", required=False,
            type=numpy.float64, default=10.0)
    parser.add_argument("--ymax", help="cube ymax (default =  10.0)", required=False,
            type=numpy.float64, default=10.0)
    parser.add_argument("--zmax", help="cube zmax (default =  10.0)", required=False,
            type=numpy.float64, default=10.0)

    parser.add_argument("--addgridpot", help="Import a custom grid potential [\"gridfile.txt;pofile.txt\"]", required=False,
            type=str, default="")

    args = parser.parse_args()

    for path in args.modpaths.split(";"):
        sys.path.append(path)

    modpaths = os.environ.get('PYBERTHA_MOD_PATH')

    if modpaths is not None :
        for path in modpaths.split(";"):
            sys.path.append(path)

    gridfilename = ""
    potfilename = ""

    if args.addgridpot != "":
        if len(args.addgridpot.split(";")) == 2:
            gridfilename = args.addgridpot.split(";")[0]
            potfilename = args.addgridpot.split(";")[1]
        else:
            print("ERROR in --addgridpot, you need to specify two filenames ; separated")
            exit(1)

    pberthaopt = pyberthaoption

    if args.geometry != "":

        if not os.path.isfile (args.geometry):
            print("File ", args.geometry, " does not exist")
            exit(1)

        if args.fittset == "" or args.obs == "":
            print("Need to specify fittset and basis set for each elements")
            exit(1)

        # generate input 
        import pybgen

        pygenoption = pybgen.berthainputoption
        pygenoption.inputfile = args.geometry
        pygenoption.jsonbasisfile = args.jsonbasisfile
        pygenoption.fittset = args.fittset
        pygenoption.basisset = args.obs
        pygenoption.functxc = args.func
        pygenoption.convertlengthunit = args.convertlengthunit
        pygenoption.maxit = args.berthamaxit

        #pberthaopt.inputfile = "input.inp"
        #pberthaopt.fittfile = "fitt2.inp"
       
        pberthaopt.inputfile = str(uuid.uuid4())
        pberthaopt.fittfile = str(uuid.uuid4())  

        pygenoption.berthainfname = pberthaopt.inputfile
        pygenoption.berthafittfname = pberthaopt.fittfile 

        for filename in [pberthaopt.inputfile , pberthaopt.fittfile]:

            if os.path.isfile(filename):
                print("File ", filename, " will be overwritten")
                try:
                    os.remove(filename)
                except OSError:
                    pass

        pybgen.generateinputfiles (pygenoption)

    else:
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
    pberthaopt.berthamodpath = args.modpaths
    pberthaopt.eda_nocv_info = args.eda_nocv_info
    pberthaopt.eda_nocv_frag_file = args.eda_nocv_frag_file
    pberthaopt.cube = args.cube
    pberthaopt.deltax = args.deltax
    pberthaopt.deltay = args.deltay
    pberthaopt.deltaz = args.deltaz
    pberthaopt.lmargin = args.lmargin
    pberthaopt.cubebox = args.cubebox
    pberthaopt.xmin = args.xmin
    pberthaopt.ymin = args.ymin
    pberthaopt.zmin = args.zmin
    pberthaopt.xmax = args.xmax
    pberthaopt.ymax = args.ymax
    pberthaopt.zmax = args.zmax
    pberthaopt.gridfilename = gridfilename
    pberthaopt.potfilename = potfilename

    #import resource
    #rsrc = resource.RLIMIT_STACK
    #soft, hard = resource.getrlimit(rsrc)
    #print("Limit starts as:", soft, hard)
    #val = 100000000
    #resource.setrlimit(rsrc, (val, val))
    #soft, hard = resource.getrlimit(rsrc)
    #print("Limit is now:", soft, hard)
    #sys.setrecursionlimit(10**6)

    runspbertha (pberthaopt)

    if args.geometry != "":
        for filename in [pberthaopt.inputfile , pberthaopt.fittfile]:
            if os.path.isfile(filename):
                try:
                    os.remove(filename)
                except OSError:
                    pass