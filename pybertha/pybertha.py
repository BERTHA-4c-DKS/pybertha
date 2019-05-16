import argparse
import ctypes
import numpy
import sys
import re

import os.path

from numpy.linalg import eigvalsh
from scipy.linalg import eigh

sys.path.insert(0, '../src/')
import berthamod

import time

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
parser.add_argument("--tresh", help="det treshold (default = 1.0e-12)", required=False, 
        type=numpy.float64, default=1.0e-12)
parser.add_argument("--wrapperso", help="set wrapper SO (default = ../../lib/bertha_wrapper.so)", 
        required=False, type=str, default="../lib/bertha_wrapper.so")

args = parser.parse_args()

print("Options: ")
print(args) 
print("")
print("")

if not os.path.isfile(args.wrapperso):
    print("SO File ", args.wrapperso, " does not exist")
    exit(1)

bertha = berthamod.pybertha(args.wrapperso)

fittcoefffname = args.fitcoefffile
vctfilename = args.vctfile
ovapfilename = args.ovapfile

fnameinput = args.inputfile
if not os.path.isfile(fnameinput):
    print("File ", fnameinput, " does not exist")
    exit(1)

fittfname = args.fittfile
if not os.path.isfile(fittfname):
    print("File ", fittfname , " does not exist")
    exit(1)

verbosity = args.verbosity
dumpfiles = int(args.dumpfiles)

bertha.set_fittcoefffname(fittcoefffname)
bertha.set_ovapfilename(ovapfilename)
bertha.set_vctfilename(vctfilename)
bertha.set_fnameinput(fnameinput)
bertha.set_fittfname(fittfname)
bertha.set_tresh(args.tresh)

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
cstart = time.clock()

ovapm, eigem, fockm, eigen = bertha.run()

end = time.time()
cend = time.clock()

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

print("")
print("total electronic energy  = %20.8f"%(etotal-(sfact*nocc)))
print("nuclear repulsion energy = %20.8f"%(erep))
print("total energy             = %20.8f"%(etotal+erep-(sfact*nocc)))

occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)

iocc = 0
for i in range(ndim):
    if i >= nshift and iocc < nocc:
        for j in range(ndim):
            occeigv[j, iocc] = eigem[j, i]
        iocc = iocc + 1

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

bertha.realtime_init()
normalise = 1
eps = bertha.get_eps (0.0, 0.0, 0.0)

print (eps)

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
