import ctypes
import numpy
import sys
import re

import os.path

from numpy.linalg import eigvalsh
from scipy.linalg import eigh

sys.path.insert(0, '../src/')
import berthamod

bertha = berthamod.pybertha("../../lib/bertha_wrapper.so")

fittcoefffname = "fitcoeff.txt"
vctfilename = "vct.txt" 
ovapfilename = "ovap.txt"
fnameinput = "input.inp"
fittfname = "fitt2.inp"

verbosity = -1
dumpfiles = 0


bertha.set_fittcoefffname(fittcoefffname)
bertha.set_ovapfilename(ovapfilename)
bertha.set_vctfilename(vctfilename)
bertha.set_fnameinput(fnameinput)
bertha.set_fittfname(fittfname)

bertha.set_verbosity(verbosity)
bertha.set_dumpfiles(dumpfiles)

bertha.init()

ndim = bertha.get_ndim()
nshift = bertha.get_nshift()
nocc = bertha.get_nocc()
sfact = bertha.get_sfact()
nopen = bertha.get_nopen()

print "Verbosity       : ", verbosity
print "Dumpfiles       : ", dumpfiles
print ""
print "Matrix dimension: ", ndim
print "            nocc: ", nocc
print "          nshift: ", nshift
print "           nopen: ", nopen
print "     level shift: ", sfact
print ""

ovapm, eigem, fockm, eigen = bertha.run()

if (fockm is None) or (eigen is None) or (fockm is None) \
        or (eigen is None):
    print "Error in bertha run"
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

print ""
print "Final results "
for i in range(nocc+nopen):
    print "eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact)
    
print "      lumo       %20.8f"%(eigen[i+nshift+1])

erep = bertha.get_erep()
etotal = bertha.get_etotal()

print ""
print "total electronic energy  = %20.8f"%(etotal-(sfact*nocc))
print "nuclear repulsion energy = %20.8f"%(erep)
print "total energy             = %20.8f"%(etotal+erep-(sfact*nocc))

bertha.finalize()

occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex64)

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

print ""
print "Compute density matrix "
density = numpy.matmul(occeigv, numpy.conjugate(occeigv.transpose()), out=None)
density = numpy.matmul(density.transpose(), ovapm)
print "Trace  "
trace = density.trace()
print "(%20.10f, %20.10fi)"%(trace.real, trace.imag)

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
