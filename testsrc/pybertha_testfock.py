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
dumpfiles = 1


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
print("fock[0,10] = ", fockm[0,10])
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

print("")
print("total electronic energy  = %20.8f"%(etotal-(sfact*nocc)))
print("nuclear repulsion energy = %20.8f"%(erep))
print("total energy             = %20.8f"%(etotal+erep-(sfact*nocc)))

bertha.realtime_init()

normalise = 1
direction = 1

vextbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
vextbuffer = numpy.ascontiguousarray(vextbuffer, dtype=numpy.double)

vextm = bertha.get_realtime_dipolematrix(direction, normalise)
#build density
occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)
iocc = 0
for i in range(ndim):
    if i >= nshift and iocc < nocc:
        for j in range(ndim):
            occeigv[j, iocc] = eigem[j, i]
        iocc = iocc + 1
density=numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))
fockm_test = bertha.get_realtime_fock(density.T)
trace=numpy.trace(numpy.matmul(density,ovapm))
print("trace of density dot ovap " , trace)
print("fock_test[0,10] = ", fockm_test[0,10])
diff=fockm-fockm_test
print("max value of diff" ,numpy.max(diff))
diff_trace=numpy.trace(numpy.matmul(density.T,diff))
print('trace diff_dot_dens: ', diff_trace)
bertha.finalize()

