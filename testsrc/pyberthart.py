import ctypes
import numpy
import sys
import re
import os.path

import rtutil

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
print "fock[0,10] = ", fockm[0,10]
if (fockm is None) or (eigen is None) or (fockm is None) \
        or (eigen is None):
    print "Error in bertha run"
    exit(-1)

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

bertha.realtime_init()

D_0=numpy.zeros((ndim,ndim),dtype=numpy.complex128)
for num in range(nocc):
    D_0[num,num]=1.0+0.0j
dt=0.1
#containers
ene_list = []
dip_list = []
imp_list = []
C = eigem
print type(eigem)
#C_inv used to backtransform D(AO)
try: 
    C_inv = numpy.linalg.inv(eigem)
except LinAlgError:
    print "error" 
test=numpy.matmul(C_inv,eigem)
test1=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(ovapm,C))
#print test1
print numpy.allclose(numpy.eye((ndim),dtype=numpy.complex128),test1)
#build density in ao basis
occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)

iocc = 0

for i in range(ndim):
    if i >= nshift and iocc < nocc:
        for j in range(ndim):
            occeigv[j, iocc] = eigem[j, i]
        iocc = iocc + 1

Da=numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))

direction = 2
normalise = 1

dipz_mat = bertha.get_realtime_dipolematrix(direction, normalise)


fock_mid_init = rtutil.mo_fock_mid_forwd_eval(bertha,Da,fockm,0,numpy.float_(dt),\
	dipz_mat,C,C_inv,ovapm,ndim)


bertha.finalize()

