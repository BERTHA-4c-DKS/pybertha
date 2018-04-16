import ctypes
import numpy

import os.path

soname = './bertha_wrapper.so'

if (not os.path.isfile(soname) ):
    print "SO ", soname, " does not exist "
    exit()

bertha = ctypes.cdll.LoadLibrary(soname)

fittcoefffname = "fitcoeff.txt"
vctfilename = "vct.txt" 
ovapfilename = "ovap.txt"
fnameinput = "input.inp"
fittfname = "fitt2.inp"
in_fittcoefffname = ctypes.c_char_p(fittcoefffname)
in_vctfilename = ctypes.c_char_p(vctfilename)
in_ovapfilename = ctypes.c_char_p(ovapfilename)
in_fnameinput = ctypes.c_char_p(fnameinput)
in_fittfname = ctypes.c_char_p(fittfname)

verbosity = 0
bertha.init(in_fnameinput, ctypes.c_int(verbosity))
ndim = bertha.get_ndim()
nshift = bertha.get_nshift()
nocc = bertha.get_nocc()

print "Verbosity: ", verbosity
print "Matrix dimension: ", ndim
print "            nocc: ", nocc
print "          nshift: ", nshift

eigen = numpy.zeros(ndim, dtype=numpy.double)

bertha.mainrun(in_fittcoefffname, in_vctfilename, \
        in_ovapfilename, in_fittfname, \
        ctypes.c_void_p(eigen.ctypes.data))

bertha.get_sfact.restype = ctypes.c_double
bertha.get_etotal.restype = ctypes.c_double

sfact = bertha.get_sfact()
nopen = bertha.get_nopen()

print "           nopen: ", nopen
print "     level shift: ", sfact

for i in range(nocc+nopen):
    print "eigenvalue %5d %12.5f"%(i+1, eigen[i+nshift]-sfact)
    
print "          lumo %12.5f"%(eigen[i+nshift+1])

"""
      write(6,*) 'total electronic energy = ',etotal-(sfact*nocc)
      write(6,*) 'nuclear repulsion energy = ',erep
      write(6,*) 'total energy = ',etotal+erep-(sfact*nocc)
"""

bertha.finalize()
