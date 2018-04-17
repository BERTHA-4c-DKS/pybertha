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

verbosity = -1
dumpfiles = 0

bertha.init(in_fnameinput, ctypes.c_int(verbosity), 
        ctypes.c_int(dumpfiles))
ndim = bertha.get_ndim()
nshift = bertha.get_nshift()
nocc = bertha.get_nocc()
bertha.get_sfact.restype = ctypes.c_double
bertha.get_etotal.restype = ctypes.c_double
bertha.get_erep.restype = ctypes.c_double
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

eigen = numpy.zeros(ndim, dtype=numpy.double)

bertha.mainrun(in_fittcoefffname, in_vctfilename, \
        in_ovapfilename, in_fittfname, \
        ctypes.c_void_p(eigen.ctypes.data))

print ""
print "Final results "
for i in range(nocc+nopen):
    print "eigenvalue %5d %12.5f"%(i+1, eigen[i+nshift]-sfact)
    
print "      lumo       %12.5f"%(eigen[i+nshift+1])

erep = bertha.get_erep()
etotal = bertha.get_etotal()

print ""
print "total electronic energy  = %12.5f"%(etotal-(sfact*nocc))
print "nuclear repulsion energy = %12.5f"%(erep)
print "total energy             = %12.5f"%(etotal+erep-(sfact*nocc))

bertha.finalize()
