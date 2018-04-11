import ctypes
import numpy

import os.path

soname = './bertha_wrapper.so'

if (not os.path.isfile(soname) ):
    print "SO ", soname, " does not exist "
    exit()

berthaw = ctypes.cdll.LoadLibrary(soname)

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
berthaw.bertha_init(in_fnameinput, ctypes.c_int(verbosity))
ndim = berthaw.get_ndim()

print "Verbosity: ", verbosity
print "Matrix dimension: ", ndim

a = numpy.zeros(ndim, dtype=numpy.double)
for i in range(ndim):
    a[i] = i

berthaw.bertha_main(in_fittcoefffname, in_vctfilename, \
        in_ovapfilename, in_fittfname, \
        ctypes.c_void_p(a.ctypes.data))

berthaw.bertha_finalize()
