import ctypes
import numpy
import sys
import re

import os.path

from numpy.linalg import eigvalsh
from scipy.linalg import eigh

###############################################################################

def complexmat_to_doublevct (inm):

    if len(inm.shape) != 2:
        return None

    if inm.shape[0] != inm.shape[1]:
        return None

    dim = inm.shape[0]
    
    cbuffer = numpy.zeros((2*dim*dim), dtype=numpy.double)
    cbuffer = numpy.ascontiguousarray(cbuffer, dtype=numpy.double)

    counter = 0
    for j in range(dim):
        for i in range(dim):
            cbuffer[  counter] = inm[i, j].real
            cbuffer[counter+1] = inm[i, j].imag
            counter = counter + 2

    return cbuffer

###############################################################################

def doublevct_to_complexmat (invector, dim):

    if (invector.size != (2*dim*dim)):
        return None

    outm = numpy.zeros((dim,dim), dtype=numpy.complex64)
    counter = 0
    for j in range(dim):
        for i in range(dim):
            outm[i, j] = complex(invector[counter], invector[counter+1])
            counter = counter + 2

    return outm

###############################################################################

def read_ovapfile (fname):
    
    fp = open(fname)
    ndim = int(fp.readline())

    ovapcmp = numpy.zeros((ndim,ndim), dtype=numpy.complex64)
    
    for j in range(ndim):
        for i in range(ndim):
            line = fp.readline()
            
            line = re.sub( '\s+', ' ', line).strip()
            sline = line.split(" ")
            ovapcmp[i, j] = complex(numpy.float64(sline[0]), numpy.float64(sline[1]))
    
    fp.close()

    return ovapcmp

###############################################################################

def read_vctfile (fname):
    
    fp = open(fname)
    line = fp.readline()
    line = re.sub( '\s+', ' ', line).strip()
    sline = line.split(" ")
 
    nocc = int(sline[0])
    ndim = int(sline[1])

    vct = numpy.zeros((ndim,nocc), dtype=numpy.complex64)
    
    for i in range(nocc):
        line = fp.readline()
        for j in range(ndim):
            line = fp.readline()
            
            line = re.sub( '\s+', ' ', line).strip()
            sline = line.split(" ")
            vct[j, i] = complex(numpy.float64(sline[0]), numpy.float64(sline[1]))
    
    fp.close()

    return vct

###############################################################################

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

ovapbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
ovapbuffer = numpy.ascontiguousarray(ovapbuffer, dtype=numpy.double)

fockbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
fockbuffer = numpy.ascontiguousarray(fockbuffer, dtype=numpy.double)

eigenvctbu = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
eigenvctbu = numpy.ascontiguousarray(eigenvctbu, dtype=numpy.double)

bertha.mainrun(in_fittcoefffname, in_vctfilename, \
        in_ovapfilename, in_fittfname, \
        ctypes.c_void_p(eigen.ctypes.data), \
        ctypes.c_void_p(ovapbuffer.ctypes.data), 
        ctypes.c_void_p(eigenvctbu.ctypes.data),
        ctypes.c_void_p(fockbuffer.ctypes.data))

ovapm = doublevct_to_complexmat (ovapbuffer, ndim)
if ovapm is None:
    print "Error in ovap matrix size"
    exit(-1)

eigem = doublevct_to_complexmat (eigenvctbu, ndim)
if eigem is None:
    print "Error in ovap matrix size"
    exit(-1)

fockm = doublevct_to_complexmat (fockbuffer, ndim)
if fockm is None:
    print "Error in ovap matrix size"
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

normalise = 1
direction = 1

vextbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
vextbuffer = numpy.ascontiguousarray(vextbuffer, dtype=numpy.double)

bertha.realtime_dipolematrix(direction, normalise, \
        ctypes.c_void_p(vextbuffer.ctypes.data))

vextm = doublevct_to_complexmat (vextbuffer, ndim)
if vextm is None:
    print "Error in ovap matrix size"
    exit(-1)

cbuffer = complexmat_to_doublevct (eigem)

bertha.realtime_fock(ctypes.c_void_p(cbuffer.ctypes.data), \
        ctypes.c_void_p(fockbuffer.ctypes.data))

bertha.realtime_finalize()

bertha.finalize()
