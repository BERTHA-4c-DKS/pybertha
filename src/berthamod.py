import threading
import ctypes
import numpy
import sys
import re

import os.path
import time

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
            cbuffer[  counter] = inm[j, i].real
            cbuffer[counter+1] = inm[j, i].imag
            counter = counter + 2

    return cbuffer

###############################################################################

def doublevct_to_complexmat (invector, dim):

    if (invector.size != (2*dim*dim)):
        return None

    outm = numpy.zeros((dim,dim), dtype=numpy.complex128)
    counter = 0
    for j in range(dim):
        for i in range(dim):
            outm[j, i] = complex(invector[counter], invector[counter+1])
            counter = counter + 2

    return outm

###############################################################################

def read_ovapfile (fname):
    
    fp = open(fname)
    ndim = int(fp.readline())

    ovapcmp = numpy.zeros((ndim,ndim), dtype=numpy.complex128)
    
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

    vct = numpy.zeros((ndim,nocc), dtype=numpy.complex128)
    
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

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class pybertha:
    
    def __init__(self, sopath="./bertha_wrapper.so"):
        soname = sopath
        if (not os.path.isfile(soname) ):
            raise Error("SO %s  does not exist" % soname)

        self.__bertha = ctypes.cdll.LoadLibrary(soname)
        
        self.__reset()

    def __reset(self):

        self.__mainruntime = 0.0
        self.__mainrunctime = 0.0

        self.__fittcoefffname = "fitcoeff.txt"
        self.__vctfilename = "vct.txt" 
        self.__ovapfilename = "ovap.txt"
        self.__fnameinput = "input.inp"
        self.__fittfname = "fitt2.inp"
        self.__verbosity = -1
        self.__dumpfiles = 0

        self.__init = False
        self.__mainrundone = False

        self.__realtime_init = False

        self.set_densitydiff (0)

    def get_mainruntime(self):
        return self.__mainruntime

    def get_mainrunctime(self):
        return self.__mainrunctime

    def set_densitydiff (self, ini):
        self.__bertha.set_densitydiff(ctypes.c_int(ini)) 

    def get_densitydiff (self):
        self.__bertha.get_densitydiff.restype = ctypes.c_int
 
        return self.__bertha.get_densitydiff()

    def set_tresh (self, ini):
        self.__bertha.set_tresh(ctypes.c_double(ini)) 

    def get_tresh (self):
        self.__bertha.get_tresh.restype = ctypes.c_double
 
        return self.__bertha.get_tresh()

    def set_fittcoefffname (self, ini):
        self.__fittcoefffname = ini

    def get_fittcoefffname (self):
        return self.__fittcoefffname

    def set_vctfilename(self, ini):
        self.__vctfilename = ini

    def get_vctfilename(self):
        return self.__vctfilename

    def set_ovapfilename(self, ini):
        self.__ovapfilename = ini

    def get_ovapfilename(self):
        return self.__ovapfilename

    def set_fnameinput (self, ini):
        self.__fnameinput = ini

    def get_fnameinput (self):
        return self.__fnameinput

    def set_fittfname (self, ini):
        self.__fittfname = ini

    def get_fittfname (self):
        return self.__fittfname

    def set_verbosity (self, ini):
        self.__verbosity = ini

    def get_verbosity (self):
        return self.__verbosity

    def set_dumpfiles (self, ini):
        self.__dumpfiles = ini

    def get_dumpfiles (self):
        return self.__dumpfiles

    def init (self):
       
        in_fnameinput = ctypes.c_char_p(self.__fnameinput)

        self.__bertha.init(in_fnameinput, ctypes.c_int(self.__verbosity), 
                ctypes.c_int(self.__dumpfiles))

        self.__init = True

    def get_ndim(self):
        if self.__init:
            return self.__bertha.get_ndim()
        
        return 0

    def get_nshift(self):
        if self.__init:
            return self.__bertha.get_nshift()
        
        return 0

    def get_nocc(self):
        if self.__init:
            return self.__bertha.get_nocc()
        
        return 0

    def get_nopen(self):
        if self.__init:
            return self.__bertha.get_nopen()
        
        return 0

    def get_sfact(self):
        self.__bertha.get_sfact.restype = ctypes.c_double

        if self.__init:
            return self.__bertha.get_sfact()
        
        return 0

    def density_to_cube(self, dens, fname, margin = 10.0, 
            drx = 0.2, dry = 0.2, drz = 0.2):

        in_fittfname = ctypes.c_char_p(self.__fittfname)
        in_fname = ctypes.c_char_p(fname)

        cbuffer = complexmat_to_doublevct (dens)

        self.__bertha.density_to_cube (ctypes.c_void_p(cbuffer.ctypes.data), \
                ctypes.c_double(margin), ctypes.c_double(drx), ctypes.c_double(dry), \
                ctypes.c_double(drz), in_fname, in_fittfname)

    def run(self):
        if self.__init:
            ndim = self.get_ndim()

            eigen = numpy.zeros(ndim, dtype=numpy.double)
            eigen = numpy.ascontiguousarray(eigen, dtype=numpy.double)

            ovapbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
            ovapbuffer = numpy.ascontiguousarray(ovapbuffer, dtype=numpy.double)
            
            fockbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
            fockbuffer = numpy.ascontiguousarray(fockbuffer, dtype=numpy.double)
            
            eigenvctbu = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
            eigenvctbu = numpy.ascontiguousarray(eigenvctbu, dtype=numpy.double)

            in_fittcoefffname = ctypes.c_char_p(self.__fittcoefffname)
            in_vctfilename = ctypes.c_char_p(self.__vctfilename)
            in_ovapfilename = ctypes.c_char_p(self.__ovapfilename)
            in_fittfname = ctypes.c_char_p(self.__fittfname)

            start = time.time()
            cstart = time.clock()

            maint = threading.Thread(target=self.__bertha.mainrun, \
                    args=[in_fittcoefffname, \
                          in_vctfilename, \
                          in_ovapfilename, \
                          in_fittfname, \
                          ctypes.c_void_p(eigen.ctypes.data), \
                          ctypes.c_void_p(ovapbuffer.ctypes.data), \
                          ctypes.c_void_p(eigenvctbu.ctypes.data), \
                          ctypes.c_void_p(fockbuffer.ctypes.data)])
            maint.daemon = True
            maint.start()
            while maint.is_alive():
                    maint.join(.1)

            end = time.time()
            cend = time.clock()

            self.__mainruntime = end - start
            self.__mainrunctime = cend - cstart
 
            #self.__bertha.mainrun(in_fittcoefffname, in_vctfilename, \
            #        in_ovapfilename, in_fittfname, \
            #        ctypes.c_void_p(eigen.ctypes.data), \
            #        ctypes.c_void_p(ovapbuffer.ctypes.data), \
            #        ctypes.c_void_p(eigenvctbu.ctypes.data), \
            #        ctypes.c_void_p(fockbuffer.ctypes.data))
            
            ovapm = doublevct_to_complexmat (ovapbuffer, ndim)
            if ovapm is None:
                raise Error("Error in ovap matrix size")

            eigem = doublevct_to_complexmat (eigenvctbu, ndim)
            if eigem is None:
                raise Error("Error in ovap matrix size")
                
            fockm = doublevct_to_complexmat (fockbuffer, ndim)
            if fockm is None:
                raise Error("Error in ovap matrix size")

            self.__mainrundone = True

            return ovapm, eigem, fockm, eigen

        return None, None, None, None

    def get_etotal(self):
        self.__bertha.get_etotal.restype = ctypes.c_double

        if self.__mainrundone:
            return self.__bertha.get_etotal()

        return 0.0

    def get_erep(self):
        self.__bertha.get_erep.restype = ctypes.c_double

        if self.__mainrundone:
            return self.__bertha.get_erep()

        return 0.0

    def realtime_init(self):
        if self.__init and self.__mainrundone:
            self.__bertha.realtime_init()
            self.__realtime_init = True
        else:
            raise Error("You should firstly init and perform the run")

    def get_realtime_dipolematrix (self, direction, normalise):
        if self.__realtime_init:

            ndim = self.get_ndim()

            vextbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
            vextbuffer = numpy.ascontiguousarray(vextbuffer, dtype=numpy.double)
            
            self.__bertha.realtime_dipolematrix(direction, normalise, \
                    ctypes.c_void_p(vextbuffer.ctypes.data))
            
            vextm = doublevct_to_complexmat (vextbuffer, ndim)
            if vextm is None:
                raise Error("Error in vextm matrix size")
            
            return vextm

        else:
            return None

    def get_realtime_fock (self, eigem):
        if self.__realtime_init:

            ndim = self.get_ndim()

            cbuffer = complexmat_to_doublevct (eigem)

            fockbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
            fockbuffer = numpy.ascontiguousarray(fockbuffer, dtype=numpy.double)

            self.__bertha.realtime_fock(ctypes.c_void_p(cbuffer.ctypes.data), \
                    ctypes.c_void_p(fockbuffer.ctypes.data))

            fockm = doublevct_to_complexmat (fockbuffer, ndim)
            if fockm is None:
                raise Error("Error in ovap matrix size")

            return fockm

        return None

    def finalize(self):
        if self.__realtime_init:
            self.__bertha.realtime_finalize()
        
        if self.__init:
            self.__bertha.finalize()

        self.__reset()
