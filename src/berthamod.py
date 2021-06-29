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

    cbuffer[0::2] = inm.flatten().real
    cbuffer[1::2] = inm.flatten().imag

    return cbuffer

###############################################################################

def doublevct_to_complexmat (invector, dim):

    if (invector.size != (2*dim*dim)):
        return None

    outm = numpy.zeros((dim,dim), dtype=numpy.complex128)

    inmtxreal = numpy.reshape(invector[0::2], (dim,dim))
    inmtximag = numpy.reshape(invector[1::2], (dim,dim))
    outm[:,:] = inmtxreal[:,:] + 1j * inmtximag[:,:]

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
        """
        param: sopath is needed to specify the 
        bertha_wrapper Shared Object file.
        """
        soname = sopath
        if (not os.path.isfile(soname) ):
            raise Error("SO %s  does not exist" % soname)

        self.__bertha = ctypes.cdll.LoadLibrary(soname)
        
        self.__reset()

    def __reset(self):

        self.__mainruntime = 0.0
        self.__mainrunctime = 0.0

        self.__focktime = 0.0
        self.__fockctime = 0.0

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

        self.__childthread = False

        self.set_densitydiff (0)

    def get_density_on_grid (self, grid):

        # in true o mainrun o fockmtx 
        # grid N x 4 array containing x,y,z,w double  
        # chiamata interno su bertha_wrapper che viene poi implemntata su bertha_ng

        density = None # vector N double

        if self.__realtime_init or self.__mainrundone:
            if isinstance(grid,(numpy.ndarray)):
                if not grid.flags['C_CONTIGUOUS']:
                    raise TypeError("get_density_on_grid: " \
                        + "input must be C_CONTIGUOUS see numpy.ascontiguousarray( ")

                if len(grid.shape) == 2:
                    if grid.shape[1] == 4:
                        npoints = grid.shape[0]

                        density = numpy.zeros(npoints, dtype=numpy.double)
                        density = numpy.ascontiguousarray(density, dtype=numpy.double)

                        # call to main function to get the density on the grid
                        # call in self.__bertha
                        self.__bertha.bertha_get_density_ongrid(ctypes.c_int(npoints), \
                            ctypes.c_void_p(grid.ctypes.data), \
                            ctypes.c_void_p(density.ctypes.data))

                    else:
                        raise TypeError("get_density_on_grid: input must be a 2D numpy.ndarray with 4 columns")
                else:
                    raise TypeError("get_density_on_grid: input must be a 2D numpy.ndarray")
            else:
                raise TypeError("get_density_on_grid: input must be a numpy.ndarray")

        return density

    def set_embpot_on_grid (self, grid, pot):

        # solo init true rivedere

        # pot vector N double

        # grid N x 4 array containing x,y,z,w double  

        # in bertha wrapper rimane e serve sia per mainrun che per get_fovk_realtime

        # questo dato va in bertha_wrapper e poi con if decide se sommare o meno duranti il ciclo SCF
        # e durante la chiamata a get_fock_realtime

        if self.__init:
            if not isinstance(grid,(numpy.ndarray)):
                raise TypeError("set_embpot_on_grid: input must be a numpy.ndarray")

            if not grid.flags['C_CONTIGUOUS']:
                raise TypeError("get_density_on_grid: " \
                        + "input grid must be C_CONTIGUOUS see numpy.ascontiguousarray( ")

            if not isinstance(pot,(numpy.ndarray)):
                raise TypeError("set_embpot_on_grid: input must be a numpy.ndarray")

            if not pot.flags['C_CONTIGUOUS']:
                raise TypeError("get_density_on_grid: " \
                        + "input pot must be C_CONTIGUOUS see numpy.ascontiguousarray( ")

            if len(grid.shape) == 2 and len(pot.shape) == 1:
                if grid.shape[1] == 4:
                    npoints = grid.shape[0]

                    if pot.shape[0] != npoints:
                        raise TypeError("set_embpot_on_grid: incompatible dimensions ")

                    # call to main function to set the embed potentil and grid
                    # and a flag to be called in scf or get reltime_fock
                    # call in self.__bertha

                else:
                    raise TypeError("set_embpot_on_grid: input must be a 2D numpy.ndarray with 4 columns")
            else:
                raise TypeError("set_embpot_on_grid: grid must be a 2D numpy.ndarray and pot a 1D ndarray")

        return

    def get_childthread(self):
        """
        Return childthread  value
        """

        return self.__childthread

    def set_childthread(self, val):

        """
        Set childthread, if True the run spawn a thread to call 
        the BERTHA mainrun. If False (to be used with Intel compiler 
        and OMP because of stackseze realted issue) a thread is not 
        spawned, in cusch a case cannot stop the main run 
        using a Ctrl-C.
        """

        if not isinstance(val, bool):
            raise TypeError("set_densitydiff: input must be a bool")

        self.__childthread = val

    def get_natoms(self):
        """
        Return natoms 
        """

        self.__bertha.get_ncent.restype = ctypes.c_int

        return self.__bertha.get_ncent()

    def get_mainruntime(self):
        """
        Returns the wall time to perfom all the SCF iterations, 
        not including he time needed to copy the arrays.
        """

        return self.__mainruntime

    def get_mainrunctime(self):
        """
        Returns the CPU time to perfom all the SCF iterations, 
        not including the time needed to copy the arrays.
        """

        return self.__mainrunctime

    def get_focktime(self):
        """
        Returns the wall time to perfom get_realtime_fock,
        not including the time needed to copy the arrays.
        """
 
        return self.__focktime

    def get_fockctime(self):
        """
        Returns the CPU time to perfom get_realtime_fock,
        not including the time needed to copy the arrays.
        """
 
        return self.__fockctime

    def set_densitydiff (self, ini):
        """
        Set densitydiff flag value.
        - If set to 1 the code will compute the maximum 
          difference (elementwise) of the density matrix 
          to the iteration n with respect to the one 
          related to the previous iteration (n-1).
        """

        if not isinstance(ini, int):
            raise TypeError("set_densitydiff: input must be an integer")
 
        self.__bertha.set_densitydiff(ctypes.c_int(ini)) 

    def get_densitydiff (self):
        """
        The *get_densitydiff* returns the **densitydiff flag** value.
        """
 
        self.__bertha.get_densitydiff.restype = ctypes.c_int
 
        return self.__bertha.get_densitydiff()

    def set_thresh (self, ini):
        """
        Set the SCF convergence threshold. 
        """
        if not isinstance(ini, float):
            raise TypeError("set_thresh: input must be a float")

        self.__bertha.set_thresh(ctypes.c_double(ini)) 

    def get_thresh (self):
        """
        The *get_thresh* returns the SCF convergence threshold.
        """

        self.__bertha.get_thresh.restype = ctypes.c_double
 
        return self.__bertha.get_thresh()

    def set_fittcoefffname (self, ini):
        """
        To specify the filename where the density fitting 
        coefficients will be written if the dumpfiles flag
        is equal to 1.
        """

        if not isinstance(ini, str):
            raise TypeError("set_fittcoefffname: input must be a string")


        self.__fittcoefffname = ini

    def get_fittcoefffname (self):
        """
        get_fittcoefffname returns the density fitting 
        coefficients filename.
        """

        return self.__fittcoefffname

    def set_vctfilename(self, ini):
        """
        To specify the filename where the eigenvectors will 
        be written if **dumpfiles flag** is equal to 1.
        """

        if not isinstance(ini, str):
            raise TypeError("set_vctfilename: input must be a string")


        self.__vctfilename = ini

    def get_vctfilename(self):
        """
        *get_vctfilename* return the eigenvectors filename
        """

        return self.__vctfilename

    def set_ovapfilename(self, ini):
        """
        To specify the filename where the overlap matrix will 
        be written if **dumpfiles flag** is equal to 1.
        """

        if not isinstance(ini, str):
            raise TypeError("set_ovapfilename: input must be a string")


        self.__ovapfilename = ini

    def get_ovapfilename(self):
        """
        get_ovapfilename returns the overlap matrix filename 
        """

        return self.__ovapfilename

    def set_fnameinput (self, ini):
        """
        To specify the BERTHA input filename. 
        """

        if not isinstance(ini, str):
            raise TypeError("set_fnameinput: input must be a string")


        self.__fnameinput = ini

    def get_fnameinput (self):
        """
        get_fnameinput returns the BERTHA input filename.
        """

        return self.__fnameinput

    def set_fittfname (self, ini):
        """
        To specify the density fitting basis set filename. 
        """

        if not isinstance(ini, str):
            raise TypeError("set_fittfname: input must be a string")


        self.__fittfname = ini

    def get_fittfname (self):
        """
        get_fittcoefffname returns the density fitting 
        basis set filename.
        """

        return self.__fittfname

    def set_verbosity (self, ini):
        """
        To set the verbosty level:
            -  0 the program will run in a silent modes
            - -1 only basic information are printed-out
            -  1 all the details about the ongoing simlutaion 
               are printed-out
        """

        if not isinstance(ini, int):
            raise TypeError("set_verbosity: input must be an integer")

        self.__verbosity = ini

    def get_verbosity (self):
        """
        get_verbosity return the verbosity level.
        """

        return self.__verbosity

    def set_dumpfiles (self, ini):
        """
        To specify the dumpfiles flag value:
            - if flag is 1  
        """

        if not isinstance(ini, int):
            raise TypeError("set_dumpfiles: input must be an integer")

        self.__dumpfiles = ini

    def get_dumpfiles (self):
        """
        *get_dumpfiles* returns the dumpfiles flag value.
        """

        return self.__dumpfiles

    def init (self):
        """
        Initialize method, it initializes all the variables 
        and basix data. The user need to call this method 
        before he/she can peform the main run.
        """
       
        in_fnameinput = ctypes.c_char_p(self.__fnameinput.encode('utf-8'))

        self.__bertha.init(in_fnameinput, ctypes.c_int(self.__verbosity), 
                ctypes.c_int(self.__dumpfiles))

        self.__init = True

    def get_ndim(self):
        """
        get_ndim returns the matrix dimension.
        """

        if self.__init:
            return self.__bertha.get_ndim()
        
        return 0

    def get_nshift(self):
        """
        get_nshift returns the nshift value.
        """

        if self.__init:
            return self.__bertha.get_nshift()
        
        return 0

    def get_nocc(self):
        """
        get_nocc returns the number of occupied spinors.
        """

        if self.__init:
            return self.__bertha.get_nocc()
        
        return 0

    def get_nopen(self):
        """
        *get_nopen* return number of unoccupied spinors 
        """

        if self.__init:
            return self.__bertha.get_nopen()
        
        return 0

    def get_sfact(self):
        """
        *get_sfact* returns the level shift parameter.
        """

        self.__bertha.get_sfact.restype = ctypes.c_double

        if self.__init:
            return self.__bertha.get_sfact()
        
        return 0

    def density_to_cube_limit(self, dens, fname, ri, rf, 
            drx = 0.2, dry = 0.2, drz = 0.2):
        """
        density_to_cube_limit generates a cube file (named: fname) 
        with the specified density (dens).
        """

        if not (len(ri) == 3):
             raise TypeError("density_to_cube_limit: ri must be list size 3")

        if not (len(rf) == 3):
             raise TypeError("density_to_cube_limit: rf must be list size 3")

        if not isinstance(dens, numpy.ndarray):
            raise TypeError("density_to_cube_limit: input must be a numpy array")

        if not isinstance(fname, str):
            raise TypeError("density_to_cube_limit: input must be a string")

        if not isinstance(drx, float):
            raise TypeError("density_to_cube_limit: input must be a float")

        if not isinstance(dry, float):
            raise TypeError("density_to_cube_limit: input must be a float")

        if not isinstance(drz, float):
            raise TypeError("density_to_cube_limit: input must be a float")

        if self.__init:
            in_fittfname = ctypes.c_char_p(self.__fittfname.encode('utf-8'))
            in_fname = ctypes.c_char_p(fname.encode('utf-8'))
            
            cbuffer = complexmat_to_doublevct (dens)
            
            self.__bertha.density_to_cube_limit (ctypes.c_void_p(cbuffer.ctypes.data), \
                ctypes.c_double(ri[0]), ctypes.c_double(ri[1]), ctypes.c_double(ri[2]), \
                ctypes.c_double(rf[0]), ctypes.c_double(rf[1]), ctypes.c_double(rf[2]), \
                ctypes.c_double(drx), ctypes.c_double(dry), \
                ctypes.c_double(drz), in_fname, in_fittfname)

    def density_to_cube(self, dens, fname, margin = 10.0, 
            drx = 0.2, dry = 0.2, drz = 0.2):
        """
        density_to_cube generates a cube file (named: fname) 
        with the specified density (dens).
        """

        if not isinstance(dens, numpy.ndarray):
            raise TypeError("density_to_cube: input must be a numpy array")

        if not isinstance(fname, str):
            raise TypeError("density_to_cube: input must be a string")

        if not isinstance(margin, float):
            raise TypeError("density_to_cube: input must be a float")

        if not isinstance(drx, float):
            raise TypeError("density_to_cube: input must be a float")

        if not isinstance(dry, float):
            raise TypeError("density_to_cube: input must be a float")

        if not isinstance(drz, float):
            raise TypeError("density_to_cube: input must be a float")

        if self.__init:
            in_fittfname = ctypes.c_char_p(self.__fittfname.encode('utf-8'))
            in_fname = ctypes.c_char_p(fname.encode('utf-8'))
            
            cbuffer = complexmat_to_doublevct (dens)
            
            self.__bertha.density_to_cube (ctypes.c_void_p(cbuffer.ctypes.data), \
                ctypes.c_double(margin), ctypes.c_double(drx), ctypes.c_double(dry), \
                ctypes.c_double(drz), in_fname, in_fittfname)

    def run(self):
        """
        This is the method to perform the SCF computation.
        """

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

            in_fittcoefffname = ctypes.c_char_p(self.__fittcoefffname.encode('utf-8'))
            in_vctfilename = ctypes.c_char_p(self.__vctfilename.encode('utf-8'))
            in_ovapfilename = ctypes.c_char_p(self.__ovapfilename.encode('utf-8'))
            in_fittfname = ctypes.c_char_p(self.__fittfname.encode('utf-8'))

            start = time.time()
            cstart = time.process_time()

            if (self.__childthread):
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
            else:
                self.__bertha.mainrun (in_fittcoefffname, \
                    in_vctfilename, \
                    in_ovapfilename, \
                    in_fittfname, \
                    ctypes.c_void_p(eigen.ctypes.data), \
                    ctypes.c_void_p(ovapbuffer.ctypes.data), \
                    ctypes.c_void_p(eigenvctbu.ctypes.data), \
                    ctypes.c_void_p(fockbuffer.ctypes.data))

            end = time.time()
            cend = time.process_time()

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
        """
        get_etotal returns etotal so that etotal-(sfact*nocc) 
        is equal to the total electronic energy.
        """

        self.__bertha.get_etotal.restype = ctypes.c_double

        if self.__mainrundone:
            return self.__bertha.get_etotal()

        return 0.0

    def get_erep(self):
        """
        It returns nuclear repulsion energy.
        """

        self.__bertha.get_erep.restype = ctypes.c_double

        if self.__mainrundone:
            return self.__bertha.get_erep()

        return 0.0

    def get_eecoul(self):
        """
        It returns coulomb electronic repulsion energy from fitting.
        """

        self.__bertha.get_eecoul.restype = ctypes.c_double

        if self.__mainrundone:
            return self.__bertha.get_eecoul()

        return 0.0

    def get_eexc(self):
        """
        It returns exchange-correlation energy.
        """

        self.__bertha.get_eexc.restype = ctypes.c_double

        if self.__mainrundone:
            return self.__bertha.get_eexc()

        return 0.0



    def realtime_init(self):
        """
        realtime_init initializes all the data needed by the 
        get_realtime_fock method.
        """

        if self.__init and self.__mainrundone:
            self.__bertha.realtime_init()
            self.__realtime_init = True

            ndim = self.get_ndim()

            self.__fockbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
            self.__fockbuffer = numpy.ascontiguousarray(self.__fockbuffer, \
                    dtype=numpy.double)


        else:
            raise Error("You should firstly init and perform the run")

    def get_realtime_dipolematrix (self, direction, normalise):
        """
        get_realtime_dipolematrix computes the dipolematrix, 
        the direction can be:
            - 0 = all three directions
            - 2 = z direction
            - 3 = y direction
            - 4 = x direction
        """

        if not isinstance(direction, int):
            raise TypeError("get_realtime_dipolematrix: input must be an integer")

        if (direction != 0) and (direction != 2) and (direction != 3) and \
                (direction != 4):
            raise TypeError("get_realtime_dipolematrix: input must be 0 or 2 or 3 or 4")

        if not isinstance(normalise, int):
            raise TypeError("get_realtime_dipolematrix: input must be an integer")

        if self.__realtime_init:

            ndim = self.get_ndim()

            vextbuffer = numpy.zeros((2*ndim*ndim), dtype=numpy.double)
            vextbuffer = numpy.ascontiguousarray(vextbuffer, dtype=numpy.double)

            if direction == 0:

                self.__bertha.realtime_dipolematrix(4, normalise, \
                    ctypes.c_void_p(vextbuffer.ctypes.data))

                vextmx = doublevct_to_complexmat (vextbuffer, ndim)
                if vextmx is None:
                    raise Error("Error in vextm matrix size")

                self.__bertha.realtime_dipolematrix(3, normalise, \
                    ctypes.c_void_p(vextbuffer.ctypes.data))

                vextmy = doublevct_to_complexmat (vextbuffer, ndim)
                if vextmy is None:
                    raise Error("Error in vextm matrix size")

                self.__bertha.realtime_dipolematrix(2, normalise, \
                    ctypes.c_void_p(vextbuffer.ctypes.data))

                vextmz = doublevct_to_complexmat (vextbuffer, ndim)
                if vextmz is None:
                    raise Error("Error in vextm matrix size")

                return (vextmx, vextmy, vextmz)

            else:
                self.__bertha.realtime_dipolematrix(direction, normalise, \
                    ctypes.c_void_p(vextbuffer.ctypes.data))

                vextm = doublevct_to_complexmat (vextbuffer, ndim)
                if vextm is None:
                    raise Error("Error in vextm matrix size")
            
            return vextm

        else:
            return None


    def get_eps (self, x, y, z):
        """
        """

        if not isinstance(x, float):
            raise TypeError("get_eps: input must be a float")

        if not isinstance(y, float):
            raise TypeError("get_eps: input must be a float")

        if not isinstance(z, float):
            raise TypeError("get_eps: input must be a float")

        if self.__init:

            self.__bertha.eps.restype = ctypes.c_double

            eps = self.__bertha.eps(ctypes.c_double(x), ctypes.c_double(y), \
                    ctypes.c_double(z))
           
            return eps

        else:
            return None


    def get_coords (self, i):

        if not isinstance(i, int):
            raise TypeError("get_coord: input must be an integer")

        if self.__init:

            vec = numpy.zeros(4, dtype=numpy.double)
            vec = numpy.ascontiguousarray(vec, \
                    dtype=numpy.double)

            self.__bertha.get_coord(i, ctypes.c_void_p(vec.ctypes.data))

            return vec[0], vec[1], vec[2], vec[3]

        return None


    def get_realtime_fock (self, densm):
        """
        get_realtime_fock returns the Fock matrix given the 
        density matrix densm.
        """

        if not isinstance(densm, numpy.ndarray):
            raise TypeError("get_realtime_fock: input must be a numpy array")

        if self.__realtime_init:

            ndim = self.get_ndim()

            cbuffer = complexmat_to_doublevct (densm)

            start = time.time()
            cstart = time.process_time()

            self.__bertha.realtime_fock(ctypes.c_void_p(cbuffer.ctypes.data), \
                    ctypes.c_void_p(self.__fockbuffer.ctypes.data))

            end = time.time()
            cend = time.process_time()

            self.__focktime = end - start
            self.__fockctime = cend - cstart

            fockm = doublevct_to_complexmat (self.__fockbuffer, ndim)
            if fockm is None:
                raise Error("Error in ovap matrix size")

            return fockm

        return None

    def finalize(self):
        """
        To finalize and free all the allocated memory.
        """

        if self.__realtime_init:
            self.__bertha.realtime_finalize()
        
        if self.__init:
            self.__bertha.finalize()

        self.__reset()
