import numpy as np
import os

class PyEmbError(Exception):
    """Base class for exceptions in this module."""
    pass


class pyemb:

    # TODO we should use instead a generic molecular data stracuture  
    def __init__ (self, activefname = "", envirofname = ""):
        self.__activefname = activefname
        self.__envirofname = envirofname

        self.__init = False

    def set_active_fname (self, activefname):
        if not isinstance(activefname, str):
            raise TypeError("set_fittcoefffname: input must be a string")

        if (not os.path.isfile(activefname) ):
            raise PyEmbError("File  %s  does not exist" % activefname)

        self.__activefname = activefname
        self.finalize () # we should rerun the init

    def set_enviro_fname (self, envirofname):
        if not isinstance(envirofname, str):
            raise TypeError("set_fittcoefffname: input must be a string")

        if (not os.path.isfile(envirofname) ):
            raise PyEmbError("File  %s  does not exist" % envirofname)

        self.__envirofname = envirofname
        self.finalize () # we should rerun the init

    def get_active_fname (self):
        return self.__activefname

    def get_enviro_fname (self, envirofname):
        return self.__envirofname

    def initialize (self):
        # Here we can inser some init code alla the besic 
        # needed for the proer get_potential start with checking
        # filename maybe we need to specify here the grid to be used ? 
        
        # TODO 

        self.__init = True

    def get_potential (self, density, grid):
        # TODO given the density on the grid get the potential 

        if not isinstance(density,(np.ndarray)):
            raise TypeError("input must be a numpy.ndarray")

        if not isinstance(grid,(np.ndarray)):
            raise TypeError("input must be a numpy.ndarray")
        
        npoints = 0
        
        if len(grid.shape) == 2 and len(density.shape) == 1:
            if grid.shape[1] == 4:
                npoints = density.shape[0]

                if (npoints != density.shape[0]):
                    raise PyEmbError ("incomaptible grid dimension")
            else:
                raise TypeError("input must be a numpy.ndarray npoints X 4")
        else:
            raise TypeError("input must be a numpy.ndarray npoints and npoints X 4 ")
                
        pot = np.zeros(npoints, dtype=np.double)
        pot = np.ascontiguousarray(pot, dtype=np.double)
        
        # TODO 

        return pot

    def finalize (self):
        self.__init = False
        self.__activefname = ""
        self.__envirofname = ""

        # TODO clean all the allocated data 