from ctypes import cdll

import os.path

soname = './bertha_wrapper.so'

if (not os.path.isfile(soname) ):
    print "SO ", soname, " does not exist "
    exit()

berthaw = cdll.LoadLibrary(soname)

berthaw.__bertha_wrapper_MOD_bertha_main()
