import ctypes

import os.path

soname = './bertha_wrapper.so'

if (not os.path.isfile(soname) ):
    print "SO ", soname, " does not exist "
    exit()

berthaw = ctypes.cdll.LoadLibrary(soname)

fittcoefffname = "fitcoeff.txt"
vctfilename = "vct.txt" 
ovapfilename = "ovap.txt"
in_fittcoefffname = ctypes.c_char_p(fittcoefffname)
in_vctfilename = ctypes.c_char_p(vctfilename)
in_ovapfilename = ctypes.c_char_p(ovapfilename)

berthaw.__bertha_wrapper_MOD_bertha_main(in_fittcoefffname, in_vctfilename, \
        in_ovapfilename, len(fittcoefffname), len(vctfilename), len(ovapfilename))
