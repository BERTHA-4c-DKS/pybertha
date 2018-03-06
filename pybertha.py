from ctypes import cdll

berthaw = cdll.LoadLibrary('./bertha_wrapper.so')

berthaw.__bertha_wrapper_MOD_bertha_main()
