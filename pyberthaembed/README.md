rm -rf resultfiles/ 

 python3 pyberthaembed.py -gA H2O.xyz -gB NH3.xyz  \
   --berthamodpaths="/home/redo/BERTHA/pybertha/pyemb;/home/redo/BERTHA/xcfun/build/lib/python;/home/redo/BERTHA/pybertha/src;/home/redo/BERTHA/pyadf/src;/home/redo/BERTHA/berthaingen/pybgen;/home/redo/BERTHA/pybertha/psi4rt" \
    --basisset="" --fittset=""

