$ conda create -n p4env psi4 -c psi4/label/dev
$ conda activate p4env

$ conda install -c anaconda swig
$ conda install cmake

$ git clone git@github.com:lstorchi/xcfun.git
$ cd xcfun
$ python ./setup --pybindings 
$ cd build
$ make

$ git clone git@gitlab.pyadf.org:loriano.storchi/pyadf.git

$ . /usr/local/adf2019.307/adfbashrc.sh 

# Dowload and compile bertha_ng, download berthaingen nd pybertha

$ export LD_LIBRARY_PATH+=:/home/redo/BERTHA/bertha_ng/lib

# now export path of the various modules

$ export PYBERTHA_MOD_PATH="/home/redo/BERTHA/pybertha/pyemb;/home/redo/BERTHA/xcfun/build/lib/python;/home/redo/BERTHA/pybertha/src;/home/redo/BERTHA/pyadf/src;/home/redo/BERTHA/berthaingen/pybgen;/home/redo/BERTHA/pybertha/psi4rt;/home/redo/BERTHA/pybertha/psi4rt"

python psi4embedrt.py -gA H2O.xyz  -gB environment_cut2.xyz  -d --grid_opts 2 -a 2 --acc_int 4.0 -i 

# no realtime:

python3 psi4embedrt.py -gA H2O.xyz -gB NH3.xyz  --jumprt

# to activate other env, if needed.

$ conda activate p4env_oldversion



