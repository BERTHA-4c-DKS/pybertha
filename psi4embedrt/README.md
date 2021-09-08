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

# Download and compile bertha_ng

$ git clone git@github.com:BERTHA-4c-DKS/bertha_ng.git
$ cd bertha_ng 
$ make 

# download berthaingen

$ git clone git@github.com:BERTHA-4c-DKS/berthaingen.git

# in pybertha edit the config.mk to spevify the bertha_ng dir BERTHAROOT=/home/redo/BERTHA/bertha_ng
# you need also to export the LD_LIBRARY_PATH to the SOs whithin bertha_ng

$ export LD_LIBRARY_PATH+=:/home/redo/BERTHA/bertha_ng/lib

# now export path of the various modules

$ export PYBERTHA_MOD_PATH="/home/redo/BERTHA/pybertha/pyemb;/home/redo/BERTHA/xcfun/build/lib/python;/home/redo/BERTHA/pybertha/src;/home/redo/BERTHA/pyadf/src;/home/redo/BERTHA/berthaingen/pybgen;/home/redo/BERTHA/pybertha/psi4rt;/home/redo/BERTHA/pybertha/psi4rt"

# no realtime:

python3 psi4embedrt.py -gA H2O.xyz -gB NH3.xyz  --jumprt

# for realtime you need also an input.inp file 

python psi4embedrt.py -gA H2O.xyz  -gB environment_cut2.xyz  -d --grid_opts 2 -a 2 --acc_int 4.0 -i 

# to activate other env, if needed.

$ conda activate p4env_oldversion

