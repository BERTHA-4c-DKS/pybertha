# basic internal info to install and run using anaconda python and GNU fortran

$ git clone git@bitbucket.org:lstorchi/bertha_ng.git
$ cd bertha_ng
$ make


$ git clone git@github.com:lstorchi/pybertha.git
$ cd pybertha

# set correct BERTHAROOT in config.mk

$ make

# install psi4 

$  curl "http://vergil.chemistry.gatech.edu/psicode-download/Psi4conda-1.3.2-py37-Linux-x86_64.sh" \
    -o Psi4conda-1.3.2-py37-Linux-x86_64.sh --keepalive-time 2
$ bash Psi4conda-1.3.2-py37-Linux-x86_64.sh -b -p $HOME/BERTHAEmb/psi4conda

$ echo $'. $HOME/BERTHAEmb/psi4conda/etc/profile.d/conda.sh\nconda activate' > $HOME/.psi4activate.sh 
$ . $HOME/.psi4activate.sh 
$ psi4 --test

# compile local py3 fork of xcfun

$ git clone git@github.com:lstorchi/xcfun.git
$ cd xcfun 
$ python ./setup --pybindings --cmake-executable=/usr/local/cmake_3.6/bin/cmake
$ cd build 
$ make 

# my py3 fork of the pyadf code

$ git clone git@gitlab.pyadf.org:loriano.storchi/pyadf.gi


# need to edit the source psi4embedrt.py to specify some env variables, like the following:

export PYBERTHA_MOD_PATH="/home/redo/BERTHA/pybertha/pyemb;/home/redo/BERTHA/xcfun/build/lib/python;/home/re
do/BERTHA/pybertha/src;/home/redo/BERTHA/pyadf/src;/home/redo/BERTHA/berthaingen/pybgen;/home/redo/BERTHA/pybertha/psi4rt;/home/redo/BERTHA/pybertha/psi4rt"

python psi4embedrt.py -gA H2O.xyz  -gB environment_cut2.xyz  -d --grid_opts 2 -a 2 --acc_int 4.0 -i -f --sscf



