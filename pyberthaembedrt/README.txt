$ conda create -n p4env psi4 -c psi4/label/dev
$ conda activate p4env

$ conda install -c anaconda swig
$ conda install cmake

$ pip3 install mendeleev
$ pip3 install scipy

$ git clone git@github.com:lstorchi/xcfun.git
$ cd xcfun
$ python ./setup --pybindings 

# the pythi it is finding must be the same as the one used fpr the run

$ cd build
$ make

# instead you can use 
$ python ./setup --pybindings --show 
# to see the cmake script and afterwards add lso the PYTHON library path and bin
$ mkdit build
$ cd ./build 
$ cmake -DCMAKE_CXX_COMPILER=g++ -DEXTRA_CXXFLAGS="''" -DCMAKE_C_COMPILER=gcc -DEXTRA_CFLAGS="''" -DCMAKE_Fortran_COMPILER=gfortran -DEXTRA_FCFLAGS="''" -DENABLE_FC_SUPPORT=OFF -DENABLE_CODE_COVERAGE=False -DSTATIC_LIBRARY_ONLY=False -DXCFun_XC_MAX_ORDER="3" -DENABLE_PYTHON_INTERFACE=True -DCMAKE_BUILD_TYPE=debug -G"Unix Makefiles" -DPYTHON_LIBRARIES=/gpfswork/rech/uhf/ums99fa/anaconda3/envs/p4env/lib -DPYTHON_PREFIX=/gpfswork/rech/uhf/ums99fa/anaconda3/envs/p4env/bin/ -DPYTHON_EXECUTABLE=/gpfswork/rech/uhf/ums99fa/anaconda3/envs/p4env/bin/python3 ../ 
$ make 


$ git clone git@gitlab.pyadf.org:loriano.storchi/pyadf.git

$ . /usr/local/adf2019.307/adfbashrc.sh 

# Download and compile bertha_ng

$ git clone git@github.com:BERTHA-4c-DKS/bertha_ng.git
$ cd bertha_ng 
$ make 

# download berthaingen

$ git clone git@github.com:BERTHA-4c-DKS/berthaingen.git

# to compile pyberth need to specify the bertha_ng root dir:

export BerthaRootPath=/home/redo/BERTHA/

# and to run you need also to export the LD_LIBRARY_PATH to the SOs whithin bertha_ng

$ export LD_LIBRARY_PATH+=:/home/redo/BERTHA/bertha_ng/lib

# can use also PYBERTHA_MOD_PATH env variable to specify all PATHS:

$ export PYBERTHA_MOD_PATH="/home/redo/BERTHA/pybertha/pyemb;/home/redo/BERTHA/xcfun/build/lib/python;/home/redo/BERTHA/pybertha/src;/home/redo/BERTHA/pyadf/src;/home/redo/BERTHA/berthaingen/pybgen;/home/redo/BERTHA/pybertha/psi4rt;/home/redo/BERTHA/pybertha/pyberthaemb;/home/redo/BERTHA/xcfun/build/lib64/python/xcfun;/home/redo/BERTHA/pybertha/pyberthaembed"

# now run

python3 pyberthaembedrt.py  -gA H2O.xyz -gB NH3.xyz --act_fittset "H:aux4,N:aux4,O:aux4" --act_obs "H:aug-cc-pvdz,O:aug-cc-pvdz,N:aug-cc-pvdz" --convertlengthunit=1.889725988

