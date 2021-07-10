$ conda create -n p4env psi4 -c psi4/label/dev
$ conda activate p4env

$ conda install -c anaconda swig
$ conda install cmake

$ git clone git@github.com:lstorchi/xcfun.git
$ cd xcfun
$ python ./setup --pybindings --cmake-executable=/usr/local/cmake_3.6/bin/cmake
$ cd build
$ make

$ git clone git@gitlab.pyadf.org:loriano.storchi/pyadf.git

# need to edit the source psi4embedrt.py to specify some env variables

