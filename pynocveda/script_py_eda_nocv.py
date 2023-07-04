#!/bin/bash

BASEDIR=/home/belp/BERTHAProject
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BASEDIR/bertha_ng/lib

# Set the number of threads to 1 to ensure serial
ulimit -s 4048576
export KMP_STACKSIZE=3400m

# Set the number of threads to 1 to ensure serial
export OMP_NUM_THREADS=4
export OMP_SCHEDULE=dynamic

python3 $BASEDIR/pybertha_private/pynocv/py_eda_nocv.py --thresh  1.0e-11 --wrapperso="$BASEDIR/pybertha_private/lib/bertha_wrapper.so" --berthamodpath="$BASEDIR/pybertha_private/src" -np 4



