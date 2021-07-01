#!/bin/bash

BASEDIR=/home/belp/BERTHAProject
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BASEDIR/bertha_ng/lib

python3 $BASEDIR/pybertha_private/pynocv/py_eda_nocv.py --thresh  1.0e-11 --wrapperso="$BASEDIR/pybertha_private/lib/bertha_wrapper.so" -np 4



