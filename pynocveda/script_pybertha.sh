#!/bin/bash

BASEDIR=/home/belp/BERTHAProject
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BASEDIR/bertha_ng/lib

python3 ~/BERTHAProject/pybertha_private/pybertha/pybertha.py --thresh 1.0e-11 --wrapperso="$BASEDIR/pybertha_private/lib/bertha_wrapper.so" --eda_nocv_info --eda_nocv_frag_file info_eda_nocv_fragX.json




