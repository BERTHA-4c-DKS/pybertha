#!/bin/bash

ulimit -s unlimited

export OMP_STACKSIZE=1024M
export OMP_NUM_THREADS=6
export OMP_SCHEDULE=dynamic
export LD_LIBRARY_PATH+=:/home/redo/BERTHA/bertha_ng/lib

export SYSYEMNAME="AuCN"
export ATOMA="Au"
export ATOMB="CN"

export BASISSET="Au:dyall_vtz,C:aug-cc-pVTZ-DK,N:aug-cc-pVTZ-DK"
export FITSET="Au:b20,C:fittset,N:fittset"

. ./edafunctions.sh --source-only

#runeda
posteda


