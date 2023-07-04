#!/bin/bash

ulimit -s unlimited

export OMP_STACKSIZE=1024M
export OMP_NUM_THREADS=6
export OMP_SCHEDULE=dynamic
export LD_LIBRARY_PATH+=:/home/redo/BERTHA/bertha_ng/lib

. ./edafunctions.sh --source-only

export SYSYEMNAME="AuHg+"
export ATOMB="Hg"
export ATOMA="Au+"
export BASISSET="Au:dyall_vtz,Hg:dyall_vtz"
export FITSET="Au:b20,Hg:hg"

runeda

export SYSYEMNAME="AuCn+"
export ATOMB="Cn"
export ATOMA="Au+"
export BASISSET="Au:dyall_vtz,Cn:dyall_vtz"
export FITSET="Au:b20,Cn:cn"

runeda

export SYSYEMNAME="AuPb+"
export ATOMB="Pb"
export ATOMA="Au+"
export BASISSET="Au:dyall_vtz,Pb:dyall_vtz"
export FITSET="Au:b20,Pb:pb"

runeda


export SYSYEMNAME="AuFl+"
export ATOMB="Fl"
export ATOMA="Au+"
export BASISSET="Au:dyall_vtz,Fl:dyall_vtz"
export FITSET="Au:b20,Fl:fl"

runeda


export SYSYEMNAME="AuRn+"
export ATOMB="Rn"
export ATOMA="Au+"
export BASISSET="Au:dyall_vtz,Rn:dyall_vtz"
export FITSET="Au:b20,Rn:rn"

runeda


export SYSYEMNAME="AuOg+"
export ATOMB="Og"
export ATOMA="Au+"
export BASISSET="Au:dyall_vtz,Og:dyall_vtz"
export FITSET="Au:b20,Og:og"

runeda


