#!/bin/bash

ulimit -s unlimited

export OMP_STACKSIZE=1024M
export OMP_NUM_THREADS=6
export OMP_SCHEDULE=dynamic
export LD_LIBRARY_PATH+=:/home/redo/BERTHA/bertha_ng/lib

callfunctions () {
  python3 full_eda_nocv.py --berthamodpaths "../pybertha;../../berthaingen/pybgen"  --fragA ./SHE/$SYSYEMNAME/$ATOMA.xyz  \
     --fragB ./SHE/$SYSYEMNAME/$ATOMB.xyz  --molecule ./SHE/$SYSYEMNAME/$SYSYEMNAME.xyz \
     --basisset "$BASISSET" --fittset "$FITSET" \
     --energyconverter 627.50961 --npairs=12 --convertlengthunit 1.8897259886 --cube \
     --lmargin 10.0 --deltax 0.15 --deltay 0.15 --deltaz 0.15 --externalproces | tee $SYSYEMNAME"out.txt" 

  mv $SYSYEMNAME"out.txt"  nocv_eigv.txt *.cube ./SHE/$SYSYEMNAME/
}

export SYSYEMNAME="AuHg+"
export ATOMB="Hg"
export ATOMA="Au+"
export BASISSET="Au:dyall_vtz,Hg:dyall_vtz"
export FITSET="Au:b20,Hg:hg"

callfunctions


export SYSYEMNAME="AuCn+"
export ATOMB="Cn"
export ATOMA="Au+"
export BASISSET="Au:dyall_vtz,Cn:dyall_vtz"
export FITSET="Au:b20,Cn:cn"

callfunctions

