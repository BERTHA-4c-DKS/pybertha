#!/bin/bash

BASEDIR=/home/belp/BERTHAEDA
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BASEDIR/bertha_ng/lib

# Set the number of threads to 1 to ensure serial
ulimit -s 4048576
export KMP_STACKSIZE=3400m

# Set the number of threads to 1 to ensure serial
export OMP_NUM_THREADS=1
export OMP_SCHEDULE=dynamic

# Set GEOMETRY files FRAG1, FRAG2 and MOLECULE 
frag1xyz=w_1.xyz
frag2xyz=w_2.xyz
moleculexyz=w_dimer.xyz

FRAG1DIR=W1
FRAG2DIR=W2
NOCVDIR=NOCV
excfunct=LXCBP86


# Set BASIS SET 

GspinoBasis="O:aug-cc-pvtz,H:aug-cc-pvtz"
fittingbasis="O:aux4,H:aux4"
basisset=$BASEDIR/berthainputbasis/fullsets.json

# Gen INPUT Frag1 #
mkdir $FRAG1DIR
cd $FRAG1DIR

python3 $BASEDIR/berthaingen/pybgen/pybgen.py -f ../$frag1xyz  -b $GspinoBasis -t $fittingbasis --convertlengthunit 1.8897259886 --functxc $excfunct -j $basisset 

echo "Run pybertha in $FRAG1DIR" 
python3 $BASEDIR/pybertha_private/pybertha/pybertha.py --wrapperso=$BASEDIR/pybertha_private/lib/bertha_wrapper.so --berthamodpath=/home/belp/BERTHAEDA/pybertha_private/src --thresh=0.00000000001 --eda_nocv_info

cd ..
mkdir $FRAG2DIR
cd $FRAG2DIR

python3 $BASEDIR/berthaingen/pybgen/pybgen.py -f ../$frag2xyz  -b $GspinoBasis -t $fittingbasis --convertlengthunit 1.8897259886 --functxc $excfunct -j $basisset


echo "Run pybertha in $FRAG2DIR" 
python3 $BASEDIR/pybertha_private/pybertha/pybertha.py --wrapperso=$BASEDIR/pybertha_private/lib/bertha_wrapper.so --berthamodpath=/home/belp/BERTHAEDA/pybertha_private/src --thresh=0.00000000001 --eda_nocv_info

cd ..
mkdir $NOCVDIR 
cd $NOCVDIR
python3 $BASEDIR/berthaingen/pybgen/pybgen.py -f ../$moleculexyz -b $GspinoBasis -t $fittingbasis --convertlengthunit 1.8897259886 --functxc $excfunct -j $basisset 

cp ../$FRAG1DIR/info_eda_nocv_fragX.json info_eda_nocv_fragA.json
cp ../$FRAG2DIR/info_eda_nocv_fragX.json info_eda_nocv_fragB.json

echo 'Run pyeda in NOCV'
python3  $BASEDIR/pybertha_private/pynocv/py_eda_nocv.py --wrapperso=/home/belp/BERTHAEDA/pybertha_private/lib/bertha_wrapper.so --berthamodpath=/home/belp/BERTHAEDA/pybertha_private/src -np 12 --thresh  0.00000000001 > eda.out

python3 ../extractintE.py eda.out  > ../tabrisu.tex


