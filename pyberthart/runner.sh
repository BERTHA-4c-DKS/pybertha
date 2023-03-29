#!/bin/sh

export OMP_SCHEDULE=dynamic
export OMP_STACKSIZE=500M
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

ulimit -s unlimited

export BerthaRootPath=/home/redo/BERTHA/

export PYBERTHA_MOD_PATH=$BerthaRootPath"pybertha/pyemb;"$BerthaRootPathxcfun"xcfun/build/lib/python;"$BerthaRootPathx"pybertha/src;"$BerthaRootPath"pyadf/src;"$BerthaRootPath"berthaingen/pybgen;"$BerthaRootPathxcfun"pybertha/psi4rt;"$BerthaRootPath"pybertha/pyberthaemb;"$BerthaRootPath"pybertha/pyberthaembed"

export LD_LIBRARY_PATH=${PWD%/*}"/lib":$LD_LIBRARY_PATH

python ./pyberthart.py  --select "-99; 109; 110; 111; 112; 113; 114; 115; 116 & 117;118; 119;120; 121;122" --berthamaxit 500 -d -m 0.08 -g geom.xyz  --restartfile restart_rt.json  --dumprestartnum 3000 --direction 2 --wrapperso $BerthaRootPath/pybertha/lib/bertha_wrapper.so --iterations -T 1000.0  --thresh 1.0e-16 --pulseFmax 0.00005 --pulse analytic --obs "Pb:dyall_2zp,Cl:dyall_2zp" --fittset "Pb:dyall_2zp_autoGENA3spdfg+,Cl:dyall_2zp_autoGENA3spdfg+" --convertlengthunit=1.88972598 -j $BerthaRootPath/berthainputbasis/fullsets.json  1> Pb_res.out 2> Pb_res.err 


python ./pyberthart.py  --select "-99; 109; 110; 111; 112; 113; 114; 115; 116 & 117;118; 119;120; 121;122" --berthamaxit 500 -d -m 0.08 -g geom.xyz  --restartfile restart_rt.json  --dumprestartnum 3000 --direction 3 --wrapperso $BerthaRootPath/pybertha/lib/bertha_wrapper.so --iterations -T 1000.0  --thresh 1.0e-16 --pulseFmax 0.00005 --pulse analytic --obs "Pb:dyall_2zp,Cl:dyall_2zp" --fittset "Pb:dyall_2zp_autoGENA3spdfg+,Cl:dyall_2zp_autoGENA3spdfg+" --convertlengthunit=1.88972598 -j $BerthaRootPath/berthainputbasis/fullsets.json  1> Pb_res.out 2> Pb_res.err 


python ./pyberthart.py  --select "-99; 109; 110; 111; 112; 113; 114; 115; 116 & 117;118; 119;120; 121;122" --berthamaxit 500 -d -m 0.08 -g geom.xyz  --restartfile restart_rt.json  --dumprestartnum 3000 --direction 4 --wrapperso $BerthaRootPath/pybertha/lib/bertha_wrapper.so --iterations -T 1000.0  --thresh 1.0e-16 --pulseFmax 0.00005 --pulse analytic --obs "Pb:dyall_2zp,Cl:dyall_2zp" --fittset "Pb:dyall_2zp_autoGENA3spdfg+,Cl:dyall_2zp_autoGENA3spdfg+" --convertlengthunit=1.88972598 -j $BerthaRootPath/berthainputbasis/fullsets.json  1> Pb_res.out 2> Pb_res.err 

# limit must be < number of point in the dipole file

python /home/redo/BERTHA/Exanalysis/misc/pade/pade_transform.py --limit 12500 --gamma 1.0e-3 -f dipole.txt --dw 0.002 --fmin 1.5 --frequency 7.5 -o zz_w.txt > pade.txt
