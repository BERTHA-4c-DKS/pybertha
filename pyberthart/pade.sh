#!/bin/sh

export OMP_SCHEDULE=dynamic
export OMP_STACKSIZE=500M
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

ulimit -s unlimited

export BerthaRootPath=/home/redo/BERTHA/

export PYBERTHA_MOD_PATH=$BerthaRootPath"pybertha/pyemb;"$BerthaRootPathxcfun"xcfun/build/lib/python;"$BerthaRootPathx"pybertha/src;"$BerthaRootPath"pyadf/src;"$BerthaRootPath"berthaingen/pybgen;"$BerthaRootPathxcfun"pybertha/psi4rt;"$BerthaRootPath"pybertha/pyberthaemb;"$BerthaRootPath"pybertha/pyberthaembed"

export LD_LIBRARY_PATH=${PWD%/*}"/lib":$LD_LIBRARY_PATH

# limit must be < number of point in the dipole file

export n=1
for gamma in  1.0e-2 2.0e-2 1.0e-3 1.5e-3 2.0e-3
do
   for dw in  0.02 0.002 0.01 0.1 
   do
	echo $gamma $dw
	python /home/redo/BERTHA/Exanalysis/misc/pade/pade_transform.py --limit 25000 --gamma $gamma -f ./omp.mklparallel.quad.2/dipole.txt --dw $dw --fmin 1.5 --frequency 10.0 -o zz_w.txt > pade.txt
	sed "s/dw/$dw/" plot.p > tmp
	sed "s/gamma/$gamma/" tmp > plot.p$n
	gnuplot < plot.p$n
	mv 1.pdf  omp_mklparallel_quad_2_"$n".pdf
	rm -f  plot.p$n tmp
	n=$((n+1))
   done
done
