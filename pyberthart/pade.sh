#!/bin/bash

export OMP_SCHEDULE=dynamic
export OMP_STACKSIZE=500M
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

ulimit -s unlimited

export BerthaRootPath=/home/redo/BERTHA/

export PYBERTHA_MOD_PATH=$BerthaRootPath"pybertha/pyemb;"$BerthaRootPathxcfun"xcfun/build/lib/python;"$BerthaRootPathx"pybertha/src;"$BerthaRootPath"pyadf/src;"$BerthaRootPath"berthaingen/pybgen;"$BerthaRootPathxcfun"pybertha/psi4rt;"$BerthaRootPath"pybertha/pyberthaemb;"$BerthaRootPath"pybertha/pyberthaembed"

export LD_LIBRARY_PATH=${PWD%/*}"/lib":$LD_LIBRARY_PATH

export XDIPDIR=./omp.mklparallel.quad.x/dipole.txt
export YDIPDIR=./omp.mklparallel.quad.y/dipole.txt
export ZDIPDIR=./omp.mklparallel.quad.z/dipole.txt
export PADEPY=/home/redo/BERTHA/Exanalysis/misc/pade/pade_transform.py

 #limit must be < number of point in the dipole file

plotting_func () {
  paste ./xx_w.txt ./yy_w.txt ./zz_w.txt > paste.txt
  sed "s/dw/$dw/" plot.p > tmp
  sed "s/dump/$dump/" tmp > tmp1
  sed "s/gamma/$gamma/" tmp1 > plot.p$n
  gnuplot < plot.p$n
  mv 1.pdf  omp_mklserial_quad_2_"$n"_"$dump".pdf
  rm tmp tmp1 plot.p$n paste.txt
}

export n=1
for dump in gauss 
do
  for gamma in  1.0e-6 1.5e-5 2.0e-6 2.0e-5 1.0e-7 
  do
    for dw in  0.02 0.002 0.01 0.1 
    do
      echo $gamma $dw $dump
      python $PADEPY --damping $dump --limit 25000 --gamma $gamma \
	      -f $XDIPDIR --dw $dw --fmin 2.5 --frequency 6.5 -o xx_w.txt > padex.txt
      python $PADEPY --damping $dump --limit 25000 --gamma $gamma \
	      -f $YDIPDIR --dw $dw --fmin 2.5 --frequency 6.5 -o yy_w.txt > padey.txt
      python $PADEPY --damping $dump --limit 25000 --gamma $gamma \
	      -f $ZDIPDIR --dw $dw --fmin 2.5 --frequency 6.5 -o zz_w.txt > padez.txt
      plotting_func
      n=$((n+1))
    done
  done
done

for dump in exp 
do
  for gamma in  1.0e-3 1.5e-3 2.0e-3 1.0e-2 2.0e-2 
  do
    for dw in  0.02 0.002 0.01 0.1 
    do
      echo $gamma $dw $dump
      python $PADEPY --damping $dump --limit 25000 --gamma $gamma \
	      -f $XDIPDIR --dw $dw --fmin 2.5 --frequency 6.5 -o xx_w.txt > padex.txt
      python $PADEPY --damping $dump --limit 25000 --gamma $gamma \
	      -f $YDIPDIR --dw $dw --fmin 2.5 --frequency 6.5 -o yy_w.txt > padey.txt
      python $PADEPY --damping $dump --limit 25000 --gamma $gamma \
	      -f $ZDIPDIR --dw $dw --fmin 2.5 --frequency 6.5 -o zz_w.txt > padez.txt
      plotting_func
      n=$((n+1))
    done
  done
done
