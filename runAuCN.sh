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

callfunctionspost () {

  cp ./SHE/plot.p ./SHE/$SYSYEMNAME/
  cp ./SHE/extractintE.py ./SHE/$SYSYEMNAME/
  cp load_cube.py ./SHE/$SYSYEMNAME/
  cd ./SHE/$SYSYEMNAME/


  i=1
  while [ $i -ne 13 ]
  do
	  echo "$i"

	  python3 /home/redo/BERTHA/pycubescd/pycd.py -f pair"$i".cube

	  pymol  -c load_cube.py  pair"$i".cube  $SYSYEMNAME.xyz
	  mv this.png pair"$i".png 

	  pymol  -c load_cube.py  "nocv-"$i".cube"  $SYSYEMNAME.xyz
	  mv this.png "nocv-"$i".png"

	  pymol  -c load_cube.py  "nocv+"$i".cube"  $SYSYEMNAME.xyz
	  mv this.png "nocv+"$i".png"

	  i=$(($i+1))
  done

  pymol  -c load_cube.py diff_tot_ortho.cube $SYSYEMNAME.xyz
  mv this.png diff_tot_ortho.png
  python3 /home/redo/BERTHA/pycubescd/pycd.py -f diff_tot_ortho.cube

  pymol  -c load_cube.py diff_tot.cube $SYSYEMNAME.xyz
  mv this.png diff_tot.png
  python3 /home/redo/BERTHA/pycubescd/pycd.py -f diff_tot.cube 

  gnuplot plot.p
  epstopdf cd.ps 

  cp ../results.tex .

  python3 extractintE.py $SYSYEMNAME"out.txt" >> results.tex

  pdflatex results.tex

  cd -
}

export SYSYEMNAME="AuCN"
export ATOMA="Au"
export ATOMB="CN"

export BASISSET="Au:dyall_vtz,C:aug-cc-pVTZ-DK,N:aug-cc-pVTZ-DK"
export FITSET="Au:b20,C:fittset,N:fittset"

callfunctions
callfunctionspost


