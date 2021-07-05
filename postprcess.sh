#!/bin/bash


callfunctions () {

  cp load_cube.py ./SHE/$SYSYEMNAME/
  cd ./SHE/$SYSYEMNAME/


  i=1
  while [ $i -ne 13 ]
  do
	  echo "$i"

	  python3 /home/redo/BERTHA/pycubescd/pycd.py -f pair"$i".cube

	  pymol  -c load_cube.py  pair"$i".cube  $SYSYEMNAME.xyz

	  mv this.png pair"$i".png 

	  i=$(($i+1))
  done

  pymol  -c load_cube.py diff_tot_ortho.cube $SYSYEMNAME.xyz
  mv this.png diff_tot_ortho.png
  python3 /home/redo/BERTHA/pycubescd/pycd.py -f diff_tot_ortho.cube

  pymol  -c load_cube.py diff_tot.cube $SYSYEMNAME.xyz
  mv this.png diff_tot.png
  python3 /home/redo/BERTHA/pycubescd/pycd.py -f diff_tot.cube 
  
  cd -
}

export SYSYEMNAME="AuHg+"
callfunctions

export SYSYEMNAME="AuCn+"
callfunctions

export SYSYEMNAME="AuPb+"
callfunctions

export SYSYEMNAME="AuFl+"
callfunctions

export SYSYEMNAME="AuRn+"
callfunctions

export SYSYEMNAME="AuOg+"
callfunctions

