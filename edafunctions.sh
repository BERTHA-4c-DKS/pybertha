
runeda () {

  #TODO dump also density 

  python3 full_eda_nocv.py --berthamodpaths "../pybertha;../../berthaingen/pybgen"  --fragA ./SHE/$SYSYEMNAME/$ATOMA.xyz  \
     --fragB ./SHE/$SYSYEMNAME/$ATOMB.xyz  --molecule ./SHE/$SYSYEMNAME/$SYSYEMNAME.xyz \
     --basisset "$BASISSET" --fittset "$FITSET" \
     --energyconverter 627.50961 --npairs=12 --convertlengthunit 1.8897259886 --cube \
     --lmargin 10.0 --deltax 0.15 --deltay 0.15 --deltaz 0.15 --externalproces | tee $SYSYEMNAME"out.txt" 

  mv $SYSYEMNAME"out.txt"  nocv_eigv.txt *.cube ./SHE/$SYSYEMNAME/
}

posteda () {

  cp ./SHE/plot.p ./SHE/$SYSYEMNAME/
  cp ./SHE/extractintE.py ./SHE/$SYSYEMNAME/
  cp load_cube.py ./SHE/$SYSYEMNAME/
  cd ./SHE/$SYSYEMNAME/

  # TODO
  # need to call pydens_iso.py an so to cd give the isodensity point 
  export ISODENVAL=$(python3 /home/redo/BERTHA/pycubescd/PY3/pydens_iso.py \
	  -f1 fragmentAdensity.cube  -f2 fragmentBdensity.cube \
	  -axis $AXIS -iseed $SEEDVAL | grep isodensity_point | awk '{print $3}')
  
  echo "isodensity point " $ISODENVAL

  > CT_iso.dat

  i=1
  while [ $i -ne 13 ]
  do
	  echo "RUN: " $i

	  python3 /home/redo/BERTHA/pycubescd/PY3/pycd.py -f pair"$i".cube --axis $AXIS -i $ISODENVAL

	  pymol  -c load_cube.py  pair"$i".cube  $SYSYEMNAME.xyz 0.0014
	  mv this.png pair"$i".png 

	  pymol  -c load_cube.py  "nocv-"$i".cube"  $SYSYEMNAME.xyz 0.0008
	  mv this.png "nocv-"$i".png"

	  pymol  -c load_cube.py  "nocv+"$i".cube"  $SYSYEMNAME.xyz 0.0008
	  mv this.png "nocv+"$i".png"

	  i=$(($i+1))
  done

  pymol  -c load_cube.py diff_tot_ortho.cube $SYSYEMNAME.xyz 0.0014
  mv this.png diff_tot_ortho.png
  python3 /home/redo/BERTHA/pycubescd/PY3/pycd.py -f diff_tot_ortho.cube --axis $AXIS -i $ISODENVAL

  pymol  -c load_cube.py diff_tot.cube $SYSYEMNAME.xyz 0.0014
  mv this.png diff_tot.png
  python3 /home/redo/BERTHA/pycubescd/PY3/pycd.py -f diff_tot.cube --axis $AXIS -i $ISODENVAL

  gnuplot plot.p
  epstopdf cd.ps 

  cp ../results.tex .

  python3 extractintE.py $SYSYEMNAME"out.txt" >> results.tex

  sed -i "s/Bond\ analysis\ SHE/Bond\ analysis\ SHE\ $SYSYEMNAME/g" results.tex

  pdflatex results.tex

  cp results.pdf ../../results$SYSYEMNAME.pdf

  cd -
}


