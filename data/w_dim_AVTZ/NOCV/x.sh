for fname in  pair1.cube   ; do

  cp $fname this.cube

  sed -e "s/_xyzfile/geometry/g" \
  -e "s/_colorp/gray90/g" \
  -e "s/_colorm/red/g" \
  -e "s/_label/$fname/g" \
  -e "s/_iso/0.001/g" bp_cube.py > cube.py
  pymol cube.py
  mv this.png $fname.png
  rm this.cube

done
