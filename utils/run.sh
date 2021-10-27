export i=-1
while [ $i -ne 15 ]
do         
  i=$(($i+1));         
  echo "$i" ; 
  python3 valuealongaxis.py  fullpot19_pos.dx $i > out_"$i".txt 
done
