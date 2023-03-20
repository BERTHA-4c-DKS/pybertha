set terminal pdf 
set output "1.pdf"
plot "zz_w.txt"  u ($1*27.2114):($1*$3) w l lw 2 ti "DW=dw GAMMA=gamma"

