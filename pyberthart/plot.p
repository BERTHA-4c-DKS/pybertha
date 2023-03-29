set terminal pdf 
set output "1.pdf"
plot "paste.txt"  u ($1*27.2114):($1*($3+$7+$11)) w l lw 2 ti "DW=dw GAMMA=gamma DUMPING=dump"
