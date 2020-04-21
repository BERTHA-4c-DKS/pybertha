set terminal postscript eps enhanced color font "Helvetica" 18
set output 'water_zz.eps'

set xlabel "{/Palatino-Italic Energy} /{/Palatino-Italic eV}"
set ylabel "{/Palatino-Italic S_z(w)}/arb. units"
set key samplen 1.5
plot 'ref/zz_w.txt' u ($1*27.2114):($1*$3) w lp pt 6 ps 0.4 ti 'ref', 'ref/zz_w.txt' u ($1*27.2114):($1*$3) w l noti
