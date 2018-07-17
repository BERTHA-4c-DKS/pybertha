set terminal x11 enhanced
k=0.0001
set xrange[0:70]
plot "zw.dat" u ($1*27.2114):($1*4.0/3.0/137.0*pi/k*($3)) w l ti 'bertha-ls-fit-dga2', '../../a1/ls/zw.dat'  u ($1*27.2114):($1*4.0/3.0/137.0*pi/k*($3)) w l ti 'bertha-ls-fit-dga1' , '../../../../h2/zw.dat' u ($1*27.2114):($1*4.0/3.0/137.0*pi/k*($3)) w l ti 'psi4' , '../../../fitt_std/a2/ls/zw.dat' u ($1*27.2114):($1*4.0/3.0/137.0*pi/k*($3)) w l ti 'bertha-ls-fit-std-a2', '../nols/zw.dat' u ($1*27.2114):($1*4.0/3.0/137.0*pi/k*($3)) w l ti 'bertha-nols-fit-dga2',
