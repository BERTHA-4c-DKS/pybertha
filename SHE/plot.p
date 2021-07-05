


set terminal postscript landscape enhanced color dashed lw 1 "DejaVuSans" 12

set output "cd.ps"


plot "diff_tot.cube_cdz.txt" u 1:2 w l lw 4, "pair1.cube_cdz.txt"  u 1:($2*2) w l lw 4, "pair3.cube_cdz.txt" u 1:(-2*$2) w l lw 4, "pair5.cube_cdz.txt" u 1:(-2*$2) w l lw 4, "pair7.cube_cdz.txt" u 1:(-2*$2) w l lw 4, "pair9.cube_cdz.txt" u 1:(-2*$2) w l lw 4
