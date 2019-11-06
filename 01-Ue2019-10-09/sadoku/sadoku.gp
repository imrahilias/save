reset
unset key
set logscale x
set xlabel 'kbT / a.u.'
set ylabel 'flips / pixel / sweep'
set term png
set output "flips.png"
plot 'stats.dat' u 2:($3/81) w l

reset
unset key
set logscale x
set xlabel 'kbT / a.u.'
set ylabel 'energy'
set term png
set output "energy.png"
plot 'stats.dat' u 2:($4) w l