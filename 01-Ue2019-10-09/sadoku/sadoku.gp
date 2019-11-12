reset
unset key
set logscale x
set xlabel 'kbT / a.u.'
set ylabel 'flips'
set term png
set output "flips.png"
plot 'stats.dat' u 2:3 w l

reset
unset key
set logscale x
set xlabel 'kbT / a.u.'
set ylabel 'energy'
set term png
set output "energy.png"
plot 'stats.dat' u 2:4 w l

reset
unset key
set xlabel 'kbT / a.u.'
set ylabel 'flips '
set term png
set output "flipsraw.png"
plot 'stats.dat' u 1:3 w l

reset
unset key
set xlabel 'kbT / a.u.'
set ylabel 'energy'
set term png
set output "energyraw.png"
plot 'stats.dat' u 1:4 w l