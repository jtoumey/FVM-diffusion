set terminal pngcairo enhanced font 'Verdana,10'
set output 'difsrc_temp_distr.png'

set title "Temperature Distribution"
set xlabel "Length Along Domain [m]"
set ylabel "Temperature [*C]"
set key center top
plot 'temp_distr.dat' with linespoints