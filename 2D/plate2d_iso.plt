
set contour 
unset surface 

set xrange [0:0.3]
set yrange [0:0.4]
set title "Temperature Distribution [*C] in 2-D Plate\n30 x 40 Grid"
set xlabel 'Distance in x-direction [m]'
set ylabel 'Distance in y-direction [m]'

set view map
unset key
unset surface
splot "plate_temp.dat" using 1 with lines

