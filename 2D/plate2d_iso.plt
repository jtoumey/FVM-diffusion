set contour base
set cntrparam level incremental 120, 10, 260
unset surface
set table 'cont.dat'
splot 'plate_temp.dat'
unset table 

reset 
set xrange [0:0.3]
set yrange [0:0.4]
set title "Temperature Distribution [*C] in 2-D Plate\n30 x 40 Grid"
set xlabel 'Distance in x-direction [m]'
set ylabel 'Distance in y-direction [m]'
#unset key
set key title
set palette rgbformulae 33,13,10
p 'plate_temp.dat' with image, 'cont.dat' w l lt -1 lw 1.5

