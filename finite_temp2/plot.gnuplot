#!/usr/bin/gnuplot

set title 'Temperature Dependence of Dielectric Function' font 'hack,12'
set xlabel 'q/{/Symbol L}_0' font 'hack,11'
set ylabel '{/Symbol e}_{{/Symbol L}=0}(q)' font 'hack,11' rotate
set yrange [1:2]

# plot 'velon0_0.00.dat' title 't=0.00' w l
plot 'epson0.dat' title 't=0' w l
replot 'epson0_0.0001.dat' title 't=0.0001' w l
replot 'epson0_0.001.dat' title 't=0.001' w l
replot 'epson0_0.01.dat' title 't=0.01' w l
replot 'epson0_0.02.dat' title 't=0.02' w l
replot 'epson0_0.03.dat' title 't=0.03' w l
replot 'epson0_0.04.dat' title 't=0.04' w l
replot 'epson0_0.05.dat' title 't=0.05' w l
replot 'epson0_0.06.dat' title 't=0.06' w l
replot 'epson0_0.07.dat' title 't=0.07' w l
replot 'epson0_0.08.dat' title 't=0.08' w l
replot 'epson0_0.09.dat' title 't=0.09' w l
replot 'epson0_0.10.dat' title 't=0.10' w l
pause -1
