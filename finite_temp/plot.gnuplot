#!/usr/bin/gnuplot

set title 'Temperature Dependence of Velocity' font 'hack,12'
set xlabel 'k/{/Symbol L}_0' font 'hack,11'
set ylabel 'v_{{/Symbol L}=0}(k)/v_F' font 'hack,11' rotate

plot 'velon0_0.00.dat' title 't=0.00' w l
replot 'velon0_0.01.dat' title 't=0.01' w l
replot 'velon0_0.02.dat' title 't=0.02' w l
replot 'velon0_0.03.dat' title 't=0.03' w l
replot 'velon0_0.04.dat' title 't=0.04' w l
replot 'velon0_0.05.dat' title 't=0.05' w l
replot 'velon0_0.06.dat' title 't=0.06' w l
replot 'velon0_0.07.dat' title 't=0.07' w l
replot 'velon0_0.08.dat' title 't=0.08' w l
replot 'velon0_0.09.dat' title 't=0.09' w l
replot 'velon0_0.10.dat' title 't=0.10' w l
pause -1
