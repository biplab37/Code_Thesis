set output 'frequency.png'

set title 'Frequency Dependence of Dielectric Function' font ',12'
set xlabel 'q/{/Symbol L}_0' font 'hack,11'
set ylabel '{/Symbol e}_{{/Symbol L}=0}(q)' font 'hack,11' rotate
set tics font ',9'

#plot 'epson0_0.0.dat' title 'omega=0' w l
plot 'epson0_0.1.dat' title 'omega=0.1' w l
replot 'epson0_0.2.dat' title 'omega=0.2' w l
replot 'epson0_0.3.dat' title 'omega=0.3' w l
replot 'epson0_0.4.dat' title 'omega=0.4' w l
replot 'epson0_0.5.dat' title 'omega=0.5' w l
replot 'epson0_0.6.dat' title 'omega=0.6' w l
replot 'epson0_0.7.dat' title 'omega=0.7' w l
replot 'epson0_0.8.dat' title 'omega=0.8' w l
replot 'epson0_0.9.dat' title 'omega=0.9' w l
replot 'epson0_1.0.dat' title 'omega=1.0' w l
pause -1
