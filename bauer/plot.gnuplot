#!/usr/bin/gnuplot

set terminal png
set output 'dielectric_function3d.png'

set pm3d
set hidden3d
set title 'Dielectric Function' font 'hack,12'
set xlabel '{/Symbol L}/{/Symbol L}_0' font 'hack,11'
set ylabel '  q/{/Symbol L}_0' font 'hack,11'
set zlabel '{/Symbol e}_{{/Symbol L}}(q)' font 'hack,11' rotate
set tics font ',9'

splot 'eps.dat' title 'dielectric function' w l


set output 'velocity3d.png'

set title 'Velocity' font 'hack,12'
set xlabel '  {/Symbol L}/{/Symbol L}_0' font 'hack,11'
set ylabel ' k/{/Symbol L}_0' font 'hack,11'
set zlabel 'v_{/Symbol L}(k)/v_F' font 'hack,11' rotate
set tics font ',9'
set view 60,120,1

splot 'vel.dat' title 'velocity' w l

unset pm3d
unset hidden3d

set output 'renormalised_velocity.png'

f(x) = a + b*log(x)
fit f(x) 'velon0.dat' via a,b
set title 'Renormalised Velocity' font 'hack,12'
set xlabel 'k/{/Symbol L}_0' font 'hack,11'
set ylabel 'v_{{/Symbol L}=0}(k)/v_F' font 'hack,11'
set tics font ',9'

plot 'velon0.dat' title 'velocity', f(x) title 'fitted curve'


set output 'renormalised_dielectric.png'
set title 'Renormalised Velocity' font 'hack,12'
set xlabel 'q/{/Symbol L}_0' font 'hack,11'
set ylabel '{/Symbol e}_{{/Symbol L}=0}(q)' font 'hack,11'
set tics font ',9'

plot 'epson0.dat' title 'dielectric function'
