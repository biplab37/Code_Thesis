#!/usr/bin/gnuplot
set terminal png
set output 'zero_momentum_temp_dependence1.png'

set title 'Temperature Dependence of Velocity at zero momentum' font 'hack,12'
set xlabel 'k_BT/{v_F{/Symbol L}_0}' font 'hack,11'
set ylabel 'v_{{/Symbol L}=0}(k)/v_F' font 'hack,11' rotate
set tics font ',9'

f(x) = a+b*log(x)
g(x) = c+d*log(x)+e*log(x)*log(x)
fit f(x) 'zero_momentum_temp_dependence.dat' via a,b
fit g(x) 'zero_momentum_temp_dependence.dat' via c,d,e

plot 'zero_momentum_temp_dependence.dat' title 'velocity at k->0',f(x) title 'a + b ln(x)',g(x) title 'c + d ln(x) + e ln(x)^2'
