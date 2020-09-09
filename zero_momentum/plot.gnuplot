#!/usr/bin/gnuplot
set terminal png
set output 'zero_momentum_temp_dependence1.png'

set title 'Temperature Dependence of Velocity at zero momentum' font 'hack,12'
set xlabel 'k_BT/{v_F{/Symbol L}_0}' font 'hack,11'
set ylabel 'v_{{/Symbol L}=0}(k)/v_F' font 'hack,11' rotate
set tics font ',9'

e(x) = a + b*(exp(-c/x) -exp(-d/x))
fit e(x) 'zero_momentum_temp_dependence.dat' via a,b,c,d

plot 'zero_momentum_temp_dependence.dat' title 'velocity at k->0',e(x) title 'a + b ln(x)',g(x) title 'c + d ln(x) + e ln(x)^2'
