#!/usr/bin/zsh
set -x

gfortran bauer_velocity.f90 -o velocity
./velocity

start=0
end=1
filename="frequency_dependence.f90"

for omega in $(seq $start 0.1 $end);
do
	sed -i 's/.*omega =.*/\tomega = '$omega'/' ${filename}
	sed -i 's/.*file="epson0.*/\topen(10,file="epson0_'$omega'.dat")/' ${filename}
	gfortran ${filename} -o output
	./output
done

gnuplot 'plot.gnuplot'
