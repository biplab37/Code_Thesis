#!/usr/bin/zsh
set -x

start=0.01
end=0.1
filename="finite2.f90"

	sed -i 's/.*temp =.*/\ttemp = 0.0001/' ${filename}
	sed -i 's/.*file="epson0.*/\topen(3,file="epson0_0.0001.dat")/' ${filename}
	sed -i 's/.*file="velon0.*/\topen(4,file="velon0_0.0001.dat")/' ${filename}
	gfortran ${filename} -o output
	./output

	sed -i 's/.*temp =.*/\ttemp = 0.001/' ${filename}
	sed -i 's/.*file="epson0.*/\topen(3,file="epson0_0.001.dat")/' ${filename}
	sed -i 's/.*file="velon0.*/\topen(4,file="velon0_0.001.dat")/' ${filename}
	gfortran ${filename} -o output
	./output

for temp in $(seq $start 0.01 $end);
do
	sed -i 's/.*temp =.*/\ttemp = '$temp'/' ${filename}
	sed -i 's/.*file="epson0.*/\topen(3,file="epson0_'$temp'.dat")/' ${filename}
	sed -i 's/.*file="velon0.*/\topen(4,file="velon0_'$temp'.dat")/' ${filename}
	gfortran ${filename} -o output
	./output
done

