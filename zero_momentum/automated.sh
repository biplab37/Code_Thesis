#!/usr/bin/zsh
set -x

start=0
end=0.1
filename="finite.f90"

for temp in $(seq $start 0.01 $end);
do
	sed -i 's/.*temp =.*/\ttemp = '$temp'/' ${filename}
	sed -i 's/.*file="epson0.*/\topen(1,file="epson0_'$temp'.dat")/' ${filename}
	sed -i 's/.*file="velon0.*/\topen(2,file="velon0_'$temp'.dat")/' ${filename}
	gfortran ${filename} -o output
	./output
done

