#!/usr/bin/zsh

start=0.00
end=0.10
difference=0.01

rm tempvsvel.dat

for temp in $(seq $start 0.01 $end);
do
        filename=velon0_$temp.dat
		echo $temp '\t' `head -1 $filename | awk '{print$2}'` >> tempvsvel.dat
done
		
