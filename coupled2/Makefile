run: coupled2.f90
	gfortran coupled2.f90 -o coupled2
	./coupled2

DATAFILE=$(shell echo "*.dat")

plot: $(DATAFILE) plot.gnuplot
	gnuplot plot.gnuplot

clean:
	rm coupled2
	rm fit.log