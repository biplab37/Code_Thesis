## This is a Makefile which automates the process of compiling, running of codes as well as producing plots
PROGRAM=gfortran
DATAFILE=$(shell echo *.dat)\

.PHONY: all plot plotgnu

all: bauer

## bauer   : Run the FORTRAN code to generate the data
bauer: bauer.f90
	${PROGRAM} $^ -o $@
	./bauer

## plot    : Plotting command using matplotlib create pdf files suitable for pdflatex.
plot: plot.py $(DATAFILE)
	python3 $<

## plotgnu : Plotting using Gnuplot. 3d plots are easier to handle here
plotgnu: plot.gnuplot $(DATAFILE)
	gnuplot $<

.PHONY: clean

## clean   : Remove extra files
clean:
	rm bauer
	rm fit.log

.PHONY : help

## help    : Shows help.
help : Makefile
	@sed -n 's/^##//p' $<