## Makefile
## -----------------------------------------------------
## MPI Version of the Spherical Acoustic Sun Simulator.
## Copyright 2006, Shravan Hanasoge
##
## Hansen Experimental Physics Laboratory
## 455 Via Palou way, Stanford
## CA 94305, USA
## Email: shravan@stanford.edu
## -----------------------------------------------------
##

OBJS1=   analysis.o	dbyd2.o	mtridag.o	dbyd1.o	tridag.o	data_analysis.o	frequencies.o	ziggurat.o
##	bspline90_22.o
###        kernels.o

OBJS2=	postprocess.o

OBJS3=	setup.o

OBJS4=	main.o	kernels.o	dbyd2.o	mtridag.o	dbyd1.o	tridag.o	frequencies.o

FC=	gfortran
FC77=   gfortran

##FC=	gfortran	
##FC77=   gfortran
current_dir = $(shell pwd)
FFLAGS= -O3 -DDOUBLE_PRECISION## -fbounds-check -p -g ##-mcmodel=large ##!-p -g ##-check all ##-fpe0 -traceback -debug #-check bounds
LIBS1=	-L/home/shravan/fftw-3.3.4/lib -lfftw3 -L./ -lcfitsio -L/home/apps/lapack-3.5.0 -llapack -lrefblas -lcurl -L$(current_dir)/optlib -lf90getopt
LIBS2=	-L$(current_dir)/optlib -lf90getopt

COMMAND1=	analyze
COMMAND2=	process
COMMAND3=	setup
COMMAND4=	kernels

INCLUDE= $(current_dir)/opt_lib

$(COMMAND1): $(OBJS1) 
	$(FC) -I$(INCLUDE) $(FFLAGS) -o $(COMMAND1) $(OBJS1) $(LIBS1) 

$(COMMAND2): $(OBJS2) 
	$(FC) -I$(INCLUDE) $(FFLAGS) -o $(COMMAND2) $(OBJS2) $(LIBS1) 

$(COMMAND3): $(OBJS3) 
	$(FC) -I$(INCLUDE) $(FFLAGS) -o $(COMMAND3) $(OBJS3) $(LIBS2)

$(COMMAND4): $(OBJS4) 
	$(FC) -I$(INCLUDE) $(FFLAGS) -o $(COMMAND4) $(OBJS4) $(LIBS1)


%.o : %.f
	$(FC77) -I$(INCLUDE) $(FFLAGS) -c $< 

%.o : %.f90
	$(FC) -I$(INCLUDE) $(FFLAGS) -c $< 

clean:
	@find . \( -name $(COMMAND1) -o -name $(COMMAND2) \) -delete
	@find . \( -name $(COMMAND3) -o -name $(COMMAND4) \) -delete
	@find \( -name "*.o" -o -name "*.mod" \) -delete


data_analysis.o:	params.i	ziggurat.o
kernels.o:	params.i	dbyd2.o	frequencies.o
main.o:       kernels.o
analysis.o:	data_analysis.o	normal.o	f90getopt.mod
frequencies.o:	params.i
dbyd2.o:        mtridag.o
dbyd1.o:        tridag.o
postprocess.o:
setup.o:	f90getopt.mod
normal.o:
mtridag.o:
tridag.o:	
