FFLAGS=-c -ffixed-line-length-0 -fpic
MEXFLAGS=-v -f ./mexopts.sh -largeArrayDims
CCFLAGS=-c -fpic
FC=gfortran
CC=gcc
MEX=mex

FFILES=anirec.spike.F90
FOBJ=anirec.spike.o

CPRE=refft.c
COBJ=refft.o

CFILES=anirec_gateway.c

MEXBIN=anirec.mexmaci64

all:

	$(CC) $(CCFLAGS) $(CPRE)
	$(FC) $(FFLAGS) $(FFILES)
	$(MEX) $(MEXFLAGS) $(CFILES) $(FOBJ) $(COBJ)

make_f:

	$(FC) $(FFLAGS) $(FFILES)

make_mex:

	$(MEX) $(MEXFLAGS) $(CFILES) $(FOBJ)

clean:

	rm -rf $(FOBJ) $(MEXBIN)

clean_obj:

	rm -rf $(FOBJ)

#end
