F2PYC=f2py
F2PYFLAGS=--fcompiler=gnu95 --f90flags="-fopenmp -lgomp"
FC=gfortran
FFLAGS=-g # instrument for debugging
FFLAGS+=-Wall 	# turn on all warnings 
FFLAGS+=-fbounds-check 		# check array bounds
FFLAGS+=-ffpe-trap=invalid,zero,overflow # trap "serious" FP exceptions


%.so: %.f90
	$(F2PYC) $(F2PYFLAGS) -c -m $* $<    

all: fastcube.so

check: fastcube.f90
	$(FC) $(FFLAGS) -c $<

clean: 
	rm *.so
