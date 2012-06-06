F2PYC=f2py
F2PYFLAGS=--fcompiler=gnu95
FC=gfortran
FFLAGS=-g # instrument for debugging
FFLAGS+=-Wall 	# turn on all warnings 
FFLAGS+=-fbounds-check 		# check array bounds
FFLAGS+=-ffpe-trap=invalid,zero,overflow # trap "serious" FP exceptions


%.so: %.f90
	$(F2PYC) $(F2PYFLAGS) --f90flags="$(FFLAGS)" -c -m $* $<    

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: fastcube.so fastcube_gauss.so

check: fastcube.o fastcube_gauss.o
	rm *.o

clean: 
	rm *.so
