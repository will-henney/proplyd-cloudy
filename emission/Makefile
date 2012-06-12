# config.mk should be a symbolic link to either gfortran-config.mk or ifort-config.mk
include config.mk

%.so: %.f90
	$(F2PYC) $(F2PYFLAGS) --f90flags="$(FFLAGS)" -c -m $* $<    

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: fastcube.so fastcube_gauss.so

check: fastcube.o fastcube_gauss.o
	rm *.o

clean: 
	rm *.so
