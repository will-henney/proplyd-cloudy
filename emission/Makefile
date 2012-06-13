# config.mk should be a symbolic link to either gfortran-config.mk or ifort-config.mk
include config.mk
# optionally read in afdditional local configuration
-include local.mk

%.so: %.f90
	$(F2PYC) $(F2PYFLAGS) --f90flags="$(FFLAGS)" -c -m $* $<    

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: fastcube.so fastcube_gauss.so

check: fastcube.o fastcube_gauss.o test_fastcube.o
	rm *.o

fortran_tests: test_fastcube
	echo 1000 180 | time ./test_fastcube

f2py_tests: fastcube.so fastcube_gauss.so
	export PYTHONPATH=$$PWD/../src; cd test_data;  echo $$PYTHONPATH; time python ../model_profile_fast.py --inc 45 --rebin linear-101 --onlylinesin hlines-only.dat

clean: 
	rm *.so *.o *.mod

test_fastcube: test_fastcube.o fastcube.o
	$(FC) $(FFLAGS) -o $@ $^


