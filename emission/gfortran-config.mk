F2PYC=f2py
F2PYFLAGS=--fcompiler=gnu95
FC=gfortran
FFLAGS=-O3
# FFLAGS=-g # instrument for debugging
FFLAGS+=-Wall 	# turn on all warnings 
# FFLAGS+=-fbounds-check 		# check array bounds
FFLAGS+=-ffpe-trap=invalid,zero,overflow # trap "serious" FP exceptions
FFLAGS+=-funroll-loops -ffast-math
#FFLAGS+=-fopenmp -lgomp

