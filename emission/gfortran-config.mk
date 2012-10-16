F2PYC=f2py
F2PYFLAGS=--fcompiler=gnu95
FC=gfortran

FAST_FFLAGS=-O3
FAST_FFLAGS+=-funroll-loops -ffast-math

DEBUG_FFLAGS=-g # instrument for debugging
DEBUG_FFLAGS+=-fbounds-check 		# check array bounds

GENERAL_FFLAGS=-Wall 	# turn on all warnings 
GENERAL_FFLAGS+=-ffpe-trap=invalid,zero,overflow # trap "serious" FP exceptions
GENERAL_FFLAGS+=-pedantic -std=f95

OMPFLAGS=-fopenmp -lgomp

FFLAGS=$(GENERAL_FFLAGS) $(FAST_FFLAGS)
