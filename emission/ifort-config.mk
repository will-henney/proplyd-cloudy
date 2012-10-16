F2PYC=f2py
F2PYFLAGS=--fcompiler=intelem
FC=ifort

FAST_FFLAGS=-O3
FAST_FFLAGS+=-no-prec-div -xHost -ipo
FAST_FFLAGS+=-fp-model fast=2 -fp-speculation=fast # potentially unsafe optimizations

DEBUG_FLAGS=-g # instrument for debugging
DEBUG_FFLAGS+=-C 		# check array bounds
DEBUG_FFLAGS+=-traceback 

GENERAL_FFLAGS+=-warn all 	# turn on all warnings 
GENERAL_FFLAGS+=-fpe0
GENERAL_FFLAGS+=-std95

OMPFLAGS+=-openmp -openmp-report2

FFLAGS=$(GENERAL_FFLAGS) $(FAST_FFLAGS)
