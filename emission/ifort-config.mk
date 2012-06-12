F2PYC=f2py
F2PYFLAGS=--fcompiler=intelem
FC=ifort
FFLAGS=-O3
FFLAGS=-g # instrument for debugging
FFLAGS+=-warn all 	# turn on all warnings 
FFLAGS+=-C 		# check array bounds
FFLAGS+=-fpe0
#FFLAGS+=-fopenmp -lgomp

