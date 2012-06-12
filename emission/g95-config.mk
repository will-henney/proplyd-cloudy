F2PYC=f2py
F2PYFLAGS=--fcompiler=g95
FC=g95
FFLAGS=-O3
FFLAGS=-g # instrument for debugging
FFLAGS+=-Wall 	# turn on all warnings 
FFLAGS+=-fbounds-check		# check array bounds
#FFLAGS+=-fopenmp -lgomp

