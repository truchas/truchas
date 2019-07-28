#!/bin/bash -v
mpif90 -c test.f90
mpicc -c scorpiof.c
mpif90 -o test test.o scorpiof.o -L.. -lscorpio -L/soft/hdf5/default -lhdf5 -lz
