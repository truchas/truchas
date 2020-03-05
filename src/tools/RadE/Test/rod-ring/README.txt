Golden output generated with NAG 6.2 chk build:
* mpirun -np 4 genre vf1.inp vf1.re => vf1.001.re vf1.002.re vf1003.re
* mpirun -np 4 genre vf2.inp vf2.re => vf2.cccec10.re, vf2.d3399b7.re, vf2.8bf6188.re

Note that corresponding pairs of output, e.g., vf1.001.re and vf2.cccec10.re,
should be identical when run with the same executable and same number of
MPI processes. 
