COMSOL MICROWAVE OVEN EXAMPLE
-----------------------------

This is a Truchas model of the EM portion of COMSOL's microwave oven example;
see https://www.comsol.com/model/microwave-oven-1424.

Cubit input files for generating meshes at a range of resolutions are in
the meshes directory. Only the coarsest 10 points-per-wavelength mesh is
included. This one runs quickly.

To run:

  mpiexec -np 16 ./fwfd_solver mw-oven.inp

Output in out.vtkhdf
