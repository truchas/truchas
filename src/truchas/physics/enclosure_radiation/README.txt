Here is a brief description of each of the source files contained in this
directory.  This should provide a rough overview of how the code and its
functionality is organized.

Neil Carlson <nnc@lanl.gov>
May 2015

================================================================================

rad_system_type.F90
  This defines the rad_system derived type that stores the fixed part of the
  algebraic radiosity system (the view factors), and provides a collection of
  computational procedures associated with the radiosity system.  In these the
  variable data (temperatures, emissivities) are provided as input.

rad_encl_type.F90
  This defines the rad_encl derived type that stores the distributed radiation
  enclosure surface mesh.

rad_encl_func_type.F90
  This defines the rad_encl_func derived type which describes a time-dependent
  function on the enclosure mesh.  The emissivity is an object of this type,
  providing a (possibly different and possibly time-dependent) value for each
  face of the mesh.

rad_solver_type.F90
  This defines the rad_solver derived type that bundles together the radiosity
  system (rad_system) and enclosure surface mesh (rad_encl) along with the
  emissivity and ambient temperature data.  This forms a complete description
  of the radiosity problem, with surface temperature distribution and surface
  radiosity as input and output.  It provides the same collection of
  computational procedures as rad_system (used when embedding the system in a
  larger problem), and a procedure for solving the system for the radiosity.
  This is interface needed for a standalone radiosity problem.

rad_problem_type.F90
  This defines the rad_problem derived type, which provides a bridge between
  the larger 3D heat conduction problem and a radiosity problem that couples
  surface heat fluxes with surface temperatures.  It bundles together a full
  radiosity system (rad_solver) together with data that links enclosure
  surface mesh faces to faces, and exposes the computational procedures
  provided by rad_solver but, with the inputs and outputs in the context of
  the 3D heat transport problem.  The fully coupled 3D heat transport problem
  may involve multiple rad_problem objects, one for each independent enclosure.

ER_file.F90
  This module provides procedures for interacting with the NetCDF radiation
  enclosure disk file, for both reading and writing data.

ER_input.F90
  This module provides the interface to input data.  It provides procedures
  for reading the ENCLOSURE_RADIATION and ENCLOSURE_SURFACE namelists from
  an input file, storing the data as private module data, and then a collection
  of procedures for accessing that data at a later time.

rad_encl_gmv.F90
  This module provides procedures for writing a radiation enclosure surface
  mesh, and face-based fields over that mesh to a GMV-format graphics file.
  It is intended for ad hoc use in standalone test programs and debugging
  situations.

rad_solver_gmv.F90
rad_problem_gmv.F90
  These modules extend the procedures provided by rad_encl_gmv to accept
  rad_solver and rad_problem objects, using the rad_encl object contained
  within.
