Available Examples
===================================

Below is a list of examples demonstrating the use of IRL. In each example directory, a `run_example.sh` bash script exists to run the example. In this script, the path to the compiler to use will need to be given. A test case can then be run by providing the source file, as in

```
./run_example.sh cut_timing.cpp
```
to perform the C++ example `cut_timing`.
  
Each individual example code has a preceeding comment section explaining the example
and its purpose. The list of examples are given below, and can be found in IRL/examples/*.

# Examples:

## Fortran Examples (examples/fortran)
- @ref cutting_methods.f90 : Finding volume moments for intersected volumes using each of the three methods provided in IRL. 

- @ref localized_separator_link.f90 : Calculate volume moments for a polyhedron distributed across a mesh.

- @ref reconstruction_elvira.f90 : Setup an ELVIRANeighborhood object and use it to perform the ELVIRA interface reconstruction.  

- @ref reconstruction_mof.f90 : Perform a MoF reconstruction with the default and user-provided weights.

- @ref separate_dodecahedron.f90 : Construct a dodecahedron object and calculate the moments internal and external to a thin-sheet intersecting the dodecahedron.

- @ref separate_polygon.f90 : Demonstration of volume moment calculations on polygons.  


\example fortran_interface_example.f90
\example tet_cutting.f90
\example rectilinear_cell_cutting.f90
\example polygon.f90
