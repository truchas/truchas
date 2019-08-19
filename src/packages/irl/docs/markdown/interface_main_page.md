C / Fortran IRL Interface {#mainpage}
================

This document contains a description of all C and Fortran interface functions to IRL. It is primarily written in C with Fortran wrappers provided. To compile the C / Fortran interface, use the command `make opt_interfaces` for an optimized build, or `make debug_interfaces` to use the `DBGFLAGS` defined in `Makefile.in`. 

In order to provide this interface, C functions are mapped to C++ equivalents, often times
with several C functions being used to mimic the behavior of a class. The source files for
the C interface is stored in `src/c_interface`, and each corresponds to a C++ header/implementation
file accept for the prefixed `c_`. Functions that work to mimic the behavior of a C++ class
are formatted so that the class name appears after the prefixed `c_` and before the name
of the C++ method, e.g. `c_ClassName_getId()`. Use of the C interface is handled through
the inclusion of `IRL_c_interface.h` and linking to `libirl.a` and `libirl_c.a`.
Examples using the C interface can be found in the `examples/C` directory.

To use the Fortran IRL interface, the module `irl_fortran_interface`, which is placed in 
`IRL/include`, must be `use`d in the application code.  
This single module provides access to the entire IRL Fortran interface.
Linking must then be performed with `libirl.a`, `libirl_c.a`, and `libirl_fortran.a`.
Fortran derived types are used to wrap C pointers representing the IRL C++ objects. Fortran 2003's `final` keyword is used to provide some form of RAII and help prevent memory leaks. Additionally, Fortran wrappers are written for the C functions in order to provide type and bounds checking before calling the C functions. 
Examples using the Fortran interface can be found in the `examples/fortran` directory.


