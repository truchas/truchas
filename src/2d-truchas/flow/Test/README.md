# Two-dimensional flow tests

The integrated tests are organized by the physics they exercise:

- [`stokes`](stokes/README.md) contains tests for the `flow_2d` Stokes
  executable.
- [`ns`](ns/README.md) contains tests for the `ns_2d` Navier--Stokes
  executable.

The dormant `poiseuille_solid_wall` case is an analytic viscous-channel test
for the mainline-style representation in which pure SOLID cells form the
walls.  It is not registered with CTest because the current solid/liquid
no-slip discretization places the effective wall at the solid-cell center,
rather than at the material interface.

The Fortran unit tests for shared flow components remain in this directory.
