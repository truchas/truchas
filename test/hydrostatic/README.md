HYDROSTATIC FLUID TESTS
=======================

A collection of tests that exercise inviscid flow with gravitational body force
in various static equilibrium configurations of fluid/void/solid materials. The
tests check the computed initial pressure field for correctness, and check that
the velocity field remains nearly zero with time stepping (no generation of
spurious currents).

Tests
-----

### Single fluid, gravity aligned with mesh
* test1a/hydrostatic-1a. Coordinate-aligned mesh
* test1b/hydrostatic-1b. Rotated 45 degrees about x
* test1c/hydrostatic-1c. Rotated 45 degrees about y
* test1d/hydrostatic-1d. Rotated 45 degrees about z

### Single fluid, gravity directed diagonally across mesh
* test2a/hydrostatic-2a. Coordinate-aligned mesh
* test2b/hydrostatic-2b. Rotated 45 degrees about x
* test2c/hydrostatic-2c. Rotated 45 degrees about y
* test2d/hydrostatic-2d. Rotated 45 degrees about z

### Two fluids, gravity aligned with mesh, lighter fluid on top
* test3a/hydrostatic-3a. Coordinate-aligned mesh
* test3b/hydrostatic-3b. Rotated 45 degrees about x
* test3c/hydrostatic-3c. Rotated 45 degrees about y
* test3d/hydrostatic-3d. Rotated 45 degrees about z

### Fluid/void, gravity aligned with mesh, void on top
* test4a/hydrostatic-4a. Coordinate-aligned mesh
* test4b/hydrostatic-4b. Rotated 45 degrees about x
* test4c/hydrostatic-4c. Rotated 45 degrees about y
* test4d/hydrostatic-4d. Rotated 45 degrees about z

Meshes
------

### Mesh1
* [-5,5]x[-5,5] domain; 10x10 cells, 1 cell thick in z direction
* Side sets 1, 2, 3, 4 on left, right, bottom, top sides, respectively.
* Side set 5 is the symmetry z=const sides
* Domain split into lower half (y<0), block 1; and upper half (y>0) block 2.

mesh1.gen: the base Cartesian mesh aligned with coordinate axes
mesh1-rotx.gen: mesh1.gen rotated 45 degrees about the x axis
mesh1-roty.gen: mesh1.gen rotated 45 degrees about the y axis
mesh1-rotz.gen: mesh1.gen rotated 45 degrees about the z axis
