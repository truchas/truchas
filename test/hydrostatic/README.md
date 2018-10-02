HYDROSTATIC FLUID TESTS
=======================

A collection of tests that exercise inviscid flow with gravitational body force
in various static equilibrium configurations of fluid/void/solid materials. The
tests check the computed initial pressure field for correctness, and check that
the velocity field remains nearly zero with time stepping (no generation of
spurious currents).

Tests
-----

### 2D Single fluid, gravity aligned with mesh (mesh1)
* test1a/hydrostatic-1a. Coordinate-aligned mesh
* test1b/hydrostatic-1b. Rotated 45 degrees about x
* test1c/hydrostatic-1c. Rotated 45 degrees about y
* test1d/hydrostatic-1d. Rotated 45 degrees about z

### 2D Single fluid, gravity directed diagonally across mesh (mesh1)
* test2a/hydrostatic-2a. Coordinate-aligned mesh
* test2b/hydrostatic-2b. Rotated 45 degrees about x
* test2c/hydrostatic-2c. Rotated 45 degrees about y
* test2d/hydrostatic-2d. Rotated 45 degrees about z

### 2D Two fluids, gravity aligned with mesh, lighter fluid on top (mesh1)
* test3a/hydrostatic-3a. Coordinate-aligned mesh
* test3b/hydrostatic-3b. Rotated 45 degrees about x
* test3c/hydrostatic-3c. Rotated 45 degrees about y
* test3d/hydrostatic-3d. Rotated 45 degrees about z

### 2D Fluid/void, gravity aligned with mesh, void on top (mesh1)
* test4a/hydrostatic-4a. Coordinate-aligned mesh
* test4b/hydrostatic-4b. Rotated 45 degrees about x
* test4c/hydrostatic-4c. Rotated 45 degrees about y
* test4d/hydrostatic-4d. Rotated 45 degrees about z

### 2D Two fluids, gravity aligned with mesh, with mixed material cells (mesh2)
* test5a/hydrostatic-5a. Coordinate-aligned mesh
* test5b/hydrostatic-5b. Rotated 45 degrees about x
* test5c/hydrostatic-5c. Rotated 45 degrees about y
* test5d/hydrostatic-5d. Rotated 45 degrees about z

### 2D Fluid/void, gravity aligned with mesh, with mixed material cells (mesh2)
* test6a/hydrostatic-6a. Coordinate-aligned mesh
* test6b/hydrostatic-6b. Rotated 45 degrees about x
* test6c/hydrostatic-6c. Rotated 45 degrees about y
* test6d/hydrostatic-6d. Rotated 45 degrees about z

### 3D multimaterial, mixed material cells, diagonally directed gravity (mesh3)
* test7a/hydrostatic-7a. Two fluids
* test7b/hydrostatic-7b. Fluid/void

### 2D Solid-Fluid, mesh-aligned gravity, perpendicular to solid-fluid interface
* test8a/hydrostatic-8a. Coordinate-aligned mesh
* test8b/hydrostatic-8b. Rotated 45 degrees about x
* test8c/hydrostatic-8c. Rotated 45 degrees about y
* test8d/hydrostatic-8d. Rotated 45 degrees about z

### 2D Solid-Fluid, mesh-aligned gravity, parallel to solid-fluid interface
* test9a/hydrostatic-9a. Coordinate-aligned mesh
* test9b/hydrostatic-9b. Rotated 45 degrees about x
* test9c/hydrostatic-9c. Rotated 45 degrees about y
* test9d/hydrostatic-9d. Rotated 45 degrees about z

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

### Mesh2
* Like mesh1 but 3 cells thick in the z direction
mesh2.gen: the base Cartesian mesh aligned with coordinate axes
mesh2-rotx.gen: mesh1.gen rotated 45 degrees about the x axis
mesh2-roty.gen: mesh1.gen rotated 45 degrees about the y axis
mesh2-rotz.gen: mesh1.gen rotated 45 degrees about the z axis

### Mesh3
* 5^3 domain centered on origin, but rotated so that one of its diagonals is
  aligned with the z-axis.
* 5x5x5 Cartesian mesh (unit cell size)
