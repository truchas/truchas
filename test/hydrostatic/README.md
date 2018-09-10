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
* test1a/hydrostatic-old-1a. Coordinate-aligned mesh
* test1b/hydrostatic-old-1b. Rotated 45 degrees about x
* test1c/hydrostatic-old-1c. Rotated 45 degrees about y
* test1d/hydrostatic-old-1d. Rotated 45 degrees about z

### 2D Single fluid, gravity directed diagonally across mesh (mesh1)
* test2a/hydrostatic-old-2a. Coordinate-aligned mesh
* test2b/hydrostatic-old-2b. Rotated 45 degrees about x
* test2c/hydrostatic-old-2c. Rotated 45 degrees about y
* test2d/hydrostatic-old-2d. Rotated 45 degrees about z

### 2D Two fluids, gravity aligned with mesh, lighter fluid on top (mesh1)
* test3a/hydrostatic-old-3a. Coordinate-aligned mesh
* test3b/hydrostatic-old-3b. Rotated 45 degrees about x
* test3c/hydrostatic-old-3c. Rotated 45 degrees about y
* test3d/hydrostatic-old-3d. Rotated 45 degrees about z

### 2D Fluid/void, gravity aligned with mesh, void on top (mesh1)
* test4a/hydrostatic-old-4a. Coordinate-aligned mesh
* test4b/hydrostatic-old-4b. Rotated 45 degrees about x
* test4c/hydrostatic-old-4c. Rotated 45 degrees about y
* test4d/hydrostatic-old-4d. Rotated 45 degrees about z

### 2D Two fluids, gravity aligned with mesh, with mixed material cells (mesh2)
* test5a/hydrostatic-old-5a. Coordinate-aligned mesh
* test5b/hydrostatic-old-5b. Rotated 45 degrees about x
* test5c/hydrostatic-old-5c. Rotated 45 degrees about y
* test5d/hydrostatic-old-5d. Rotated 45 degrees about z

### 2D Fluid/void, gravity aligned with mesh, with mixed material cells (mesh2)
* test6a/hydrostatic-old-6a. Coordinate-aligned mesh
* test6b/hydrostatic-old-6b. Rotated 45 degrees about x
* test6c/hydrostatic-old-6c. Rotated 45 degrees about y
* test6d/hydrostatic-old-6d. Rotated 45 degrees about z

### 3D multimaterial, mixed material cells, diagonally directed gravity (mesh3)
* test7a/hydrostatic-old-7a. Two fluids
* test7b/hydrostatic-old-7b. Fluid/void

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
