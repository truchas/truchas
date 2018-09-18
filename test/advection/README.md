VOF ADVECTION TESTS
===================
A collection of tests that exercise fluid vof advection with a prescribed
velocity field. The tests check the initial and final volume fractions for
correctness.

Tests
-----

### 2D single fluid, with diagonally-directed velocity
* test0/advection-0. Coordinate-aligned mesh

### 2D two-fluid, velocity mesh-aligned, perpendicular to interface (mesh1)
* test1a/advection-1a. Coordinate-aligned mesh
* test1b/advection-1b. Rotated 45 degrees about x
* test1c/advection-1c. Rotated 45 degrees about y
* test1d/advection-1d. Rotated 45 degrees about z

### 2D with coordinate aligned mesh (mesh2)
* test2a/advection-2a. Two-fluid, mesh-aligned velocity, oblique interface
* test2b/advection-2b. Two-fluid, diagonally-directed velocity, disk
* test2c/advection-2c. Three-fluid, diagonally-directed velocity, split disk

Meshes
------

### Mesh1
* 31x31x3 domain centered at origin; unit cells

mesh1.gen: the base Cartesian mesh aligned with coordinate axes
mesh1-rotx.gen: mesh1.gen rotated 45 degrees about the x axis
mesh1-roty.gen: mesh1.gen rotated 45 degrees about the y axis
mesh1-rotz.gen: mesh1.gen rotated 45 degrees about the z axis

### mesh2.gen
* 31x31x1 domain, centered at origin; unit cells
* Cartesian mesh aligned with coordinate axes
