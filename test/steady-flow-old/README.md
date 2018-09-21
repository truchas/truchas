STEADY INVISCID MULTI-MATERIAL FLOW TESTS
=========================================
A collection of tests that exercise fluid vof advection with steady inviscid
flow problems. The tests check the initial and final volume fractions for
correctness, as well as the pressure and velocity fields.

Tests
-----

### 2D two-fluid, equal density, velocity mesh-aligned, perpendicular to interface (mesh1a)
* test1a/steady-flow-old-1a. Coordinate-aligned mesh
* test1b/steady-flow-old-1b. Rotated 45 degrees about x
* test1c/steady-flow-old-1c. Rotated 45 degrees about y
* test1d/steady-flow-old-1d. Rotated 45 degrees about z

### 2D two-fluid, 1000:1 density ratio, heavy fluid trailing, velocity mesh-aligned,
    perpendicular to interface (mesh1a)
* test2a/steady-flow-old-2a. Coordinate-aligned mesh
* test2b/steady-flow-old-2b. Rotated 45 degrees about x
* test2c/steady-flow-old-2c. Rotated 45 degrees about y
* test2d/steady-flow-old-2d. Rotated 45 degrees about z

### 2D two-fluid, 200:1 density ratio, heavy fluid leading, velocity mesh-aligned,
    perpendicular to interface (mesh1a)
* test3a/steady-flow-old-3a. Coordinate-aligned mesh
* test3b/steady-flow-old-3b. Rotated 45 degrees about x
* test3c/steady-flow-old-3c. Rotated 45 degrees about y
* test3d/steady-flow-old-3d. Rotated 45 degrees about z

### 2D two-fluid, equal density, velocity mesh-aligned, oblique interface (mesh2a)
* test4a/steady-flow-old-4a. Coordinate-aligned mesh
* test4b/steady-flow-old-4b. Rotated 45 degrees about x
* test4c/steady-flow-old-4c. Rotated 45 degrees about y
* test4d/steady-flow-old-4d. Rotated 45 degrees about z

### 2D two-fluid, 1000:1 density ratio, heaving fluid trailing, velocity mesh-aligned,
    oblique interface (mesh2a)
* test5a/steady-flow-old-5a. Coordinate-aligned mesh
* test5b/steady-flow-old-5b. Rotated 45 degrees about x
* test5c/steady-flow-old-5c. Rotated 45 degrees about y
* test5d/steady-flow-old-5d. Rotated 45 degrees about z

### test6a/steady-flow-old-6a (mesh3a)
* 2D, velocity directed diagonally across mesh
* Two-fluid, 100:1 density ratio.
* heavy fluid in disk; light fluid in complement

### test6b/steady-flow-old-6b (mesh3a)
* 2D, velocity directed diagonally across mesh
* Three-fluid, 100:1 density ratio.
* Two heavy fluids in either half of disk; light fluid in complement

Meshes
------

### Mesh1a
* 13x13x1 domain centered at origin; unit cells

mesh1a.gen: the base Cartesian mesh aligned with coordinate axes
mesh1a-rotx.gen: mesh1.gen rotated 45 degrees about the x axis
mesh1a-roty.gen: mesh1.gen rotated 45 degrees about the y axis
mesh1a-rotz.gen: mesh1.gen rotated 45 degrees about the z axis

### Mesh2a
* 19x13x1 domain, centered at origin; unit cells

mesh2a.gen: the base Cartesian mesh aligned with coordinate axes
mesh2a-rotx.gen: mesh1.gen rotated 45 degrees about the x axis
mesh2a-roty.gen: mesh1.gen rotated 45 degrees about the y axis
mesh2a-rotz.gen: mesh1.gen rotated 45 degrees about the z axis

### mesh3a.gen. 31x31x1 domain, centered at origin, unit cells
