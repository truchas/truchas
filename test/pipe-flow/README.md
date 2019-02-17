PIPE FLOW TESTS
===============

## Tests

Viscous 2D pipe flow (flow between infinite plates). Check late-time
velocity profile against analytic parabolic steady state solution. (mesh1a)

* test1a/pipe-flow-1a. Coordinate-aligned problem
* test1b/pipe-flow-1b. Rotated 45 degrees about x
* test1c/pipe-flow-1c. Rotated 45 degrees about y
* test1d/pipe-flow-1d. rotated 45 degrees about z

Same problem but the domain includes a 2-cell layer of solid material along
the wall side, which replaces the explicit no-slip BC on the side. (mesh2a)

* test2a/pipe-flow-2a. Coordinate-aligned problem
* test2b/pipe-flow-2b. Rotated 45 degrees about x
* test2c/pipe-flow-2c. Rotated 45 degrees about y
* test2d/pipe-flow-2d. rotated 45 degrees about z

These are variations of the coordinate-aligned problem test2a that use differing
methods for specifying the no-slip BC along the solid wall.

* pipe-flow-3a. No explicit BC; just what happens at a mesh-conforming
  solid/liquid interface. This is pipe-flow-2a with modified solver parameters.
* pipe-flow-3b. No-slip BC on internal side set.
* pipe-flow-3c. Internal interface created on solid/liquid side set and no-slip
  BC specified on both sides.
* pipe-flow-3d. Internal interface created on solid/liquid side set and no-slip
  BC specified on fluid side of interface

Meshes
------

### Mesh1a
* [-0.5,0.5] x [0,1] x [-0.3,0.3] domain
* Side set 1: x = -0.5 (pipe inlet)
* Side set 2: x = 0.5 (pipe outlet)
* Side set 3: y = 1 (pipe wall)
* Side set 4: z = constant, y = 0 sides (symmetry)

mesh1a.gen: Cartesian 5x10x3 hex mesh aligned with coordinate axes
mesh1a-rotx.gen: mesh1a.gen rotated 45 degrees about the x axis
mesh1a-roty.gen: mesh1a.gen rotated 45 degrees about the y axis
mesh1a-rotz.gen: mesh1a.gen rotated 45 degrees about the z axis

### Mesh2a
* [-0.5,0.5] x [0,1.2] x [-0.3,0.3] domain
* Block 1 is fluid: [-0.5,0.5] x [0,1] x [-0.3,0.3]
* Block 2 is solid: [-0.5,0.5] x [1,1.2] x [-0.3,0.3]
* Side set 1: x = -.0.5 (pipe inlet)
* Side set 2: x = 0.5 (pipe outlet)
* Side set 3: y = 1 (pipe wall, INTERNAL TO THE MESH)
* Side set 4: z = constant, y = 0 sides (symmetry)
* Side set 5: y = 1, relative to the fluid body

mesh2a.gen: Cartesian 5x12x3 hex mesh aligned with coordinate axes
mesh2a-rotx.gen: mesh2a.gen rotated 45 degrees about the x axis
mesh2a-roty.gen: mesh2a.gen rotated 45 degrees about the y axis
mesh2a-rotz.gen: mesh2a.gen rotated 45 degrees about the z axis
