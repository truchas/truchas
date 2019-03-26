PIPE FLOW TESTS
===============

## Tests

Viscous 2D pipe flow (flow between infinite plates). Check late-time
velocity profile against analytic parabolic steady state solution. (mesh1a)

* test1a/pipe-flow-old-1a. Coordinate-aligned problem
* test1b/pipe-flow-old-1b. Rotated 45 degrees about x
* test1c/pipe-flow-old-1c. Rotated 45 degrees about y
* test1d/pipe-flow-old-1d. rotated 45 degrees about z

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
