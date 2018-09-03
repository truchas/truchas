PIPE FLOW TESTS
===============

## Tests

Inviscid flow in pipe with square cross section. Free-slip BC on pipe walls.
Fluid accelerates from rest (constant acceleration). Check uniform velocity
value after a couple steps and at final time. (mesh1a)

* test1/pipe-flow-old-inviscid-1. Coordinate-aligned problem
* test2/pipe-flow-old-inviscid-2. Rotated 45 degrees about x
* test3/pipe-flow-old-inviscid-3. Rotated 45 degrees about y
* test4/pipe-flow-old-inviscid-4. rotated 45 degrees about z

Usual viscous 2D pipe flow (flow between infinite plates). Check late-time
velocity profile against analytic parabolic steady state solution. (mesh2a)

* test5/pipe-flow-old-1. Coordinate-aligned problem
* test6/pipe-flow-old-2. Rotated 45 degrees about x
* test7/pipe-flow-old-3. Rotated 45 degrees about y
* test8/pipe-flow-old-4. rotated 45 degrees about z

Meshes
------

### Mesh1a
* [-0.5,0.5]^3 domain
* Side set 1: x = -0.5 (pipe inlet)
* Side set 2: x = 0.5 (pipe outlet)
* Side set 3: remaining sides (pipe walls)

mesh1a.gen: Cartesian 5x5x5 hex mesh aligned with coordinate axes.
mesh1a-rotx.gen: mesh1a.gen rotated 45 degrees about the x axis
mesh1a-roty.gen: mesh1a.gen rotated 45 degrees about the y axis
mesh1a-rotz.gen: mesh1a.gen rotated 45 degrees about the z axis

### Mesh2a
* [-0.5,0.5] x [0,1] x [-0.3,0.3] domain
* Side set 1: x = -0.5 (pipe inlet)
* Side set 2: x = 0.5 (pipe outlet)
* Side set 3: y = 1 (pipe wall)
* Side set 4: z = constant, y = 0 sides (symmetry)

mesh2a.gen: Cartesian 5x10x3 hex mesh aligned with coordinate axes
mesh2a-rotx.gen: mesh2a.gen rotated 45 degrees about the x axis
mesh2a-roty.gen: mesh2a.gen rotated 45 degrees about the y axis
mesh2a-rotz.gen: mesh2a.gen rotated 45 degrees about the z axis
