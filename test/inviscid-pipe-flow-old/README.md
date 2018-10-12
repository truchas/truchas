INVISCID PIPE FLOW TESTS
========================

## Tests

Inviscid flow in pipe with square cross section. Free-slip BC on pipe walls.
Fluid accelerates from rest (constant acceleration). Check uniform velocity
value after a couple steps and at final time. (mesh1a)

* test1a/inviscid-pipe-flow-old-1a. Coordinate-aligned problem
* test1b/inviscid-pipe-flow-old-1b. Rotated 45 degrees about x
* test1c/inviscid-pipe-flow-old-1c. Rotated 45 degrees about y
* test1d/inviscid-pipe-flow-old-1d. rotated 45 degrees about z

Same problem but the domain includes a 1-cell layer of solid material around
the wall sides, which replace the explicit free-slip BC there. (mesh2a)

* test2a/inviscid-pipe-flow-old-2a. Coordinate-aligned problem
* test2b/inviscid-pipe-flow-old-2b. Rotated 45 degrees about x
* test2c/inviscid-pipe-flow-old-2c. Rotated 45 degrees about y
* test2d/inviscid-pipe-flow-old-2d. rotated 45 degrees about z

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

### Mesh2sa
* [-0.5,0.5]x[-0.7,0.7]^2 domain
* Block 1 is fluid: [-0.5,0.5]^2
* Block 2 is solid: complement
* Side set 1: x = -0.5 (pipe inlet)
* Side set 2: x = 0.5 (pipe outlet)

mesh2a.gen: Cartesian 5x5x5 hex mesh aligned with coordinate axes.
mesh2a-rotx.gen: mesh2a.gen rotated 45 degrees about the x axis
mesh2a-roty.gen: mesh2a.gen rotated 45 degrees about the y axis
mesh2a-rotz.gen: mesh2a.gen rotated 45 degrees about the z axis
