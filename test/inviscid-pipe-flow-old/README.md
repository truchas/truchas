INVISCID PIPE FLOW TESTS
========================

## Tests

Inviscid flow in pipe with square cross section. Free-slip BC on pipe walls.
Fluid accelerates from rest (constant acceleration). Check uniform velocity
value after a couple steps and at final time. (mesh1a)

* test1a/pipe-flow-inviscid-old-1a. Coordinate-aligned problem
* test1b/pipe-flow-inviscid-old-1b. Rotated 45 degrees about x
* test1c/pipe-flow-inviscid-old-1c. Rotated 45 degrees about y
* test1d/pipe-flow-inviscid-old-1d. rotated 45 degrees about z

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
