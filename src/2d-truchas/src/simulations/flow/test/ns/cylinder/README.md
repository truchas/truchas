# Cylinder-flow test

This four-process reference regression computes laminar flow at Reynolds
number 20 past a solid cylinder.  It uses the upper half of a channel, with a
free-slip symmetry boundary through the cylinder center, a velocity inlet, a
pressure outlet, and a free-slip outer wall.  The inlet velocity is ramped
from zero to its final value over the first time unit, and the mesh is graded
from approximately 2:1 cells at the inlet through 1:1 cells around the
cylinder to 4:1 cells over the right half of the domain.

There is no analytic solution.  The test compares the velocity, pressure, and
material volume fractions at `t=0`, `0.25`, `0.5`, `1`, and `2` with the
stored VTKHDF reference.  It also checks that the final velocity field contains a
nontrivial wake and that the cylinder and fluid fractions sum to one.
