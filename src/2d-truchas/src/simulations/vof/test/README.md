# 2D volume-tracking tests

- `advection`: translates a circular region on a regular quadrilateral mesh with a constant velocity.
- `vortex`: advects a circular region in a time-dependent, divergence-free vortex field.
- `axisymmetric/regular`: advects a circular region in an axisymmetric `r-z` domain on a quadrilateral mesh.
- `axisymmetric/mixed`: the same axisymmetric problem on a mesh containing triangles and quadrilaterals.
- `axial_flow`: axisymmetric advection in a radial/axial velocity field on a perturbed mesh.
- `generic`: exercises the generic VOF simulation with a stationary circular region.
- `truncated_volumes`: directly checks geometric truncation volumes for planar and axisymmetric cells.

The VOF cases run the generic `truchas-2d --simulation vof` driver and compare partitioned VTKHDF output with serial gold data. `truncated_volumes` is a direct executable test of the geometric utility.
