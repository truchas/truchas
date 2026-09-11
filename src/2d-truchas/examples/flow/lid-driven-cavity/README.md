# Ghia–Ghia–Shin lid-driven cavity benchmark

This example is the standard nondimensional lid-driven cavity problem at
Reynolds number 1000, following Ghia, Ghia, and Shin, “High-Re Solutions for
Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method,”
*Journal of Computational Physics* 48 (1982), 387–411.

The domain is the unit square. The fluid is initially at rest, the upper lid
moves with velocity `(1, 0)`, and the other three walls are no-slip. Density
and the characteristic length and speed are one; viscosity is `1e-3`, giving
`Re = 1000`. The input uses a 129-by-129 quadrilateral mesh and integrates to
nondimensional time 40, with output at five-unit intervals.

The usual benchmark comparison is against the horizontal and vertical
centerline velocity profiles tabulated by Ghia et al. This is a longer-running
example for visual and quantitative benchmarking, not a CTest regression.

Run it with, for example:

```text
mpiexec -n 16 truchas-2d --simulation flow \
  --output-dir lid-driven-cavity-output --force input.json
```

The VTKHDF and log files produced by a run are artifacts and are not part of
the example source.
