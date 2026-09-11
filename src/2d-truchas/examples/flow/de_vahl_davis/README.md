# de Vahl Davis natural-convection benchmark

This example is the classic natural-convection problem in a differentially
heated square cavity, associated with the de Vahl Davis benchmark. The left
wall is held hot, the right wall cold, the horizontal walls are insulated, and
all four walls are no-slip. The fluid starts from rest with the linear
conductive temperature profile.

This is the classic de Vahl Davis (1982) benchmark problem.

The problem is nondimensional with:

- `Pr = 0.71`;
- `Ra = 10000`;
- unit cavity width and temperature difference.

The input achieves these values with density, conductivity, and specific heat
equal to one, viscosity `0.71`, thermal expansion coefficient one, and a
downward body acceleration of `7100`. It uses a 129-by-129 quadrilateral mesh
and integrates from `t=0` to `t=0.3`, writing output at `t=0.1`, `0.2`, and
`0.3`. The `t=0.2` and `t=0.3` fields are visually indistinguishable, so the
final output is a useful near-steady benchmark state, although the run is not
intended to establish formal convergence.

The benchmark is commonly used to compare centerline velocity and temperature
profiles, extrema, and heat-transfer measures against published results. This
is a longer-running flow/thermal example and is not registered as a CTest
regression.

Run it with, for example:

```text
mpiexec -n 16 truchas-2d --simulation flow_thermal \
  --output-dir de-vahl-davis-output --force input.json
```

The VTKHDF and log files produced by a run are artifacts and are not part of
the example source.
