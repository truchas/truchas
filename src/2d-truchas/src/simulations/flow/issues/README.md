# Known flow issues

This directory preserves flow configurations that expose known numerical
weaknesses. They are diagnostic cases, not active CTest regressions. Each input
file contains a short description of its physical setup; this file records the
reason for preserving the case and the observed behavior.

## Noisy hydrostatic free surface

These cases use a unit-square water/VOID configuration with downward gravity,
a 16x16 mesh, a 10% seeded mesh perturbation, and geometric material-region
initialization at refinement level 6. The intended solution is static, with
zero fluid velocity and hydrostatic pressure.

`hydrostatic_interface_aligned.json` places the nominal interface at `y=0.5`,
on a regular-mesh cell plane. Perturbing the mesh then creates very small water
and VOID fragments near the interface. With the current face-based VOID
projection path, this case develops substantial spurious motion: the observed
maximum fluid speed at `t=1` is approximately `1.4e-1`.

`hydrostatic_interface_cell_centered.json` moves the interface by half a cell,
to `y=0.5+1/32`. This avoids the near-zero fragments on the nominal mesh and
reduces the observed maximum speed by roughly a factor of six. The comparison
isolates sensitivity to small material fractions from general mesh
non-orthogonality.

The corresponding mainline 3D experiment shows qualitatively similar
behavior. These cases should become useful regression tests after the
free-surface projection, cutoff, and normalization treatment is improved.

## Solid/fluid no-slip discretization

`poiseuille_solid_wall.json` is a viscous pressure-driven channel whose
physical walls are represented by pure SOLID cells.  The companion checker,
`test_poiseuille_solid_wall.py`, compares the result with the analytic
Poiseuille profile.  It is intentionally not a CTest regression: the current
solid/fluid viscous discretization places the effective no-slip wall at an
inaccurate location, even though the material interface is mesh-aligned and
there are no mixed cells.  The checker is retained as a manual diagnostic and
is expected to report the current velocity-profile error until the interface
discretization is improved.

## Pressure-fed closed column

`trapped_void_column.json` is a 1-by-10 column on a 1-by-10 mesh. The bottom
two cells initially contain inviscid, unit-density water; the remaining eight
contain VOID. The bottom has unit pressure and water inflow composition;
the sides and top are free-slip. Gravity and initial velocity are zero.
Geometric tracking and two-sided trapped-void compliance are enabled, with
pressure-time-scale 1 and reaction-cap 100.

For the ideal sharp-interface column before impact, the liquid height obeys
h''=1/h with h(0)=2 and h'(0)=0. Integrating to h=10 gives an estimated filling
time of 6.948. The input runs to t=30, outputs every 0.25, and limits the
timestep to 0.01. Timestep growth allows recovery after a step is shortened
to land on an output time.

A one-rank NAG run reached t=30 in 3072 accepted steps. The first saved fully
filled state was t=7; all ten cells remained full thereafter. At t=30 the top
pressure was 0.99999775 and the maximum cell vertical speed was 2.20e-7.
No positive-divergence warning was emitted. This coarse column filled without
a persistent residual void in that run; it is a diagnostic, not a demonstration
that collapse is required for this configuration.

Run with the flow simulation, for example:

```sh
mpiexec -n 1 /path/to/truchas-2d --simulation flow --output-dir column-output trapped_void_column.json
```
