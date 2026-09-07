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
