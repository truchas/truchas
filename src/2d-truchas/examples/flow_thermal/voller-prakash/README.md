# Voller–Prakash solidification benchmark

This example is one of the classic natural-convection-driven solidification
problems introduced by Voller and Prakash (1987). It is a differentially
heated unit-square cavity containing a single material with solid and liquid
phases.

The reference is V. R. Voller and C. Prakash, “A fixed grid numerical
modelling methodology for convection-diffusion mushy region phase-change
problems,” *International Journal of Heat and Mass Transfer*, 30(8),
1709–1719 (1987), [doi:10.1016/0017-9310(87)90317-6](https://doi.org/10.1016/0017-9310(87)90317-6).

The material is initially liquid and isothermal at `T=0.5`. The left wall is
cooled to `T=-0.5`, the right wall is held at `T=0.5`, and the top and bottom
walls are adiabatic. All flow walls are no-slip, and the body acceleration is
downward. Solidification begins at the cold wall and interacts with the
buoyancy-driven flow.

The reference nondimensional parameters are:

- `Ra = 10^4`;
- `Pr = 1000`;
- `Gr = 10`;
- `Ste = 5` under the modern convention.

The input uses a 64-by-64 quadrilateral mesh and backward-Euler integration.
Output is written at the listed times through `t=1000`; the initial state is
also written automatically.

## More vigorous-flow variant

For a more vigorous flow, reduce the material viscosity from `1` to `0.1`
and increase its conductivity from `1e-3` to `1e-2`. The resulting case has
`Pr = 10` and `Gr = 1000`, with the other input parameters unchanged.

Run the reference case from the 2D source tree with, for example:

```sh
mpiexec -n 4 truchas-2d --simulation flow_thermal \
  --output-dir voller-prakash-output --force \
  examples/flow_thermal/voller-prakash/input.json
```
