# Trapped-void pressure compliance prototype

Enable the model by adding a `void-collapse` sublist to `flow-model` (or the
corresponding flow-model sublist in coupled flow/thermal input):

```json
"void-collapse": {
  "pressure-time-scale": 1.0,
  "reaction-cap": 100.0,
  "capillary-coefficient": 0.0
}
```

`pressure-time-scale` is the required positive product tau0*p0, with units
pressure times time in the simulation's unit system. The example value is
illustrative, not calibrated. Omitting the sublist disables collapse.
`reaction-cap` is a positive dimensionless cap, default 100. Geometric volume
tracking is required.

`capillary-coefficient` is the optional nonnegative effective surface tension
Sigma = c_sigma*sigma, with units pressure times length. Its default is zero,
which retains unbiased compliance. A positive value adds the capillary bias
from section 6.2 of the rev3 research note:

```
Pi = Sigma*(1-alpha)^3/sqrt(cell_area)
div(u_new) = -C*(p_new + Pi)
```

The equilibrium liquid pressure is -Pi: pressure above this value causes
collapse; pressure below it permits expansion. The cubic shape and
cell length sqrt(cell_area) are fixed. Both C and Pi use the current material
fractions, after transport during stepping and from the initial distribution
during the initial physical-pressure solve. Eligibility is unchanged.

The model uses the existing inexpensive classification: a fluid/void mixed
cell with no face neighbor classified as `void_t`. Such a neighbor contains
no real fluid, but can contain solid if its nonsolid fraction is at least the
solid-cell cutoff. This is a local heuristic and can include underresolved
exterior-connected void. Cells classified as VOID retain their existing
pressure treatment.

**Limitation:** collapse is disabled in cells containing solid, including
three-component solid/liquid/VOID cells. This is a prototype restriction,
not a fundamental exclusion from pressure compliance. In eligible cells,
alpha is the void fraction and C = alpha / ((1-alpha)*tau0*p0). The nominal
void pressure is zero gauge. A proposed extension normalizes liquid and VOID
fractions over the nonsolid portion and uses its area to define the cell
length. It remains deferred pending rigorous testing, especially for evolving
solid and near-cutoff nonsolid regions; see the
[solid/liquid/VOID extension TODO](../../../TODO.md#volume-tracking-and-interfaces).

The linear, two-sided law is div(u_new) = -C*p_new. With p_new=p_old+phi,
the integrated pressure matrix receives r=V*C/dt on its diagonal and the
right-hand side receives -r*(p_old+Pi). The reaction is capped at `reaction-cap`
times the local diffusion diagonal, and the same capped r is used in the
right-hand side, including the bias contribution. The bias changes no matrix
coefficients and adds no momentum force or separate volume-fraction update.
A cell with no diffusion stencil retains its reaction.
Positive reaction supplies a pressure reference, suppressing the existing
artificial pressure pin. At zero void the reaction vanishes.

No pressure clipping or active-set iteration is performed. A warning reports
the maximum measured positive divergence in reacting cells with updated
pressure below the shifted equilibrium, above a roundoff threshold. Negative
pressure alone no longer triggers the warning. This diagnostic is retained
for prototype testing and does not change the solution.

The existing transport-first timestep is retained: geometric transport uses
the previously committed face velocity; the subsequent pressure solve uses
the transported fractions and supplies the next face velocity. Collapse is
also included in the initial physical pressure solve, but not in the initial
kinematic velocity projection. Initialization retains the prescribed, projected
velocity and discards the artificial pressure step's velocity update. Thus a
stationary initial condition starts replacing void on the second transport
step. The geometric tracker's existing overfill
handling removes void. No additional direct edit of the volume fractions
is introduced. For the uncapped continuum law, fixed positive pressure gives
exponential decay; the cap and transport cutoff modify endpoint behavior.

Regression coverage includes the incremental pressure offset, both signs of
pressure loading relative to the capillary equilibrium and equilibrium itself,
reaction capping and pressure pinning, and a pressure-fed
stationary mixed layer against a closed wall. The coupled test checks liquid
gain against the boundary inflow at each transport step on one and four MPI
ranks. Timestep, mesh, and cap sensitivity remain prototype validation work;
the regression is not a calibration of the model.
