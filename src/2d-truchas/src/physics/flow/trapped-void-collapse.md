# Trapped-void pressure compliance prototype

Enable the model by adding a `void-collapse` sublist to `flow-model` (or the
corresponding flow-model sublist in coupled flow/thermal input):

```json
"void-collapse": {
  "pressure-time-scale": 1.0,
  "reaction-cap": 100.0
}
```

`pressure-time-scale` is the required positive product tau0*p0, with units
pressure times time in the simulation's unit system. The example value is
illustrative, not calibrated. Omitting the sublist disables collapse.
`reaction-cap` is a positive dimensionless cap, default 100. Geometric volume
tracking is required.

The model uses the existing inexpensive classification: a fluid/void mixed
cell with no completely void face neighbor. This is a local heuristic and
can include underresolved exterior-connected void. Pure-void cells retain
their existing pressure treatment. Cells containing solid are excluded from
this prototype. In eligible cells, alpha is the void fraction and
C = alpha / ((1-alpha)*tau0*p0). The nominal void pressure is zero gauge.

The linear, two-sided law is div(u_new) = -C*p_new. With p_new=p_old+phi,
the integrated pressure matrix receives r=V*C/dt on its diagonal and the
right-hand side receives -r*p_old. The reaction is capped at `reaction-cap`
times the local diffusion diagonal, and the same capped r is used in the
right-hand side. A cell with no diffusion stencil retains its reaction.
Positive reaction supplies a pressure reference, suppressing the existing
artificial pressure pin. At zero void the reaction vanishes.

No pressure clipping or active-set iteration is performed. A warning reports
the maximum measured positive divergence in reacting cells with negative
updated pressure, above a roundoff threshold. This warning does not change
the solution.

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
pressure loading, reaction capping and pressure pinning, and a pressure-fed
stationary mixed layer against a closed wall. The coupled test checks liquid
gain against the boundary inflow at each transport step on one and four MPI
ranks. Timestep, mesh, and cap sensitivity remain prototype validation work;
the regression is not a calibration of the model.
