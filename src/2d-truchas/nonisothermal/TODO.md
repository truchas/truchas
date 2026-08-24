# Non-isothermal flow TODOs

## Improve hydrostatic well-balance

The initial-condition solver now uses the same temperature-dependent
gravity-head path as the normal pressure correction. The hot-top/cold-bottom
test shows systematic spatial convergence, but still has a discretization
residual at finite resolution. A future operator refinement should reduce
that residual while preserving the current initial-pressure construction.
