# Lid-driven cavity examples

`lid_driven_cavity.json` is a developer example for the two-dimensional,
isothermal Stokes solver. `lid_driven_cavity_ns.json` is the corresponding
time-dependent Navier--Stokes example. The former is intentionally not
registered with CTest; the latter is exercised by the integrated NS test.

With a configured build, run it from a scratch directory so that the
`out.vtkhdf` file does not overwrite another run:

```sh
ml purge
ml nag mpich truchas-tpl
mkdir cavity-run
cd cavity-run
path/to/build/src/2d-truchas/flow/ns_2d \
  path/to/source/src/2d-truchas/flow/Test/lid_driven_cavity_ns.json
```

Both inputs start from zero velocity, impose unit tangential velocity on the
top wall, and impose no-slip on the other walls. The output contains the
time-dependent cell pressure and three-component VTK velocity, with the third
component zero. The three-component representation is required by ParaView's
streamline filter even though the physical problem is two-dimensional.

## Time control

Both 2D flow programs use the same `sim-control` time parameters as the
thermal simulation: `initial-time`, `initial-time-step`, `min-time-step`,
`max-time-step`, and a strictly increasing `output-times` array. The initial
state is written at `initial-time`; integration stops at the final entry of
`output-times`. `time-step-growth` defaults to `1.05` and limits step growth.
For Navier--Stokes, `courant-number` additionally limits each step by the
convective CFL condition.

`inviscid_channel.json` is a small developer/CTest case for the inviscid
Navier--Stokes path. It omits viscosity and uses free-slip walls.
