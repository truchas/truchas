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
