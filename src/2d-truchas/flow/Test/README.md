# Standalone lid-driven cavity

`lid_driven_cavity.json` is a developer example for the two-dimensional,
isothermal Stokes solver. It is intentionally not registered with CTest.

With a configured build, run it from a scratch directory so that the
`out.vtkhdf` file does not overwrite another run:

```sh
ml purge
ml nag mpich truchas-tpl
mkdir cavity-run
cd cavity-run
path/to/build/src/2d-truchas/flow/flow_2d \
  path/to/source/src/2d-truchas/flow/Test/lid_driven_cavity.json
```

The input starts from zero velocity, imposes unit tangential velocity on the
top wall, and imposes no-slip on the other walls. The output contains the
time-dependent cell pressure and three-component VTK velocity, with the third
component zero. The three-component representation is required by ParaView's
streamline filter even though the physical problem is two-dimensional.
