# Two-dimensional Truchas

The `truchas-2d` executable is the common driver for the two-dimensional
simulations.  Select the simulation with `-s` or `--simulation`; the input
file is a positional JSON argument.

```text
mpiexec -n 4 truchas-2d -s flow input.json
```

The supported simulations are:

| Name | Simulation |
| --- | --- |
| `flow` | Isothermal incompressible flow |
| `thermal` | Thermal transport |
| `flow_thermal` | Non-isothermal incompressible flow with thermal transport |

Use `truchas-2d --help` for the complete command-line help.  The common
options are:

```text
-s NAME, --simulation NAME   Select the simulation.
-o DIR,  --output-dir DIR    Write run products in DIR.
-f,      --force             Permit use of an existing output directory.
-h,      --help              Display help and exit.
```

If `--output-dir` is omitted, the output directory is the stem of the input
file.  The run log and simulation output are written there.  MPI initialization
and command-line diagnostics are handled by the common driver; the selected
simulation owns interpretation of the JSON input and the simulation lifecycle.
