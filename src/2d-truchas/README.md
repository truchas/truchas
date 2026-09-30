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

When the pressure-compliance void-collapse model is enabled, the `flow`
VTKHDF output includes the temporary cell diagnostic
`void_compliance`. It records the pressure-compliance coefficient
`C = (void fraction / real-fluid fraction) / impedance` in eligible cells,
in inverse pressure-time units, with NaN in inactive cells. This coefficient
is omitted when compliance is disabled or the mainline collapse model is
selected. Values use the material state at each output time, including
initialization, and precede
the projection's numerical reaction cap. Positive values identify cells
where compliance is enabled, even when the pressure is at equilibrium and
the resulting divergence is zero.

The companion `void_target_divergence` diagnostic records
`-C_eff * (pressure + pressure_bias)` from the accepted pressure correction,
using the capped compliance `C_eff`. Its units are inverse time: negative
values indicate contraction, positive values expansion, and zero equilibrium.
Cells with no compliance reaction have NaN values. The initial output uses
the temporary Stokes pressure solve's target; initialization retains the
prescribed, repaired velocity, so this target need not equal its divergence.

For the mainline model, `void_target_divergence` instead records `-q/dt`,
where `q = min(abs(alpha - alpha_previous), relaxation * alpha)` is the
requested collapse fraction. It is saved during the pressure correction,
before accepting the step resets the material history. Eligible cells have
negative values or zero when no collapse is requested; ineligible cells have
NaN. The initial mainline target is zero in eligible cells.
`void_target_divergence` is omitted when both collapse models are disabled.
Both diagnostic fields include synchronized off-process cell values.

The root-level `initial-time`, `output-times`, and optional `final-time`
define the simulation schedule.  `output-times` is a sublist with a `times`
array and an optional
`subintervals` integer array.  Each entry in `subintervals` divides the
interval to the next entry in `times`; output is written at the resulting
endpoints.  Its length is one less than the length of `times`.  If it is
omitted, each interval is divided once.  `initial-time` remains the starting time
and is output automatically.  An optional `final-time` extends the schedule
with a final output and sets the integration endpoint; without it, the last
scheduled time is the integration endpoint.  For example:

```json
"initial-time": 0.0,
"output-times": {
  "times": [0.0, 1.0, 2.0],
  "subintervals": [2, 4]
},
"final-time": 3.0
```
