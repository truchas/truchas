# Martin–Moyce dam-break benchmark

This is a nondimensional two-dimensional free-surface dam-break problem based
on the classic Martin–Moyce benchmark. A viscous fluid initially occupies a
1-by-2 column in the left side of a 5-by-3 domain; the remaining domain is
VOID. Gravity acts in the negative y direction. The walls are free-slip and
the top boundary has prescribed zero pressure.

The case uses a 250-by-150 mesh and integrates to nondimensional time 3.5.
It is intended as a longer-running example and visual benchmark, not as a
routine CTest test. A representative four-process run takes about one minute
in an optimized build and accepts 1463 time steps.

Run it with, for example:

```text
mpiexec -n 4 truchas-2d --simulation flow \
  --output-dir dam-break-output --force input.json
```

The output schedule is specified in `input.json`. The generated VTKHDF and
log files are run artifacts and are not part of the example source.
