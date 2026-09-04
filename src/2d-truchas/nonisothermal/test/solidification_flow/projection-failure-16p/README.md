# 16-process projection-failure case

This is a retained investigation artifact, not currently a CTest regression.
It exercises buoyancy-driven flow coupled to solidification and reproduces a
projection-solver failure after dynamically changing solid cells appear.
The test problem is a version of the Voller--Prakash (1987) benchmark problem.

The captured run artifacts are in this directory:

- `input.json` is the original input.
- `output-input.json` is the input copy written by the program.
- `run.log` is the 16-process NAG release log.
- `out.vtkhdf` is the corresponding output file.

Run provenance from `run.log`:

- compiler: NAG 7.2.7245, release flags, `-O3 -DNDEBUG`;
- MPI processes: 16;
- input SHA-256: `2511c1ace294a6b2d151eb4841822f4e2ed6c61ab7a14e8e93c18258aca11c47`;
- failure: projection solve at step 1687, at `t=415.228`;
- final diagnostic: `iter=43 dscg=31 amg=12 rel_res=3.55425E-07`.

The accompanying instrumented investigation found fresh ghost values and
matrix-action antisymmetry of approximately `5E-14` relative to a matrix scale
of approximately `2E+02`.  Thus the artifact is useful for studying the
remaining projection-solver failure; it should not be treated as evidence of
a demonstrably nonsymmetric assembled matrix.
