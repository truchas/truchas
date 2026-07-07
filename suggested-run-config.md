# Suggested Run Configuration

If one rank dies before or during an MPI collective, the remaining ranks may sit forever in MPI calls, making the failure look like a hang. The real failure message may have been printed only by the dying rank, or lost in mixed or buffered output.

The diagnostic goal is to make rank failure more obvious and preserve rank-specific output.

For MPICH, useful environment variables are typically:

```bash
export MPICH_ABORT_ON_ERROR=1
export MPICH_ERROR_CHECKING=1
```

Then run with rank-labeled output if the launcher supports it:

```bash
mpiexec -n 64 -l ./truchas input.inp > truchas.out 2> truchas.err
```

The `-l` option labels each output line with the rank. If rank 0 is the likely problem, a stronger version is to use a site wrapper or launcher option that splits output by rank, so rank 0 stderr is not interleaved with other ranks.

The point is to convert an apparent hang into the first real failure: an MPI error, Fortran runtime error, allocation failure, assertion, backtrace, or OS OOM-kill indication.
