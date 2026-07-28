# Truchas Binary Package for RHEL 8 and Newer

This x86_64 package was built with Intel oneAPI 2026.0, Intel MPI 2021.18.1,
OpenBLAS, and Red Hat Universal Base Image (UBI) 8.10. It includes Truchas,
its third-party libraries, Intel MPI, Hydra, libfabric, and the bundled Python
utilities.

The package does not include glibc. It requires glibc 2.28 or newer, so RHEL
8.10 and newer RHEL-family systems are intended targets. RHEL 7 is not
supported.

## Running Truchas

Use `bin/truchas` rather than invoking `bin/t-linux.x86_64.gnu` directly. The
launcher sets an unlimited stack size, which is required by some Intel Fortran
workloads, and initializes the package's libfabric provider search path.

An included Broken Dam input provides a quick serial check:

```bash
cd examples/broken-dam
../../bin/truchas broken-dam.inp
```

To run the same problem with four local MPI ranks:

```bash
cd examples/broken-dam
../../bin/mpiexec -n 4 ../../bin/truchas broken-dam.inp
```

`bin/python` and `bin/write-restart.py` provide the packaged Python runtime
and restart-file utility.

## MPI and interconnects

The package bundles Intel MPI, Hydra, libfabric, and Intel's OFI provider
plugins. The launcher uses Intel MPI's normal `shm:ofi` fabric choice unless
`I_MPI_FABRICS` is already set. It does not force a provider selection.

For a multi-node run, set `I_MPI_DEBUG=1` and inspect the rank-zero startup
output:

```bash
export I_MPI_DEBUG=1
bin/mpiexec -n 2 bin/truchas case.inp 2>&1 | tee mpi-startup.log
```

Run this across at least two nodes. `verbs` (often `verbs;ofi_rxm`) indicates
the host RDMA verbs stack; `mlx`, `psm3`, and `efa` identify their respective
fabric stacks. `tcp;ofi_rxm` means MPI is using TCP/IP rather than an RDMA
transport. Set `I_MPI_OFI_PROVIDER_DUMP=1` to print all providers visible to
the package.

The MPI runtime has been tested for local execution. Multinode and scheduler
integration must be tested by each site with its scheduler and interconnect.
