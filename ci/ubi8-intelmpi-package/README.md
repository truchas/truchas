# Red Hat Universal Base Image 8.10 Intel MPI Binary Package Workflow

This local workflow builds a self-contained Truchas builder image for the
master branch.  The image uses Red Hat Universal Base Image (UBI) 8.10 and
contains oneAPI 2026.0, Intel MPI 2021.18.1, portable OpenBLAS, Python, and
`truchas-tpl` at `fb5ef5e38`.  All MPI-aware TPLs and Truchas are built with
the Intel MPI compiler wrappers.
Portage support is disabled, so the image does not build Wonton or Portage.
The workflow applies a local MUMPS patch to select the bundled OpenBLAS rather
than the MUMPS default Netlib LAPACK lookup.

```bash
./ci/ubi8-intelmpi-package/run.sh image
./ci/ubi8-intelmpi-package/run.sh package
./ci/ubi8-intelmpi-package/run.sh test
./ci/ubi8-intelmpi-package/run.sh smoke
```

`run.sh package` deliberately only builds and packages Truchas.  The source
CTest suite and tarball smoke test remain separate subcommands, matching the
current GitLab CI job separation.

The subcommands are intended to map one-to-one to future CI jobs:

- `image` creates the reusable builder image, including the compiler, MPI, and
  all TPLs.
- `package` builds the committed source snapshot and writes the tarball,
  version, compiler, and dependency provenance to `OUTPUT_DIR`.
- `test` builds the committed source snapshot and runs the serial CTest suite,
  writing `ctest-build.log` and `ctest.log` to `OUTPUT_DIR`.
- `smoke` uses a fresh UBI 8.10 container to test the package in `PACKAGE_DIR`.

The default builder-image tag is `truchas-ubi8-intelmpi-build:2026.0`.  Use
`IMAGE_TAG`, `JOBS`, `SOURCE_DIR`, `OUTPUT_DIR`, or `PACKAGE_DIR` to override
defaults.  By default, the workflow uses all processors visible to Podman.

## glibc portability

The builder image and its Red Hat UBI 8.10 base both use glibc 2.28.  The
binary package does not bundle glibc, so it requires a host with glibc 2.28 or
newer.  This makes RHEL 8.10 an intended target and is generally compatible
with newer RHEL-family systems, but it cannot run on RHEL 7-era systems with
glibc 2.17.  This is only the glibc compatibility floor; site-specific MPI,
kernel, and hardware requirements remain separate considerations.

## Intel MPI fabrics

The package bundles Intel MPI, Hydra, libfabric, and Intel's OFI provider
plugins.  The package launcher sets `I_MPI_FABRICS=shm:ofi` only when it has
not already been set; this is Intel MPI's regular-mode default.  It does not
select an OFI provider.  Sites may set `I_MPI_OFI_PROVIDER`, `FI_PROVIDER`, or
`FI_PROVIDER_PATH` to select a provider or locate site-supplied provider
plugins.  In particular, an InfiniBand verbs provider uses the host RDMA
libraries and kernel drivers.

This workflow is separate from the MPICH-based workstation package and is not
yet a claim of validated multinode or scheduler support.  Test it on each
target scheduler and interconnect before relying on it for cluster operation.

The build does not add `-static-intel`.  It packages the dynamically linked
Intel compiler runtimes instead; static linking should be evaluated separately
for its dependency and redistribution consequences.

## Verifying interconnect use

Run a real job spanning at least two nodes with Intel MPI startup diagnostics
enabled.  A one-node run primarily uses shared memory and cannot verify the
interconnect path.

```bash
export I_MPI_DEBUG=1
bin/mpiexec -n 2 bin/truchas case.inp 2>&1 | tee mpi-startup.log
```

The rank-zero startup output identifies the selected libfabric provider.  A
`verbs` (often `verbs;ofi_rxm`) provider uses the host RDMA verbs stack;
`mlx`, `psm3`, or `efa` identify their corresponding fabric stacks.  A
`tcp;ofi_rxm` provider uses TCP/IP rather than a high-performance RDMA path.

To print every provider visible to the package's libfabric installation, add:

```bash
export I_MPI_OFI_PROVIDER_DUMP=1
export I_MPI_DEBUG=1
bin/mpiexec -n 2 bin/truchas case.inp 2>&1 | tee mpi-provider-dump.log
```

For an InfiniBand confirmation independent of MPI diagnostics, a site
administrator can compare the selected HCA's `perfquery` or `ibqueryerrors`
transmit and receive counters before and after a multi-node run.
