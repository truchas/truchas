# Red Hat Universal Base Image 8.10 oneAPI Binary Package Workflow

This local workflow builds a self-contained Truchas builder image for the
master branch.  The image uses Red Hat Universal Base Image (UBI) 8.10 and
contains oneAPI 2026.0, MPICH 5.0.1, portable OpenBLAS, Python, and
`truchas-tpl` at `fb5ef5e38`.
Portage support is disabled, so the image does not build Wonton or Portage.
The workflow applies a local MUMPS patch to select the bundled OpenBLAS rather
than the MUMPS default Netlib LAPACK lookup.

```bash
./ci/ubi8-oneapi-package/run.sh image
./ci/ubi8-oneapi-package/run.sh package
./ci/ubi8-oneapi-package/run.sh test
./ci/ubi8-oneapi-package/run.sh smoke
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

The default builder-image tag is `truchas-ubi8-oneapi-build:2026.0`.  Use
`IMAGE_TAG`, `JOBS`, `SOURCE_DIR`, `OUTPUT_DIR`, or `PACKAGE_DIR` to override
defaults.  By default, the workflow uses all processors visible to Podman.

## glibc portability

The builder image and its Red Hat UBI 8.10 base both use glibc 2.28.  The
binary package does not bundle glibc, so it requires a host with glibc 2.28 or
newer.  This makes RHEL 8.10 an intended target and is generally compatible
with newer RHEL-family systems, but it cannot run on RHEL 7-era systems with
glibc 2.17.  This is only the glibc compatibility floor; site-specific MPI,
kernel, and hardware requirements remain separate considerations.
