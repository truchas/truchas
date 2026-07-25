# RHEL 8.10 Binary Package Build

This is a local, rootless-Podman workflow for creating a `24.04-wnaa` binary
package that targets the RHEL 8.10 ABI. It does not alter GitLab CI or any
existing CI recipe.

The workflow preserves the original 24.04 compiler family: Intel Parallel
Studio XE 2020 Update 4, including `ifort`, `icc`, and `icpc` 19.1.3.304.
Portage is disabled in both the TPL and Truchas builds, avoiding the legacy
Boost dependency.

## Inputs

- `registry.access.redhat.com/ubi8/ubi:8.10`
- Intel compiler archive at
  `/archive/Library/Archive/Intel/parallel_studio_xe_2020_update4_composer_edition.tgz`
- Intel license files in `/opt/intel/licenses`
- MPICH 3.3.2
- Intel MKL bundled with Parallel Studio XE, selected through
  `BLA_VENDOR=Intel10_64lp`
- CMake 3.22.2
- `truchas-tpl` v22
- Python CMake build-system commit `11a369e29b1be82e3b3108259ab89b7095854b47`
- EPEL 8 package `patchelf` (the only package taken from EPEL)

## Build

Run all commands as your normal user:

```bash
cd ci/rhel8-ubi-package
./build-image.sh
./build-package.sh
```

`build-package.sh` builds the dependencies once in the named volume
`truchas-rhel8-24.04-mkl-dependencies`. Completed MPICH, CMake, Python, and
TPL stages have independent state markers, so a later failure resumes from the
first incomplete stage. The output directory defaults to `rhel8-package-output`
under the directory where the command is run. The dependency log is retained in
the volume as `dependency-build.log`.

Set `OUTPUT_DIR` to choose a different output directory. Set
`DEPENDENCY_VOLUME` to create a clean dependency build for a new iteration.
Builds use 32 parallel jobs by default; override this with `JOBS`.

## Validation

Before delivery, test the output tarball in a fresh UBI 8.10 container and on
the user's RHEL 8.10 system. Confirm that the package does not require a newer
glibc version:

```bash
readelf --version-info path/to/truchas | rg 'GLIBC_|GLIBCXX_'
```

One `.lic` file is mounted only while the image is built. The license directory
is then mounted read-only for compiler invocations; neither the license file
nor its directory is copied into the resulting image or package.

The build processes run as root inside a rootless user namespace. That root
identity maps to the invoking host user; it does not grant host root access.
