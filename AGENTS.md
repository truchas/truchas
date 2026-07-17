# Repository Guidelines

## Project Structure & Module Organization
Truchas is a CMake-based multi-physics simulation code. Core source lives in
`src/truchas/`, grouped by domain such as `flow`, `solver`, `distributed_mesh`,
`physics`, `electromagnetics`, and `utilities`. Python support packages are in
`src/python/`, standalone tools are in `src/tools/`, and CMake helpers are in
`cmake/` and `config/`. Documentation is under `doc/`, with user-facing build
notes in `BUILDING.md` and `SPACK.md`. Regression tests and their input decks,
meshes, and golden outputs are organized by case in `test/`. See
`src/truchas/architecture.md` for a concise overview of the core source
architecture.

## Build, Test, and Development Commands
Use an out-of-source build directory; in-source builds are rejected.

```sh
ml gcc mpich truchas-tpl       # or: ml intel mpich truchas-tpl
mkdir build
cd build
cmake .. -C ../config/linux-gcc.cmake \
  -D CMAKE_BUILD_TYPE=Debug \
  -D CMAKE_PREFIX_PATH="$TRUCHAS_TPL_DIR"
make -j32
ctest -j32 --output-on-failure
make install
```

Compiler and MPI environments are controlled by lmod. Use
`ml gcc mpich truchas-tpl` for GCC, or replace `gcc` with `intel`, `llvm`, or
`nag` for the other supported toolchains. The `truchas-tpl` module defines
`TRUCHAS_TPL_DIR`; pass its value as `CMAKE_PREFIX_PATH` so CMake can find the
required TPLs. `cmake` configures compilers, MPI, LAPACK, and dependencies.
`make` builds libraries, tools, and executables. Use all 32 local cores when
practical, e.g. `make -j32` or `ctest -j32 --output-on-failure`. `ctest` runs
unit and regression tests; use `ctest -R <name>` for a focused test.

Use the NAG toolchain for most development work:
`ml nag mpich truchas-tpl`. Before pushing, also run representative builds or
tests with `gcc`, `intel`, and `llvm` to catch portability issues. NAG accesses
a license server, so any direct or indirect `nagfor` invocation requires
network access. This includes configure, compile, and link commands launched
through CMake or Make. Running a NAG-compiled executable does not by itself
require escalated privileges, but MPI-parallel execution does for any compiler.
Agents have standing permission to request the enhanced privileges needed for
NAG compilation/linking and MPI test commands.

## Coding Style & Naming Conventions
Most production code is modern Fortran using `.F90` files, lower-case module
and file names, and descriptive suffixes such as `_type`, `_class`, `_module`,
and `_factory`. Follow surrounding indentation and declaration layout in the
file being edited. Keep new CMake targets local to the relevant subdirectory
`CMakeLists.txt`. Python test scripts generally use `test*.py` names and should
remain small drivers around a Truchas input case.

## Testing Guidelines
Register regression tests in `test/CMakeLists.txt` with `add_pytruchas_test`,
including processor counts and labels such as `DIFFUSION`, `FLOW`, or
`ELECTROMAGNETICS`. Unit-style tests live beside their implementation in
`src/truchas/*/Test/` and are registered with CTest. When changing numerical
behavior, update or add the smallest representative test case and document any
golden-output changes in the commit or merge request. Tests use MPI, so agents
have permission to request escalated privileges when running `ctest` or focused
MPI test commands.

## Commit & Pull Request Guidelines
Recent history uses short, imperative summaries, sometimes with a scope or PR
number, for example `Port to LLVM clang/flang 22.1` or `Parallel VTKHDF EM
output`. Keep commits focused and explain behavioral or numerical changes in
the body when needed. Pull requests should describe the affected physics or
infrastructure area, list the build configuration used, include relevant
`ctest` results, and link associated issues.

## Configuration Notes
Do not vendor third-party libraries into this tree. Point CMake at installed
dependencies with `TRUCHAS_TPL_DIR`, `CMAKE_PREFIX_PATH`, or documented cache
files in `config/`. Keep generated build products, module files, and test output
outside the source tree.
