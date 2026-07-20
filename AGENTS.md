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
Use an out-of-source build directory; in-source builds are rejected. Configure
the compiler, MPI wrappers, LAPACK, and required third-party libraries before
configuring Truchas. See `BUILDING.md` for prerequisites and platform-specific
instructions.

```sh
cmake -S . -B build -C config/linux-gcc.cmake \
  -D CMAKE_BUILD_TYPE=Debug \
  -D CMAKE_PREFIX_PATH=/absolute/path/to/truchas-tpl
cmake --build build --parallel
ctest --test-dir build --output-on-failure
cmake --install build
```

Example cache files for GCC, Intel, LLVM, and NAG are provided in `config/`.
Choose the cache file appropriate to the configured toolchain. Use
`ctest -R <name>` for a focused test and `ctest --parallel <jobs>` when parallel
execution is appropriate. For broad changes to shared infrastructure, run
representative builds or tests with multiple available compilers when practical.

## Coding Style & Naming Conventions
Most production code is modern Fortran using `.F90` files, lower-case module
and file names, and descriptive suffixes such as `_type`, `_class`, `_module`,
and `_factory`. Follow surrounding indentation and declaration layout in the
file being edited. Keep new CMake targets local to the relevant subdirectory
`CMakeLists.txt`. Python test scripts generally use `test*.py` names and should
remain small drivers around a Truchas input case.

## Testing Guidelines
Register regression tests in `test/CMakeLists.txt` with `add_pytruchas_test`,
including processor counts and labels such as `THERMAL`, `SPECIES`, `FLOW`, or
`ELECTROMAGNETICS`. Unit-style tests live beside their implementation in
`src/truchas/*/Test/` and are registered with CTest. When changing numerical
behavior, update or add the smallest representative test case and document any
golden-output changes in the commit or merge request. Most regression tests run
MPI programs, so a working MPI environment is required.

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
