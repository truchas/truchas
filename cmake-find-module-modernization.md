# CMake Find Module Modernization Proposal

This project carries private `Find*.cmake` modules for third-party libraries
that do not reliably provide official CMake find modules or package config
files. The goal is to modernize those modules while preserving compatibility
with the existing build.

## Goals

- Prefer target-based usage through imported targets.
- Preserve legacy result variables during the transition.
- Keep support for existing TPL layouts and environment hints.
- Make `REQUIRED`, `QUIET`, component, and version handling consistent.
- Avoid broad call-site churn until each module has compatibility aliases.

## Proposed Standard

Each private find module should expose a namespaced imported target as the
primary result:

```cmake
HYPRE::HYPRE
NetCDF::C
NetCDF::Fortran
Exodus::Exodus
PETACA::PETACA
Chaparral::Chaparral
Scorpio::Scorpio
METIS::METIS
YAJL::YAJL
MUMPS::DMUMPS
MUMPS::ZMUMPS
```

Existing unnamespaced targets such as `hypre`, `netcdf`, `exodus`, `petaca`,
`chaparral`, `scorpio`, and `metis` should remain as compatibility aliases
until all call sites have migrated.

Each module should continue setting traditional variables such as
`FOO_FOUND`, `FOO_INCLUDE_DIRS`, `FOO_LIBRARIES`, and `FOO_VERSION` where they
are already part of the local interface. These variables should be treated as
compatibility output; new CMake code should link to targets.

## Transitive Dependencies

Transitive usage requirements should be represented by target links rather than
by flattening dependency include directories and libraries into variables.

Examples:

```cmake
target_link_libraries(Exodus::Exodus INTERFACE NetCDF::C)
target_link_libraries(HYPRE::HYPRE INTERFACE MPI::MPI_C)
target_link_libraries(PETACA::PETACA INTERFACE YAJL::YAJL)
```

This makes include directories, compile definitions, and link requirements flow
through CMake's normal target usage model.

## Search Behavior

Search behavior should be normalized across modules:

- Prefer standard CMake discovery through `CMAKE_PREFIX_PATH`.
- Support `<PackageName>_ROOT` as the user-facing package root hint.
- Preserve existing environment hints such as `HYPRE_ROOT`, `YAJL_ROOT`, and
  `MUMPS_ROOT`.
- Avoid replacing `CMAKE_PREFIX_PATH` globally. If `TRUCHAS_TPL_DIR` remains,
  append or prepend it instead of overwriting user-provided paths.

When a package may provide a modern config package, the private module should
try config mode quietly first:

```cmake
find_package(Foo CONFIG QUIET)
if(Foo_FOUND)
  # Normalize repo-standard targets and compatibility variables here.
  return()
endif()
```

The important property is that Truchas call sites should not care whether the
dependency came from an upstream config package or from the fallback search.

## Module Cleanup Targets

Recommended cleanup items by module:

- `FindYAJL.cmake`: keep `YAJL::YAJL`, normalize root handling, and preserve
  legacy variables.
- `FindChaparral.cmake`: add `Chaparral::Chaparral` and keep `chaparral` as an
  alias.
- `FindScorpio.cmake`: add `Scorpio::Scorpio`, use MPI targets for parallel
  requirements, and keep `scorpio` as an alias.
- `FindMETIS.cmake`: compute and honor `METIS_VERSION`, keep component handling
  for ParMETIS, and remove top-level recreation of the `metis` target.
- `FindPETACA.cmake`: link to `YAJL::YAJL` instead of raw YAJL variables and add
  `PETACA::PETACA`.
- `FindExodus.cmake`: expose `Exodus::Exodus` and link it to `NetCDF::C`.
- `FindHYPRE.cmake`: expose `HYPRE::HYPRE`, link parallel builds to
  `MPI::MPI_C`, and keep the `HYPRE_IS_PARALLEL` suitability check.
- `FindNetCDF.cmake`: expose `NetCDF::C` and `NetCDF::Fortran`, preserve NC4
  and large-model checks, and use `nc-config`/`nf-config` as hints where useful.
- `FindMUMPS.cmake`: normalize config-mode and fallback behavior, fix imported
  target properties that currently reference unset variables, and model
  dependencies such as ScaLAPACK/LAPACK/BLAS explicitly.

## Consumer Migration

The central Truchas link line currently mixes target styles. It should migrate
from unnamespaced compatibility targets to namespaced imported targets.

Current style:

```cmake
LAPACK::LAPACK exodus hypre petaca hdf5::hdf5 fVTKHDF::fvtkhdf
METIS::METIS scorpio MPI::MPI_Fortran netcdf
```

Preferred style:

```cmake
LAPACK::LAPACK
Exodus::Exodus
HYPRE::HYPRE
PETACA::PETACA
HDF5::HDF5
fVTKHDF::fvtkhdf
METIS::METIS
Scorpio::Scorpio
MPI::MPI_Fortran
NetCDF::C
```

Compatibility aliases in the modules allow this migration to happen
incrementally.

## Suggested Order

1. Define and apply a local template for private find modules.
2. Modernize simpler modules first: `YAJL`, `Chaparral`, `Scorpio`, and
   `METIS`.
3. Modernize transitive modules: `PETACA`, `Exodus`, `HYPRE`, and `NetCDF`.
4. Modernize `MUMPS` last because it has multiple libraries, optional config
   package behavior, and more link-order sensitivity.
5. Migrate top-level and source call sites to namespaced targets.
6. Keep or later remove unnamespaced aliases depending on how much external
   compatibility the project wants to preserve.

## Validation

Add lightweight configure-time checks for the private modules where practical.
The checks should verify that:

- `find_package(... REQUIRED)` fails when required pieces are missing.
- Version checks are honored.
- Expected imported targets are created.
- Include directories propagate through target usage requirements.
- Static-library transitive dependencies remain available at link time.
