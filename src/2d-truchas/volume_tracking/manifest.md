# 2D volume tracking manifest

The `volume_tracking` directory contains three layers.

## Reusable 2D volume-tracking library

The CMake target `volume_tracking_2d` builds these modules:

- `volume_tracker_2d_class.F90`: abstract interface for initialization, VOF advancement, and inflow-material specification.
- `simple_volume_tracker_type.F90`: first-order donor-cell transport. It is diffusive and does not reconstruct interfaces.
- `geometric_volume_tracker_type.F90`: geometric VOF transport with interface-normal estimation, plane reconstruction, Brent iteration, geometric fluxes, subcycling, flux renormalization, boundedness enforcement, and boundary inflow-material handling.
- `truncation_volume_2d_type.F90`: computes the portion of a triangle or quadrilateral lying behind a plane.
- `locate_plane_os_2d_function.F90`: positions an interface plane using truncation volumes and a Brent solve.
- `plane_2d_type.F90`: basic line/plane geometry.
- `cell_geom_2d_vof_type.F90`: packs local cell geometry for geometric flux calculations.
- `gradient_2d_cc_function.F90`: cell-centered Cartesian and axisymmetric `r-z` gradients.
- `geom_axisymmetric.F90`: converts planar cell volumes and face measures into their `2*pi`-revolved axisymmetric counterparts.

## Legacy test-driver support

The `test_vof_2d` object collection contains:

- `read_inputfile.F90`: a positional-text input reader.
- `gaussian_quadrature_vofinit.F90`: quadrature and element-coordinate transformations for triangle and quadrilateral initialization.
- `vof_2d_test_driver.F90`: a standalone prescribed-velocity timestep driver. It selects the simple or geometric tracker, computes face-normal velocities, converts them into cell-oriented flux velocities, advances the VOF field, and gathers results for comparison.

This is not yet a general flow/VOF coupling layer; it is a standalone advection driver with its own state and time loop.

## Executables and reference data

The CMake target builds five programs:

- `trunc_vol`: low-level truncation-volume checks, including Cartesian and axisymmetric cases.
- `vof_advection`: Cartesian advection of a circular material distribution by a uniform velocity on a regular quadrilateral mesh.
- `vof_vortex`: Cartesian advection in a time-dependent vortex field.
- `vof_axisymmetric`: axisymmetric `r-z` advection of a circular region, supporting quadrilateral and mixed triangle/quad meshes.
- `vof_axialflow`: axisymmetric advection with a prescribed radial/axial velocity field, including a perturbed mesh case.

The `input_*.txt` files are legacy scenario configurations. They specify mesh dimensions, perturbation/probability settings, timestep count, tracker choice, and whether to generate or compare reference data.

The `circle*.txt` files are reference (“gold”) results containing final cell volume fractions in external cell order. The drivers compare computed fields with an `L-infinity` tolerance of `1e-7`; these are regression data, not analytic solutions generated during the test.

## Cartesian and axisymmetric paths

The axisymmetric implementation is a mode of the same tracker classes rather than a separate tracker family. The 2D mesh represents a meridional `r-z` cross-section. With `axisym=.true.`:

- mesh cell volumes and face areas become their `2*pi`-revolved measures;
- gradients use `gradient_rz_cc`;
- truncation volumes use revolved triangle areas;
- geometric flux calculations use axisymmetric face and cell measures.

The same simple or geometric tracker is then used.

Conceptually:

```text
2D tracker algorithm
├── Cartesian geometry
└── axisymmetric r-z geometry
```

The collection is functional but still legacy-oriented: raw text input, standalone drivers, TLS logging, gold files in the source/build workflow, and duplicated driver logic. The reusable geometric tracker is the most valuable part to preserve while the surrounding test harness and interfaces are modernized.
