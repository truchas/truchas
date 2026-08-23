# 2D heat-transport tests

The component tests are the six `test_ht_2d_*.F90` programs in this
directory. They share the material and mesh fixture in `support/` and run in
both serial and four-process configurations.

The JSON-driven integration tests are organized by problem. Most cases own an
`input.json` and a short `test.py` checker. The `linear` case contains
explicitly named quad and triangle x- and y-directed variants.

- `uniform`, `linear` (x- and y-directed cases), `source`, and `nonlinear`
  test single-material conduction and source behavior;
- `multimaterial_uniform`, `multimaterial_source`,
  `multimaterial_mixed_source`, and `multimaterial_mixed_tri_source` test
  material properties and mixed cells;
- `multimaterial_conduction` tests the steady composite-conductivity profile;
- `enthalpy` tests equivalent single-phase specific-heat and
  specific-enthalpy material specifications;
- `htc` tests an external heat-transfer-coefficient boundary condition;
- `radiation` contains transient tests of simple radiation boundary
  conditions;
- `nonlinear_inclusion` is a reference-data regression for nonlinear
  heterogeneous conduction, corresponding to mainline diffusion test DS2;
- `developer_mesh` and `developer_mesh_tri` exercise developer mesh options.

See `linear/README.md` for the x/y analytic cases and `partition/README.md`
for the standalone VTKHDF partition test.

Reference-data cases keep their VTKHDF reference under `reference/`. Shared
execution and VTKHDF helpers are in `support/ht_2d_test_util.py`.
