# 2D heat-transport tests

The component tests are the six `test_ht_2d_*.F90` programs in this
directory. They share the material and mesh fixture in `support/` and run in
both serial and four-process configurations.

The JSON-driven integration tests are organized by problem. Each case owns
an `input.json` and a short `test.py` checker:

- `uniform`, `linear`, `source`, and `nonlinear` test single-material
  conduction and source behavior;
- `multimaterial_uniform`, `multimaterial_source`,
  `multimaterial_mixed_source`, and `multimaterial_mixed_tri_source` test
  material properties and mixed cells;
- `multimaterial_conduction` tests the steady composite-conductivity profile;
- `developer_mesh` and `developer_mesh_tri` exercise developer mesh options.

The `partition_test.py` scripts compare serial and four-process VTKHDF output
in external cell and node order for selected cases. Shared execution and
VTKHDF helpers are in `support/ht_2d_test_util.py`.
